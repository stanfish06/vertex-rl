"""Gymnasium tissue controls and PettingZoo parallel cell-fate decisions."""
import copy
import math
from pathlib import Path
import tempfile

import gymnasium as gym
from gymnasium import spaces
import numpy as np
from pettingzoo import ParallelEnv

from .experiment import number
from .scoring import TargetScorer, measure_tissue, shape_moments
from .transport import Worker, WorkerProtocolError, _start_workers, positive_timeout


class SolverError(RuntimeError):
    def __init__(self, state):
        self.state, self.time = state, state["time"]
        super().__init__(f"solver error at {self.time}: {state.get('message', 'failed substep')}")


class _Environment:
    def __init__(self, experiment, scorer=None, *, step_timeout=None, readiness_timeout=30,
                 directory=None, worker_id=None, worker_command=None):
        self.experiment, self.scorer = experiment, scorer
        if isinstance(scorer, TargetScorer):
            if scorer.times[0] > 0 or scorer.times[-1] < experiment.time_limit:
                raise ValueError(f"reference [{scorer.times[0]}, {scorer.times[-1]}] must cover "
                                 f"reset 0 and episode limit {experiment.time_limit}")
            if scorer.zone_measures != experiment.zone_measures:
                raise ValueError("target scorer zone measures differ from the experiment")
        self.step_timeout = None if step_timeout is None else positive_timeout(step_timeout)
        self.readiness_timeout = positive_timeout(readiness_timeout)
        if directory is not None:
            Path(directory).mkdir(parents=True, exist_ok=True)
        # Explicit directory retains the private run after close for inspection.
        self._scratch = tempfile.TemporaryDirectory(prefix="vertex-") if directory is None else None
        self.directory = Path(self._scratch.name if self._scratch else tempfile.mkdtemp(dir=directory,prefix="vertex-"))
        command = experiment.stage(self.directory)
        self._worker = Worker(command if worker_command is None else worker_command, self.directory, worker_id)
        self._active = self._closed = False
        self._state = None
        self.time, self._interval = 0., 0

    def _reset(self, seed):
        if self._closed:
            raise RuntimeError("environment was closed")
        if isinstance(seed, np.integer):
            seed = int(seed)
        if type(seed) is not int or not 0 <= seed < 2**32:
            raise ValueError("seed must be an unsigned 32-bit integer")
        self._active = False
        _start_workers([self._worker], self.readiness_timeout)
        state = self._worker.request(self.experiment.reset_payload(seed,self.mode))
        self._accept(state)
        self.seed_used, self._interval = seed, 0
        self._active = state["status"] == "ok"
        return state

    def _accept(self, state):
        status = state["status"]
        if status == "error":
            raise ValueError(state.get("message","worker rejected request"))
        if status == "solver_error":
            self._state, self.time = state, state["time"]
            raise SolverError(state)
        if status not in ("ok","terminated","truncated") or "cells" not in state or "time" not in state:
            self._worker.close()
            raise WorkerProtocolError(f"unexpected state reply {state}")
        self._state, self.time = state, state["time"]

    def _require_active(self):
        if not self._active:
            raise RuntimeError("explicit reset required before stepping")

    def _advance(self, **action):
        self._require_active()
        self._active = False
        limit = self.experiment.time_limit
        scheduled = (self._interval+1)*self.experiment.decision_interval
        # A rounded product such as 3 * 0.3 is one ULP short of 0.9.
        deadline = limit if scheduled >= limit-math.ulp(limit) else scheduled
        state = self._worker.request({"op":"step","deadline":deadline,**action}, self.step_timeout)
        if state["status"]=="error":
            self._active = True  # The worker rejected the request before applying an action.
        self._accept(state)
        self._interval += 1
        measurements = self._measure()
        reward = 0. if self.scorer is None else number(self.scorer(measurements,self.time),"tissue reward")
        self._active = state["status"]=="ok"
        return measurements,reward,state["status"]=="terminated",state["status"]=="truncated"

    def _measure(self):
        return measure_tissue(self._state["cells"],self.experiment.zone_measures,self.experiment.dimension)

    def _info(self):
        status = self._state["status"]
        return {"time":self.time,"reason":{"truncated":"time_limit","terminated":"empty_tissue"}.get(status)}

    def close(self):
        self._active, self._closed = False, True
        self._worker.close()
        if self._scratch:
            self._scratch.cleanup()


class TissueEnv(_Environment, gym.Env):
    """One absolute action vector, in the experiment's declared control order."""
    mode = "tissue"
    metadata = {"render_modes":[]}
    render_mode = None

    def __init__(self, experiment, scorer=None, **kwargs):
        super().__init__(experiment, scorer, **kwargs)
        if not experiment.controls:
            raise ValueError("tissue mode requires at least one control")
        self.action_space = spaces.Box(np.array([c.lower for c in experiment.controls],dtype=np.float64),
                                       np.array([c.upper for c in experiment.controls],dtype=np.float64), dtype=np.float64)
        d = experiment.dimension
        size = 5 + len(experiment.zone_measures)*(4+2*d*d)
        self.observation_space = spaces.Box(-np.inf,np.inf,(size,),dtype=np.float64)
        self.measurements = None
        self._actions = []

    def _observation(self, measurements):
        self.measurements = measurements
        values = [self.time, measurements["count"], *(measurements["size_percentiles"] or [0.]*3)]
        dim = self.experiment.dimension
        for zone in measurements["zones"]:
            values.extend((zone["count"],zone["density"]))
            for key in ("orientation","arrangement"):
                tensor = zone[key]
                values.append(float(tensor is not None))
                values.extend(np.ravel(tensor) if tensor is not None else [0.]*(dim*dim))
        return np.asarray(values,dtype=np.float64)

    def reset(self, *, seed=None, options=None):
        gym.Env.reset(self, seed=seed)
        self._reset(seed if seed is not None else int(self.np_random.integers(2**32)))
        self._actions = []
        return self._observation(self._measure()), self._info()

    def _vector(self, action):
        try:
            vector = np.asarray(action)
        except ValueError as error:
            raise ValueError(f"expected action shape {self.action_space.shape}, received ragged array") from error
        if vector.shape != self.action_space.shape:
            raise ValueError(f"expected action shape {self.action_space.shape}, received {vector.shape}")
        values = []
        for c,v in zip(self.experiment.controls,vector):
            if not isinstance(v,(float,int,np.integer,np.floating)) or isinstance(v,(bool,np.bool_)) or not c.lower<=v<=c.upper:
                raise ValueError(f"{c.identifier} requested {v!r}, allowed [{c.lower}, {c.upper}]")
            values.append(float(v))
        return values

    def step(self, action):
        self._require_active()
        values = self._vector(action)
        controls = [{"zone":c.zone,"parameter":c.parameter,"value":v} for c,v in zip(self.experiment.controls,values)]
        measurements,reward,terminated,truncated = self._advance(controls=controls)
        self._actions.append(values)
        return self._observation(measurements), reward, terminated, truncated, self._info()

    def recording(self):
        """A JSON-serializable copy of accepted actions and their episode settings."""
        if self._state is None:
            raise RuntimeError("reset required before recording")
        return {"seed":self.seed_used,"decision_interval":self.experiment.decision_interval,
                "time_limit":self.experiment.time_limit,"controls":[c.identifier for c in self.experiment.controls],
                "actions":copy.deepcopy(self._actions)}

    def replay(self, recording):
        """Preflight the complete recording, reset its seed, then run its absolute actions."""
        for setting in ("decision_interval","time_limit"):
            actual, expected = recording.get(setting), getattr(self.experiment,setting)
            if actual!=expected:
                raise ValueError(f"{setting}: recorded {actual!r}, configured {expected!r}")
        actual = recording["controls"]
        expected = [c.identifier for c in self.experiment.controls]
        for i in range(max(len(actual),len(expected))):
            a = actual[i] if i<len(actual) else "<missing>"
            b = expected[i] if i<len(expected) else "<missing>"
            if a!=b:
                raise ValueError(f"control position {i}: recorded {a}, configured {b}")
        validated = []
        for i, vector in enumerate(recording["actions"],start=1):
            try:
                validated.append(self._vector(vector))
            except ValueError as error:
                raise ValueError(f"recorded interval {i}: {error}") from error
        initial = self.reset(seed=recording["seed"])
        return initial, [self.step(vector) for vector in validated]


FATES = ("grow","stop_growing","divide","die")


class CellEnv(_Environment, ParallelEnv):
    """Each observed lifetime ID takes one fate action at the next boundary."""
    mode = "cell"
    metadata = {"name":"vertex_cells_v0","render_modes":[],"is_parallelizable":True}
    render_mode = None

    def __init__(self, experiment, scorer=None, **kwargs):
        super().__init__(experiment,scorer,**kwargs)
        self.agents = []
        d = experiment.dimension
        self._cell_space = spaces.Box(-np.inf,np.inf,(6+d,),dtype=np.float64)
        bound = np.pi if d==2 else 1.
        self._fate_space = spaces.Dict({"fate":spaces.Discrete(4),
                                        "orientation":spaces.Box(-bound,bound,(1 if d==2 else 3,),dtype=np.float64)})
        self._rng = np.random.default_rng()
        self._zones = {z["id"]:z for z in experiment.zones}

    def action_space(self, agent):
        return self._fate_space

    def observation_space(self, agent):
        return self._cell_space

    def _local(self):
        observations, infos = {},{}
        dim = self.experiment.dimension
        for c in self._state["cells"]:
            if c.get("dead"):
                continue
            agent = c["id"]
            shape = shape_moments(c)
            observations[agent] = np.asarray([shape.size,c["target_area" if dim==2 else "target_volume"],
                shape.elongation,*shape.axis,c["age"],float(c["growing"]),c["zone"]],dtype=np.float64)
            zone = self._zones[c["zone"]]
            legacy = dim==2 and "rule" not in zone
            mask = [zone.get("allow_growth",True),zone.get("allow_stop",True),
                    zone.get("allow_division",not legacy or c["zone"]==1),
                    zone.get("allow_death",not legacy or c["zone"]==3)]
            infos[agent] = {**self._info(),"allowed_fates":tuple(mask)}
        return observations,infos

    def reset(self, seed=None, options=None):
        if seed is not None:
            self._rng = np.random.default_rng(seed)
        self._reset(seed if seed is not None else int(self._rng.integers(2**32)))
        obs,info = self._local()
        self.agents = list(obs) if self._active else []
        return obs,info

    def _fates(self, actions):
        if not isinstance(actions,dict) or set(actions)!=set(self.agents):
            raise ValueError(f"expected one fate per live agent {self.agents}")
        fates = []
        for agent in self.agents:
            action = actions[agent]
            if isinstance(action,str):
                action = {"fate":action}
            fate = action.get("fate") if isinstance(action,dict) else None
            if isinstance(fate,(np.integer,int)) and not isinstance(fate,(bool,np.bool_)) and 0<=fate<len(FATES):
                fate = FATES[fate]
            if fate not in FATES:
                raise ValueError(f"agent {agent}: invalid fate {fate!r}")
            row = {"id":agent,"action":fate}
            if fate=="divide":
                orientation = np.asarray(action.get("orientation"),dtype=float)
                if self.experiment.dimension==2:
                    if orientation.shape not in ((),(1,)) or not np.isfinite(orientation).all():
                        raise ValueError(f"agent {agent}: expected a finite division angle")
                    row["orientation"] = float(orientation.item())
                else:
                    if orientation.shape!=(3,) or not np.isfinite(orientation).all() or not np.any(orientation):
                        raise ValueError(f"agent {agent}: expected a nonzero finite plane normal")
                    row["orientation"] = orientation.tolist()
            fates.append(row)
        return fates

    def step(self,actions):
        self._require_active()
        actors = list(self.agents)
        _, reward, ended, limit = self._advance(actions=self._fates(actions))
        obs, info = self._local()
        participants = set(actors) | obs.keys()
        retired = set(actors) - obs.keys()
        for agent in retired:
            info[agent] = {**self._info(),"retired":True}
        terminations = {a:ended or a in retired for a in participants}
        truncations = {a:limit and a not in retired for a in participants}
        self.agents = list(obs) if self._active else []
        return obs, dict.fromkeys(actors,reward), terminations, truncations, info
