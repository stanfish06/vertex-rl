# Python training environments

Use Python 3.10 or later on Linux. From the repository root, install the package and build the workers with a C++17 compiler, OpenMP and json-c development headers. The capacity below fits the small 2D fixtures.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e .
make -C model vertex3d libvertex2d.so CPPFLAGS=-DVERTEX_ARRAY_CAPACITY=1024
```

`TissueEnv` implements the Gymnasium environment API. This experiment samples a 28 × 20 rectangle with seed 646, sets two tissue controls at each decision, and records the episode. Each environment launches its worker on first reset and reuses the process for subsequent episodes.

```python
import json
from pathlib import Path
from vertex_rl import Control, Experiment, TissueEnv

experiment = Experiment(
    dimension=2,
    geometry={"domain": {"shape": "rectangle", "min": [0, 0], "max": [28, 20]}},
    zones=[{"id": 0, "rule": {"kind": "all"}, "allow_division": False}],
    decision_interval=0.01,
    time_limit=0.025,
    controls=[Control(0, "growth_rate", 0, 2), Control(0, "target_size", 0.2, 3)],
    binary="model/libvertex2d.so",
)
env = TissueEnv(experiment, scorer=lambda measurements, time: -abs(measurements["count"] - 18))
try:
    observation, info = env.reset(seed=646)
    while True:
        observation, reward, terminated, truncated, info = env.step([0.4, 1.2])
        if terminated or truncated:
            break
    recording = env.recording()
    Path("actions.json").write_text(json.dumps(recording))
    initial, transitions = env.replay(recording)
finally:
    env.close()
```

The absolute action vector follows the declared control order and uses simulation units. Growth rate increases target area or volume per time unit until the zone's target size. Division and death rates are hazards per time unit. The default controls select growth rate and target size for each zone, in zone order. Their preset ranges are [0, 1] and [0.01, 10]. Custom bounds include both endpoints. A rejected vector leaves the worker at the same boundary. Replay checks decision interval, episode limit, ordered identifiers, and every vector against current bounds before resetting and dispatching the recording.

The tissue observation is a float64 vector. It starts with time, live cell count and the 10th, 50th and 90th size percentiles. Each configured zone then contributes count, density, an orientation-tensor availability flag and its row-major matrix, followed by the equivalent flag and matrix for arrangement. Missing percentiles and matrices use zeros. `env.measurements` holds the structured scoring input. `info["time"]` is the actual simulated end time. A terminal reason is `empty_tissue` or `time_limit`.

## Geometry and zones

A 2D domain accepts a rectangle with `min` and `max`, a circle with `center` and `radius`, an ellipse with `center` and `radii`, or a polygon with `points` in boundary order. Circles and ellipses default to 64 facets. Set `segments` to change the resolution. A polygon can also specify a `path` to a headerless two-column CSV. The legacy five-file `initial/` directory is accepted directly as `geometry`.

3D accepts a standalone experiment JSON path or the corresponding Python dictionary. Use `binary="model/vertex3d"`. Relative geometry file references in JSON resolve from the experiment file's parent. In dictionaries they resolve from the working directory where the experiment is constructed. See the [3D input format](../model/three_d/README.md).

Declare zone rules in order. `layer` uses `origin`, `axis` and optional signed-distance bounds `lower` and `upper`. `radial` uses `center` and optional radius bounds. An oriented 2D `interface` uses a polyline's `points` and `positive_side` (left of its oriented segments). A closed 2D `surface` uses a polygon and `inside`. Both accept a CSV `path` in place of points. An `all` rule covers the remaining domain. Zone settings can include `growth_rate`, `target_size`, `division_rate`, `death_rate`, `allow_growth`, `allow_stop`, `allow_division` and `allow_death`.

Pass `zone_measures={zone_id: area_or_volume, ...}` for density measurements. A single whole-domain 2D polygon or 3D box supplies its measure automatically. For partitioned domains, a curved or imported 3D wall, or legacy inputs, supply the physical measures of the configured zones. The run uses one native solver thread by default. Set `threads` in the experiment to choose a different budget.

## A shared cell policy

`CellEnv(experiment, scorer=...)` implements the PettingZoo parallel API. Reset returns observation and information dictionaries keyed by episode-local lifetime ID. Call `step(actions)` with exactly one entry per ID in `env.agents`. It returns the observation, reward, termination, truncation and information dictionaries.

```python
from vertex_rl import CellEnv

cells = CellEnv(experiment, scorer=lambda measurements, time: -abs(measurements["count"] - 18))
try:
    observations, infos = cells.reset(seed=646)
    actions = {agent: "grow" for agent in cells.agents}
    observations, rewards, terminated, truncated, infos = cells.step(actions)
finally:
    cells.close()
```

A local observation contains size, target size, elongation, unit long-axis components (zero for an isotropic cell), time since division, growing flag and zone ID. The per-agent information includes `allowed_fates` in action order. For Gym spaces, the action is a dictionary with `fate` in 0 through 3 (grow, stop growing, divide, die) and `orientation`. The orientation space is a one-element radian angle in 2D or a three-component normal in 3D. Division needs a nonzero normal in 3D. Plain fate names and dictionaries such as `{"fate":"divide", "orientation":[1,0,0]}` are also accepted. The environment sends orientation only for a division.

Each acting ID receives the same interval reward, including a retired parent or a cell that died. Daughters receive observations and first act at the next boundary. Termination flags cover retired agents, and truncation flags cover surviving agents at the time limit. Episode-local IDs generated by division make the future agent set dynamic.

## Scoring and worker lifetime

A scorer receives `(measurements, actual_time)` and returns a finite scalar. Omitting it gives a zero reward for unscored data collection. For target matching, pass `TargetScorer(samples, time_limit, zone_measures)`. Each sample has a `time` and `measurements` in the same format as `TissueEnv.measurements`. The scorer sorts whole rows, checks duplicates and coverage, and selects the latest target at or before each interval end. Reward is the negative mean error across cell count, size distribution, zone density, and the average of orientation and arrangement.

Each environment gets its own private working directory. Pass `directory=parent_directory` to retain its generated inputs and worker stderr after close. `start_batch(envs, timeout=30)` launches and pings a new group under one shared monotonic readiness deadline. It returns after all workers report ready, before any reset. A failed startup reaps the batch. `ReadinessTimeout.pending_ids` identifies missing workers. `WorkerExitError` carries the worker ID, exit status and stderr.

Steps wait for a complete JSON line with no wall-clock timeout by default. Set `step_timeout=seconds` on either environment to bound each exchange, including partial writes and reads. A timeout raises `StepTimeout` and reaps the process. An exit before the reply raises `WorkerExitError`. Explicit reset starts a replacement in either case. `SolverError.state` and `SolverError.time` give the last committed state following a numerical failure. Reset reuses that responsive worker. The failed action produces no reward, recording entry or automatic replay.
