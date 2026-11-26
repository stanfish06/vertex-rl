"""Public training APIs, atomic replay preflight and real worker pipe faults."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

try:
    import numpy as np
    from vertex_rl import (CellEnv, Control, Experiment, ReadinessTimeout, SolverError, StepTimeout,
                           TissueEnv, WorkerExitError, start_batch)
except ModuleNotFoundError as error:
    if error.name not in ("numpy","gymnasium","pettingzoo"):
        raise
    raise unittest.SkipTest(f"install vertex-rl to run training environment tests ({error.name})")

FIXTURE = [sys.executable, str(Path(__file__).with_name("training_worker_fixture.py").resolve())]


class EnvironmentTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def env(self, cls=TissueEnv, exp_args=None, **kwargs):
        options = dict(dimension=2, geometry={"shape":"rectangle","min":[0,0],"max":[10,10]},
                       time_limit=.25, decision_interval=.1,
                       controls=[Control(0,"growth_rate",0,2),Control(0,"target_size",1,5)])
        options.update(exp_args or {})
        env = cls(Experiment(**options), worker_command=FIXTURE, directory=self.tmp.name, **kwargs)
        self.addCleanup(env.close)
        return env

    @staticmethod
    def requests(env):
        return [json.loads(s) for s in (env.directory/"requests.jsonl").read_text().splitlines()]

    def test_invalid_action_and_replay_preflight(self):
        env = self.env()
        env.reset(seed=3)
        for action,error in [([[1,2]],"expected action shape"),([2,6],"zone_0.target_size"),
                             ([float("nan"),2],"zone_0.growth_rate")]:
            with self.assertRaisesRegex(ValueError,error):
                env.step(action)
        self.assertEqual([r["op"] for r in self.requests(env)],["ping","reset"])
        for action in ([0,3],[2,5],[0,3]):
            last=env.step(action)
        self.assertEqual(last[2:],(False,True,{"time":.25,"reason":"time_limit"}))
        saved=json.loads(json.dumps(env.recording()))
        for field,value,error in [("decision_interval",.2,"recorded 0.2, configured 0.1"),
                                  ("time_limit",.3,"recorded 0.3, configured 0.25"),
                                  ("controls",saved["controls"][::-1],"control position 0"),
                                  ("controls",saved["controls"][:-1],"recorded <missing>"),
                                  ("actions",[[0,3],[2,6]],"recorded interval 2.*target_size")]:
            with self.assertRaisesRegex(ValueError,error):
                env.replay({**saved,field:value})
        self.assertEqual(len(self.requests(env)),5)
        other = self.env(exp_args={"controls":[Control(0,"growth_rate",0,2),Control(0,"target_size",3,5)]})
        initial, replay = other.replay(saved)
        self.assertEqual(initial[1]["time"],0.)
        np.testing.assert_array_equal(last[0],replay[-1][0])
        self.assertEqual(last[1:],replay[-1][1:])
        self.assertEqual([x["deadline"] for x in self.requests(other) if x["op"]=="step"],[.1,.2,.25])

    def test_errors_and_recovery_in_both_modes(self):
        for cls in (TissueEnv,CellEnv):
            for fault,exception in (("partial",StepTimeout),("exit",WorkerExitError),("solver_error",SolverError)):
                with self.subTest(mode=cls.__name__,fault=fault):
                    calls=[]
                    env=self.env(cls,step_timeout=.08,scorer=lambda m,t: calls.append(t) or -1.)
                    start=env.reset(seed=646)
                    process=env._worker.process
                    action=[1.,2.] if cls is TissueEnv else {1:"grow"}
                    (env.directory/"step.fault").write_text(fault)
                    with self.assertRaises(exception) as failure:
                        env.step(action)
                    if exception is SolverError:
                        self.assertEqual(failure.exception.time,.05)
                        self.assertEqual(failure.exception.state["cells"][0]["id"],1)
                        self.assertIsNone(process.poll())
                    else:
                        self.assertIsNotNone(process.poll())
                    if exception is WorkerExitError:
                        self.assertEqual(failure.exception.status,7)
                        self.assertIn("exit diagnostic",failure.exception.stderr)
                    self.assertEqual(calls,[])
                    with self.assertRaisesRegex(RuntimeError,"explicit reset"):
                        env.step(action)
                    renewed=env.reset(seed=646)
                    self.assertEqual(process is env._worker.process, exception is SolverError)
                    fresh=self.env(cls,scorer=lambda m,t:-1.)
                    independent=fresh.reset(seed=646)
                    np.testing.assert_equal(start,renewed)
                    np.testing.assert_equal(independent,renewed)
                    for _ in range(3):
                        np.testing.assert_equal(env.step(action),fresh.step(action))
                    self.assertEqual(calls,[.1,.2,.25])
                    self.assertEqual(len([r for r in self.requests(env) if r["op"]=="step"]),4)

    def test_delayed_reply_without_step_timeout(self):
        for cls in (TissueEnv,CellEnv):
            env=self.env(cls)
            env.reset(seed=646)
            (env.directory/"step.fault").write_text("delay")
            env.step([1.,2.] if cls is TissueEnv else {1:"grow"})
            self.assertEqual(env.time,.1)

    def test_rounded_final_deadline(self):
        env=self.env(exp_args={"decision_interval":.3,"time_limit":.9})
        env.reset(seed=646)
        replies=[env.step([0,1]) for _ in range(3)]
        self.assertEqual([reply[4]["time"] for reply in replies],[.3,.6,.9])
        self.assertTrue(replies[-1][3])

    def test_readiness_batch_deadline_and_exit(self):
        for fault,exception in (("stall",ReadinessTimeout),("exit",WorkerExitError)):
            envs=[self.env(worker_id=f"child-{i}") for i in range(3)]
            for i in (0,2):
                (envs[i].directory/"ping.fault").write_text(fault if i==2 else "stall")
            with self.assertRaises(exception) as failure:
                start_batch(envs,timeout=.4)
            if fault=="stall":
                self.assertEqual(failure.exception.pending_ids,["child-0","child-2"])
            else:
                self.assertEqual(failure.exception.worker_id,"child-2")
                self.assertEqual(failure.exception.status,7)
                self.assertIn("exit diagnostic",failure.exception.stderr)
            self.assertTrue(all(e._worker.process is None for e in envs))
            start_batch(envs,timeout=2)
            self.assertTrue(all(e._worker.ready for e in envs))

    def test_python_2d_setup(self):
        folder=Path(self.tmp.name)
        curve=folder/"curve.csv"
        curve.write_text("0,0\n10,1\n")
        for shape,fields in (("ellipse",{"center":[5,8],"radii":[12,9]}),
                             ("circle",{"center":[0,0],"radius":10}),
                             ("polygon",{"points":[[0,0],[10,0],[8,8],[0,10]]})):
            zones=[{"id":7,"rule":{"kind":"interface","path":str(curve),"positive_side":True}},
                   {"id":8,"rule":{"kind":"all"}}]
            env=self.env(exp_args={"geometry":{"domain":{"shape":shape,**fields},"zones":zones},
                                   "zone_measures":{7:20,8:25},"controls":[Control(7,"growth_rate",0.,2.)]})
            request=env.experiment.reset_payload(646,"cell")
            self.assertEqual(request["zones"][0]["rule"]["points"],[[0.,0.],[10.,1.]])
            self.assertEqual(len(list((env.directory/"initial").glob("*.txt"))),5)
        with self.assertRaisesRegex(ValueError,"intersects|zero area"):
            self.env(exp_args={"geometry":{"shape":"polygon","points":[[0,0],[1,1],[0,1],[1,0]]}})


if __name__=="__main__":
    unittest.main()
