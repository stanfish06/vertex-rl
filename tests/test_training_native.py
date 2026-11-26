"""Exercise training transactions and episode recovery through an actual 3D solver."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

try:
    import numpy as np
    from vertex_rl import CellEnv, Control, Experiment, SolverError, TissueEnv
except ModuleNotFoundError as error:
    if error.name not in ("numpy","gymnasium","pettingzoo"):
        raise
    raise unittest.SkipTest(f"install vertex-rl to run training integration tests ({error.name})")

ROOT = Path(__file__).resolve().parents[1]


class NativeTraining(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(directory.cleanup)
        cls.binary = Path(directory.name)/"vertex3d"
        names = ["main","driver","config","export","solver","mesh","contact"]
        subprocess.run([os.environ.get("CXX","g++"),"-std=c++17","-O2","-fopenmp",
                        *shlex.split(os.environ.get("CPPFLAGS","")),
                        *(ROOT/"model"/"three_d"/(n+".cpp") for n in names),
                        *shlex.split(os.environ.get("LDFLAGS","")),
                        *shlex.split(os.environ.get("JSON_LIBS","-ljson-c")),"-o",cls.binary],check=True)

    def config(self):
        cfg=json.loads((ROOT/"examples"/"3d"/"cube.json").read_text())
        cfg["zones"]=[{"id":0,"rule":{"kind":"all"},"allow_death":True}]
        return cfg

    def test_real_solver_failure_and_same_process_reset(self):
        cfg=self.config()
        cfg["domain"]={"shape":"box","min":[-.52]*3,"max":[.52]*3}
        cfg["cells"][0]["target_volume"]=3
        cfg["mechanics"].update(volume_stiffness=8,surface_tension=0,pressure=0,
                                 max_step=.3,min_step=.001,max_displacement=5)
        exp=Experiment(dimension=3,geometry=cfg,decision_interval=.2,time_limit=.4,binary=self.binary,
                       controls=[Control(0,"death_rate",0,1e9)])
        for kind in (TissueEnv,CellEnv):
            with self.subTest(kind=kind.__name__):
                calls=[]
                env=kind(exp,scorer=lambda m,t:calls.append(t) or -.5)
                independent=kind(exp,scorer=lambda m,t:-.5)
                self.addCleanup(env.close)
                self.addCleanup(independent.close)
                initial=env.reset(seed=37)
                process=env._worker.process
                with self.assertRaises(SolverError) as failure:
                    env.step([0] if kind is TissueEnv else {1:"grow"})
                self.assertGreater(failure.exception.time,0)
                self.assertLess(failure.exception.time,.2)
                self.assertEqual(calls,[])
                np.testing.assert_equal(env.reset(seed=37),initial)
                self.assertIs(process,env._worker.process)
                np.testing.assert_equal(independent.reset(seed=37),initial)
                action=[1e9] if kind is TissueEnv else {1:"die"}
                transition=env.step(action)
                np.testing.assert_equal(independent.step(action),transition)
                self.assertEqual(len(calls),1)
                if kind is TissueEnv:
                    self.assertEqual(transition[2:4],(True,False))
                else:
                    self.assertEqual(transition[1:4],({1:-.5},{1:True},{1:False}))

    def test_birth_and_shared_reward(self):
        cfg=self.config()
        cfg["mechanics"].update(volume_stiffness=0,surface_tension=0,pressure=0)
        exp=Experiment(dimension=3,geometry=cfg,time_limit=.15,decision_interval=.1,binary=self.binary)
        env=CellEnv(exp,scorer=lambda tissue,t: -float(tissue["count"]))
        self.addCleanup(env.close)
        initial,_=env.reset(seed=42)
        self.assertEqual(set(initial),{1})
        obs,rewards,terminated,truncated,info=env.step({1:{"fate":"divide","orientation":[1,0,0]}})
        self.assertEqual(set(obs),{2,3})
        self.assertEqual(rewards,{1:-2.})
        self.assertEqual(terminated,{1:True,2:False,3:False})
        self.assertFalse(any(truncated.values()))
        self.assertTrue(info[1]["retired"])
        for local in obs.values():
            self.assertEqual(local.shape,(9,))
            self.assertAlmostEqual(local[0],.5)
            self.assertAlmostEqual(local[1],.5)
        obs,rewards,terminated,truncated,_=env.step({2:"die",3:"stop_growing"})
        self.assertEqual(set(obs),{3})
        self.assertEqual(rewards,{2:-1.,3:-1.})
        self.assertEqual(terminated,{2:True,3:False})
        self.assertEqual(truncated,{2:False,3:True})
        self.assertEqual(env.agents,[])


if __name__=="__main__":
    unittest.main()
