"""Exercise the native 3D CLI and persistent worker on actual simulations."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
CUBE = ROOT / "examples/3d/cube.json"


def box(lo, hi):
    return {"shape": "box", "min": lo, "max": hi}


class Solver3DTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build = tempfile.TemporaryDirectory()
        cls.addClassCleanup(build.cleanup)
        cls.binary = Path(build.name) / "vertex3d"
        sources = ["main", "driver", "config", "export", "solver", "mesh", "contact"]
        subprocess.run([os.environ.get("CXX", "g++"), "-std=c++17", "-O2", "-fopenmp",
                        *shlex.split(os.environ.get("CPPFLAGS", "")),
                        *(str(ROOT / "model/three_d" / (name + ".cpp")) for name in sources),
                        *shlex.split(os.environ.get("LDFLAGS", "")),
                        *shlex.split(os.environ.get("JSON_LIBS", "-ljson-c")), "-o", str(cls.binary)], check=True)

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.cube = json.loads(CUBE.read_text())

    def replies(self, config, *requests, threads=1):
        messages = [dict(op="ping"), dict(op="reset", seed=646, time_limit=1, experiment=config), *requests]
        result = subprocess.run([str(self.binary), "worker", "--threads", str(threads)], cwd=self.directory,
                                input="".join(json.dumps(x) + "\n" for x in messages), capture_output=True,
                                text=True, timeout=30, check=True)
        replies = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertEqual(len(replies), len(messages), result.stderr)
        self.assertEqual(replies[0], {"status": "ready"})
        return replies[1:]

    def test_native_training_modes(self):
        reset=dict(op="reset",mode="cell",seed=646,time_limit=.1,experiment=self.cube,
                   zones=[{"id":10,"allow_death":False},{"id":20,"allow_death":True}])
        step=lambda t,actions: dict(op="step",deadline=t,actions=actions)
        stop={"id":1,"action":"stop_growing"}
        divide={"id":1,"action":"divide","orientation":[1,0,0]}
        _, start, invalid, stopped, cut, died = self.replies(self.cube,reset,step(.001,[stop,stop]),
            step(.001,[stop]),step(.002,[divide]),
            step(.003,[{"id":2,"action":"die"},{"id":3,"action":"die"}]))
        self.assertEqual(invalid["status"],"error")
        self.assertFalse(stopped["cells"][0]["growing"])
        self.assertEqual([c["id"] for c in cut["cells"]],[2,3])
        self.assertTrue(all(not c["growing"] and c["age"]==.001 for c in cut["cells"]))
        self.assertAlmostEqual(sum(c["target_volume"] for c in cut["cells"]),stopped["cells"][0]["target_volume"])
        self.assertEqual([c["id"] for c in died["cells"]],[2])
        control=lambda name,value: {"parameter":name,"zone":10,"value":value}
        settings=[control("growth_rate",1),control("target_size",2)]
        reset.update(mode="tissue")
        _, initial, invalid, untouched, first, second = self.replies(self.cube,reset,
            dict(op="step",deadline=.02,controls=[*settings,control("pressure",3)]),
            dict(op="measure"),dict(op="step",deadline=.02,controls=settings),
            dict(op="step",deadline=.04,controls=settings))
        self.assertEqual(invalid["status"],"error")
        self.assertEqual(untouched,initial)
        self.assertAlmostEqual(first["cells"][0]["target_volume"],1.02)
        self.assertAlmostEqual(second["cells"][0]["target_volume"],1.04)

    def test_cube_collection(self):
        output = self.directory / "frames"
        subprocess.run([self.binary, "run", CUBE, output, "--threads", "2"], check=True, capture_output=True)
        entries = ET.parse(output / "tissue.pvd").findall(".//DataSet")
        self.assertEqual([float(node.get("timestep")) for node in entries], [0, 0.05, 0.1])
        piece = ET.parse(output / entries[1].get("file")).find(".//Piece")
        arrays = {a.get("Name"): a.text.split() for a in piece.findall("CellData/DataArray")}
        ids = list(map(int, arrays["cell_id"]))
        self.assertEqual(set(ids), {2, 3})
        self.assertEqual(len(ids), int(piece.get("NumberOfPolys")))
        for cell, zone, volume in zip(ids, map(int, arrays["zone_id"]), map(float, arrays["volume"])):
            self.assertEqual(zone, 10 if cell == 2 else 20)
            self.assertAlmostEqual(volume, 0.5, places=12)

    def test_cut_growth_and_zone_order(self):
        self.cube["zones"][0].update(allow_growth=False, growth_rate=0.2)
        self.cube["zones"][1].update(growth_rate=0.2)
        self.cube["cells"][0]["target_volume"] = 1.2
        self.cube["mechanics"]["pressure"] = 0
        self.cube["mechanics"]["volume_stiffness"] = 0
        initial, birth, after = self.replies(self.cube, dict(op="divide", cell_id=1, normal=[1,0,0]),
                                            dict(op="step", deadline=0.0123456789))
        self.assertEqual(initial["cells"][0]["zone"], 10)
        self.assertAlmostEqual(initial["cells"][0]["volume"], 1)
        self.assertAlmostEqual(initial["cells"][0]["surface_area"], 6)
        self.assertTrue(birth["division_accepted"])
        self.assertEqual([c["zone"] for c in birth["cells"]], [10, 20])
        self.assertTrue(all(c["growing"] for c in birth["cells"]))
        self.assertAlmostEqual(sum(c["target_volume"] for c in birth["cells"]), 1.2)
        self.assertEqual(after["time"], 0.0123456789)
        for c, delta in zip(after["cells"], [0, 0.2 * after["time"]]):
            self.assertAlmostEqual(c["target_volume"], 0.6 + delta)
        self.cube["zones"].reverse()
        self.cube["zones"][0]["allow_division"] = False
        reset, rejected = self.replies(self.cube, dict(op="divide", cell_id=1, normal=[1,0,0]))
        self.assertEqual(reset["cells"][0]["zone"], 20)
        self.assertFalse(rejected["division_accepted"])
        self.assertEqual(reset["cells"], rejected["cells"])

    def test_contacts_thread_replay(self):
        self.cube["cells"] = [{"mesh": box([x-.5,-.5,-.5], [x+.5,.5,.5])} for x in (-1.04,0,1.04)]
        reset, serial = self.replies(self.cube, dict(op="step", deadline=0.02))
        self.assertEqual(serial["status"], "ok")
        self.assertLess(serial["cells"][0]["center"][0], reset["cells"][0]["center"][0])
        self.assertGreater(serial["cells"][2]["center"][0], reset["cells"][2]["center"][0])
        for n in (2,4):
            self.assertEqual(self.replies(self.cube, dict(op="step", deadline=0.02), threads=n), [reset, serial])

    def test_wall_retry_rollback_and_explicit_reset(self):
        self.cube["domain"] = box([-.52]*3, [.52]*3)
        self.cube["cells"][0]["target_volume"] = 3
        self.cube["mechanics"].update(volume_stiffness=8, surface_tension=0, pressure=0,
                                      max_step=.3, min_step=.001, max_displacement=5)
        reset, failure, measured, denied, new = self.replies(
            self.cube, dict(op="step", deadline=.2), dict(op="measure"), dict(op="step", deadline=.3),
            dict(op="reset", experiment=self.cube, seed=646, time_limit=1))
        self.assertEqual(failure["status"], "solver_error")
        self.assertGreater(failure["time"], 0)
        self.assertLess(failure["time"], .2)
        self.assertGreater(failure["rejected_trials"], 0)
        self.assertEqual(failure["cells"], measured["cells"])
        self.assertEqual(measured["time"], failure["time"])
        self.assertEqual(denied["status"], "error")
        self.assertEqual(new, reset)
        _, replay = self.replies(self.cube, dict(op="step", deadline=failure["time"]))
        self.assertEqual(replay["status"], "ok")
        self.assertEqual(failure["cells"], replay["cells"])
        for c in failure["cells"]:
            self.assertTrue(all(abs(coord)<.52 for p in c["vertices"] for coord in p))

        # Export the generated wall's exact triangulation as an OBJ, then use
        # it as the fixed wall for the same crossing trial.
        shape_scene = dict(domain=box([-2]*3,[2]*3), cells=[dict(mesh=self.cube["domain"])])
        shape = self.replies(shape_scene)[0]["cells"][0]
        obj = self.directory / "wall.obj"
        with obj.open("w") as stream:
            for p in shape["vertices"]:
                print("v", *p, file=stream)
            for tri in shape["faces"]:
                print("f", *(i+1 for i in tri), file=stream)
        imported = dict(self.cube, domain={"shape":"import", "path":"wall.obj"})
        self.assertEqual(self.replies(imported, dict(op="step",deadline=.2)), [reset,failure])

    def test_curved_shapes_and_seeded_division(self):
        for wall in (dict(shape="sphere", center=[0]*3, radius=2, max_surface_error=.2),
                     dict(shape="cylinder", center=[0]*3, radius=2, height=3, max_surface_error=.2),
                     dict(shape="ellipsoid", center=[0]*3, radii=[2,2.5,3], max_surface_error=.2)):
            self.cube["domain"] = wall
            reset = self.replies(self.cube)[0]
            self.assertEqual(reset["status"], "ok", reset)
        self.cube["zones"] = [dict(id=4, rule=dict(kind="radial", center=[0]*3, upper=3),
                                   division_volume=.8, division_rate=1e6)]
        reset, split, again, replay = self.replies(
            self.cube, dict(op="step",deadline=.01), dict(op="reset", experiment=self.cube, seed=646, time_limit=1),
            dict(op="step",deadline=.01))
        self.assertEqual(len(split["cells"]), 2)
        self.assertEqual(split,replay)
        self.assertEqual(reset,again)

    def test_imported_zone_interface(self):
        (self.directory / "plane.obj").write_text(
            "v 0 -1 -1\nv 0 1 -1\nv 0 1 1\nv 0 -1 1\nf 1 2 3\nf 1 3 4\n")
        rule = dict(kind="interface", path="plane.obj")
        self.cube["zones"] = [dict(id=8,rule=dict(rule,positive_side=False)),
                              dict(id=3,rule=dict(rule,positive_side=True))]
        initial, split = self.replies(self.cube, dict(op="divide",cell_id=1,normal=[1,0,0]))
        self.assertEqual(initial["cells"][0]["zone"], 8)
        self.assertEqual([c["zone"] for c in split["cells"]], [8,3])

    def test_invalid_resets_preserve_episode(self):
        passing, *rest = self.replies(
            self.cube,
            dict(op="reset",seed=1,time_limit=1,experiment={**self.cube, "cells":self.cube["cells"]*2}),
            dict(op="reset",seed=1,time_limit=1,experiment={**self.cube, "domain":box([-.2]*3,[.2]*3)}),
            dict(op="measure"))
        self.assertEqual([r["status"] for r in rest], ["error", "error", "ok"])
        self.assertEqual(passing,rest[-1])


if __name__ == "__main__":
    unittest.main()
