"""Integration tests for the compiled 2D solver and its persistent worker."""
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "benchmarks" / "serial2d" / "initial"


class WorkerClient:
    def __init__(self, directory, library, threads):
        self.errors = tempfile.TemporaryFile(mode="w+")
        self.process = subprocess.Popen(
            [sys.executable, str(ROOT / "model/worker.py"), str(library), "--threads", str(threads)],
            cwd=directory, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=self.errors, text=True,
        )

    def call(self, op, **fields):
        self.process.stdin.write(json.dumps(dict(op=op, **fields)) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            self.errors.seek(0)
            raise RuntimeError(f"worker exited: {self.errors.read()}")
        return json.loads(line)

    def close(self):
        if self.process.poll() is None:
            self.call("close")
        self.process.stdin.close()
        self.process.stdout.close()
        self.process.wait(timeout=10)
        self.errors.close()


class WorkerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(directory.cleanup)
        cls.library = Path(directory.name) / "libvertex2d.so"
        subprocess.run([os.environ.get("CXX", "g++"), "-std=c++17", "-O3", "-fopenmp",
                        "-fPIC", "-shared", "-DVERTEX_ARRAY_CAPACITY=1024",
                        str(ROOT / "model/worker_api.cpp"), "-o", str(cls.library)], check=True)

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        shutil.copytree(FIXTURE, self.directory / "initial")
        self.worker = self.start()

    def start(self, threads=1):
        client = WorkerClient(self.directory, self.library, threads)
        self.addCleanup(client.close)
        self.assertEqual(client.call("ping"), {"status": "ready"})
        return client

    def test_training_cut_and_lifetime(self):
        reset = dict(seed=646, time_limit=1, mode="cell",
                     zones=[{"id":17, "rule":{"kind":"all"}, "allow_death":True, "growth_rate":0}])
        state = self.worker.call("reset", **reset)
        before = state["cells"]
        ids = [cell["id"] for cell in before]
        invalid = [{"id":i,"action":"die"} for i in ids[:-1]]
        self.assertEqual(self.worker.call("step",deadline=0.01,actions=invalid)["status"],"error")
        self.assertEqual(self.worker.call("measure"),state)
        grow = [{"id":i,"action":"grow"} for i in ids]
        grown = self.worker.call("step", deadline=0.001, actions=grow)
        self.assertTrue(all(c["growing"] for c in grown["cells"]))
        actions = [{"id":ids[0], "action":"divide", "orientation":0.345},
                   *[{"id":i,"action":"stop_growing"} for i in ids[1:]]]
        split = self.worker.call("step", deadline=0.001001, actions=actions)
        self.assertEqual(split["status"],"ok")
        daughters = [c for c in split["cells"] if c["id"]>max(ids)]
        self.assertEqual({c["id"] for c in daughters},{max(ids)+1,max(ids)+2})
        self.assertTrue(all(c["zone"]==17 and c["growing"] for c in daughters))
        self.assertTrue(all(abs(c["age"]-0.000001)<1e-13 for c in daughters))
        self.assertAlmostEqual(sum(c["target_area"] for c in daughters),grown["cells"][0]["target_area"],places=6)
        self.assertAlmostEqual(daughters[0]["target_area"]/daughters[1]["target_area"],
                               daughters[0]["area"]/daughters[1]["area"],places=3)
        self.assertEqual(self.worker.call("step", deadline=0.01, actions=actions)["status"],"error")
        deaths = [{"id":c["id"],"action":"die"} for c in split["cells"]]
        ended = self.worker.call("step", deadline=0.01, actions=deaths)
        self.assertEqual((ended["status"],ended["time"],ended["cells"]),("terminated",split["time"],[]))
        self.assertEqual(self.worker.call("reset",**reset),state)

    def test_training_control_atomicity_and_permissions(self):
        z = {"id":27,"rule":{"kind":"all"},"allow_division":False,"allow_death":False}
        reset = dict(seed=646,time_limit=1,zones=[z])
        start = self.worker.call("reset",mode="tissue",**reset)
        def control(name,value):
            return {"zone":27,"parameter":name,"value":value}
        bad = [control("growth_rate",0),control("kA",1)]
        self.assertEqual(self.worker.call("step",deadline=0.001,controls=bad)["status"],"error")
        self.assertEqual(self.worker.call("measure"),start)
        repeated = [control("growth_rate",2),control("target_size",4)]
        a=self.worker.call("step",deadline=0.001,controls=repeated)
        b=self.worker.call("step",deadline=0.002,controls=repeated)
        self.assertGreater(b["cells"][0]["target_area"],a["cells"][0]["target_area"])
        self.assertGreater(a["cells"][0]["target_area"],start["cells"][0]["target_area"])
        start=self.worker.call("reset",mode="cell",**reset)
        actions=[{"id":c["id"],"action":"divide","orientation":.1} for c in start["cells"]]
        unchanged=self.worker.call("step",deadline=0.001,actions=actions)
        self.assertEqual([c["id"] for c in unchanged["cells"]],[c["id"] for c in start["cells"]])
        self.assertFalse(any(c["growing"] for c in unchanged["cells"]))

    def test_early_divisions_in_translated_ellipse(self):
        # Three sub-threshold cuts produce daughters narrower than the contact shell.
        # Far edges across a disjoint neighbor must not generate penetration impulses.
        for cx, cy in ((0,0),(50,40)):
            with self.subTest(center=(cx,cy)):
                points=[(cx+24*math.cos(k*math.tau/64),cy+17*math.sin(k*math.tau/64)) for k in range(64)]
                sides=list(zip(points,points[1:]+points[:1]))
                tangents=[((b[0]-a[0])/math.dist(a,b),(b[1]-a[1])/math.dist(a,b)) for a,b in sides]
                for name,values in (("domain_coords",points),("domain_v_ts",tangents),
                                    ("domain_v_ns",[(-y,x) for x,y in tangents])):
                    (self.directory/"initial"/(name+".txt")).write_text("".join(f"{x},{y}\n" for x,y in values))
                state=self.worker.call("reset",seed=2,time_limit=.02,mode="cell",
                                       zones=[{"id":0,"rule":{"kind":"all"},"growth_rate":0.1}])
                angles={7:-.594327894,8:1.538350895,13:-.277104725}
                actions=[]
                for c in state["cells"]:
                    i=c["id"]
                    if i in angles:
                        actions.append({"id":i,"action":"divide","orientation":angles[i]})
                    else:
                        actions.append({"id":i,"action":"die" if i in (9,11,12) else "grow"})
                final=self.worker.call("step",deadline=.02,actions=actions)
                self.assertEqual(final["status"],"truncated",final.get("message"))
                self.assertFalse(set(angles)&{c["id"] for c in final["cells"]})
                self.assertTrue(all(c["area"]>0.04 for c in final["cells"]))

    def test_deadline_motion_and_export(self):
        reset = self.worker.call("reset", seed=646, time_limit=0.01)
        self.assertEqual(reset["measurements"]["cell_count"], 17)
        self.assertEqual(sum(reset["measurements"]["zone_counts"]), 17)
        for cell in reset["cells"]:
            self.assertGreater(cell["area"], 0)
        first = self.worker.call("step", deadline=0.002)
        self.assertEqual((first["status"], first["time"]), ("ok", 0.002))
        self.assertEqual(self.worker.call("measure"), first)
        final = self.worker.call("step", deadline=0.025)
        self.assertEqual((final["status"], final["time"]), ("truncated", 0.01))
        self.assertEqual(self.worker.call("step", deadline=0.03)["status"], "error")

        # Both intervals are shorter than the first natural step (~0.0106).
        # Doubling the requested duration doubles the actual first displacement.
        self.worker.call("reset", seed=646, time_limit=0.01)
        second = self.worker.call("step", deadline=0.004)
        p = reset["cells"][0]["vertices"][0]
        a = first["cells"][0]["vertices"][0]
        b = second["cells"][0]["vertices"][0]
        self.assertGreater(max(abs(u - v) for u, v in zip(a, p)), 1e-6)
        for u, v, origin in zip(a, b, p):
            self.assertAlmostEqual(v - origin, 2 * (u - origin), delta=2e-6)
        self.assertEqual(self.worker.call("export", frame_id=3)["time"], 0.004)
        frame = self.directory / "output/snapshot_3.txt"
        with frame.open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), sum(len(cell["vertices"]) for cell in second["cells"]))
        self.assertAlmostEqual(float(rows[0]["v_x"]), b[0], places=12)

    def test_reset_fresh_and_invalid_requests(self):
        settings = dict(seed=646, time_limit=0.2)
        original = self.worker.call("reset", **settings)
        expected = self.worker.call("step", deadline=0.012345678901)
        self.assertEqual(expected["time"], 0.012345678901)
        self.assertEqual(self.worker.call("step", deadline=-1)["status"], "error")
        self.assertEqual(self.worker.call("reset", **settings, parameters={"bad": 1})["status"], "error")
        self.assertEqual(self.worker.call("measure"), expected)
        self.worker.call("reset", **settings, parameters={"kA": 0.1})
        altered = self.worker.call("step", deadline=expected["time"])
        self.assertNotEqual(altered["cells"], expected["cells"])
        self.assertEqual(self.worker.call("reset", **settings), original)
        self.assertEqual(self.worker.call("step", deadline=expected["time"]), expected)
        fresh = self.start()
        self.assertEqual(fresh.call("reset", **settings), original)
        self.assertEqual(fresh.call("step", deadline=expected["time"]), expected)
        other = fresh.call("reset", seed=647, time_limit=0.2)
        self.assertNotEqual(other["cells"], original["cells"])

    def test_failure_retains_last_accepted_substep(self):
        # Closely placed cells in this smaller rectangle drive the legacy contact
        # rule unstable. Exercise a production failure after accepted substeps.
        initial = self.directory / "initial"
        (initial / "domain_coords.txt").write_text("0,0\n16,0\n16,12\n0,12\n")
        (initial / "interface_coords.txt").write_text("0,9,0,6,0,3\n16,9,16,6,16,3\n")
        settings = dict(seed=646, time_limit=2)
        original = self.worker.call("reset", **settings)
        first = self.worker.call("step", deadline=0.005)
        failure = self.worker.call("step", deadline=1)
        self.assertEqual(failure["status"], "solver_error")
        self.assertGreater(failure["time"], first["time"])
        self.assertLess(failure["time"], 1)
        self.assertEqual(self.worker.call("measure")["cells"], failure["cells"])
        self.assertEqual(self.worker.call("export", frame_id=1)["time"], failure["time"])
        self.assertEqual(self.worker.call("step", deadline=1)["status"], "error")

        fresh = self.start()
        self.assertEqual(fresh.call("reset", **settings), original)
        self.assertEqual(fresh.call("step", deadline=0.005), first)
        accepted = fresh.call("step", deadline=failure["time"])
        self.assertEqual(accepted["status"], "ok")
        for key in ("time", "cells", "measurements"):
            self.assertEqual(failure[key], accepted[key])
        # Setting CV explicitly also updates the corresponding contact radius.
        # The compiled dev defaults for those radii are independent of CV.
        explicit = dict(resting_CV=1, proliferative_CV=1.5, hypertrophic_CV=2.5, ossification_CV=2.5)
        self.worker.call("reset", **settings, parameters=explicit)
        self.assertNotEqual(self.worker.call("step", deadline=0.005)["cells"], first["cells"])
        self.assertEqual(self.worker.call("reset", **settings), original)
        self.assertEqual(self.worker.call("step", deadline=0.005), first)

    def test_thread_count_replays_through_division_and_reset(self):
        settings = dict(seed=646, time_limit=2)
        def trajectory(client):
            frames = [client.call("reset", **settings)]
            for deadline in (0.25, 1.5, 2):
                frames.append(client.call("step", deadline=deadline))
            return frames
        reference = trajectory(self.worker)
        self.assertEqual(reference[-1]["status"], "truncated")
        self.assertGreater(len(reference[-1]["cells"]), len(reference[0]["cells"]))
        for threads in (2, 4):
            with self.subTest(threads=threads):
                worker = self.start(threads)
                self.assertEqual(trajectory(worker), reference)
                self.assertEqual(trajectory(worker), reference)

    def test_empty_tissue(self):
        (self.directory / "initial/domain_coords.txt").write_text("0,0\n1,0\n1,1\n0,1\n")
        reset = self.worker.call("reset", seed=646, time_limit=1)
        self.assertEqual((reset["status"], reset["time"], reset["cells"]), ("terminated", 0, []))


if __name__ == "__main__":
    unittest.main()
