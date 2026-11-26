"""Failure isolation, launch evidence and timing gates for the paired benchmark."""
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

from compare import FrameMismatch, HERE, compare_frozen
from paired import launch, metrics, schedule


class PairedBenchmark(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def test_earlier_error_in_incomplete_trajectory(self):
        with ZipFile(HERE / "reference.zip") as archive:
            archive.extractall(self.directory)
        index = self.directory / "frame_index.csv"
        with index.open() as f:
            rows = list(csv.reader(f))[:4]
        rows[2][-1] = "0.6"
        with index.open("w") as f:
            csv.writer(f).writerows(rows)
        with self.assertRaises(FrameMismatch) as failure:
            compare_frozen(self.directory)
        self.assertEqual(failure.exception.threshold, 0.5)
        self.assertIn("actual time", str(failure.exception))
        rows[2][-1] = "0.501726508"
        with index.open("w") as f:
            csv.writer(f).writerows(rows)
        with self.assertRaises(FrameMismatch) as failure:
            compare_frozen(self.directory)
        self.assertEqual(failure.exception.threshold, 1.5)

    def test_failed_and_timed_out_processes(self):
        binary = self.directory / "solver"
        binary.write_text("#!/bin/sh\necho 'solver log' >&2\nexit 7\n")
        binary.chmod(0o755)
        failed = self.directory / "failed"
        result = launch(failed, binary, "dev", 2, 1, 5)
        self.assertEqual(result["first_failing_threshold"], 0)
        self.assertEqual(result["exit_code"], 7)
        self.assertFalse(result["valid"])
        self.assertEqual(json.loads((failed / "run.json").read_text()), result)
        with gzip.open(failed / "run.log.gz", "rt") as log:
            self.assertEqual(log.read(), "solver log\n")
        binary.write_text("#!/bin/sh\nexec sleep 10\n")
        result = launch(self.directory / "timed-out", binary, "dev", 4, 1, 0.1)
        self.assertEqual(result["exit_code"], -9)
        self.assertFalse(result["valid"])
        self.assertIn("timed out", result["errors"][0])
        self.assertGreaterEqual(result["launch_to_exit_seconds"], 0.1)

    def test_order_and_median_gates(self):
        planned = list(schedule([1, 2]))
        self.assertEqual(len(planned), 20)
        self.assertEqual([b for _, _, b in planned[:10]], ["dev", "new", "new", "dev"] * 2 + ["dev", "new"])
        runs = [dict(build=b, threads=n, repetition=r, valid=True,
                     launch_to_exit_seconds=(2 if b == "dev" else 1) * r / n) for n, r, b in planned]
        self.assertEqual(metrics([1, 2], runs)["2"], {
            "dev_median_seconds": 3.0, "new_median_seconds": 1.5,
            "dev_over_new": 2.0, "new_serial_over_new": 2.0,
        })
        runs[-2]["valid"] = False
        result = metrics([1, 2], runs)["2"]
        self.assertIsNone(result["dev_median_seconds"])
        self.assertIsNone(result["dev_over_new"])
        self.assertEqual(result["new_serial_over_new"], 2.0)
        result = metrics([1, 2], runs[:-1])["2"]
        self.assertIsNone(result["new_median_seconds"])
        self.assertIsNone(result["new_serial_over_new"])


if __name__ == "__main__":
    unittest.main()
