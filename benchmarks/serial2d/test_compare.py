"""Calibrate the comparator on modified copies of the recorded dev frames."""
import csv
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

from compare import AREA_ABSOLUTE, AREA_RELATIVE, HERE, POSITION, TIME, compare_frozen


class RectangleComparison(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.output = Path(directory.name)
        with ZipFile(HERE / "reference.zip") as archive:
            archive.extractall(self.output)

    def edit(self, name, change):
        file = self.output / name
        with file.open() as stream:
            reader = csv.DictReader(stream)
            rows = list(reader)
            columns = reader.fieldnames
        change(rows)
        with file.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    def test_reference(self):
        self.assertTrue(all(value == 0 for value in compare_frozen(self.output).values()))

    def test_forward_rotation_and_reordering_cells(self):
        def rotate(rows):
            ring = [r for r in rows if r["cell_id"] == "1"]
            rows[:] = [r for r in rows if r["cell_id"] != "1"] + ring[3:] + ring[:3]
        self.edit("snapshot_1.txt", rotate)
        self.assertEqual(compare_frozen(self.output)["vertex_displacement"], 0)

    def test_winding(self):
        def reverse(rows):
            rows[:] = list(reversed(rows))
        self.edit("snapshot_3.txt", reverse)
        with self.assertRaisesRegex(ValueError, "threshold 1.0: cell 1: .*reversed winding"):
            compare_frozen(self.output)

    def test_time_before_geometry(self):
        def shift(rows):
            rows[5]["actual_time"] = float(rows[5]["actual_time"]) + 4 * TIME
        self.edit("frame_index.csv", shift)
        (self.output / "snapshot_6.txt").unlink()
        with self.assertRaisesRegex(ValueError, "threshold 2.5: actual time"):
            compare_frozen(self.output)

    def test_position(self):
        def displace(rows):
            # Both coordinates are within the limit, but the Euclidean distance exceeds it.
            for axis in ("v_x", "v_y"):
                rows[0][axis] = float(rows[0][axis]) + 0.9 * POSITION
        self.edit("snapshot_4.txt", displace)
        with self.assertRaisesRegex(ValueError, "threshold 1.5: cell 1: vertex displacement"):
            compare_frozen(self.output)

    def check_area(self, frame, cell, fraction):
        def change(rows):
            for row in rows:
                if row["cell_id"] == str(cell):
                    area = float(row["A"])
                    row["A"] = area + fraction * (AREA_ABSOLUTE + AREA_RELATIVE * abs(area))
        self.edit(f"snapshot_{frame}.txt", change)
        if fraction > 1:
            with self.assertRaisesRegex(ValueError, f"cell {cell}: area"):
                compare_frozen(self.output)
        else:
            self.assertGreater(compare_frozen(self.output)["area"], 0)

    def test_physical_changes(self):
        def move_vertex(rows):
            rows[0]["v_x"] = float(rows[0]["v_x"]) + 0.001
        self.edit("snapshot_3.txt", move_vertex)
        with self.assertRaisesRegex(ValueError, "vertex displacement"):
            compare_frozen(self.output)
        for frame, cell in [(2, 17), (13, 3)]:
            with self.subTest(frame=frame):
                with ZipFile(HERE / "reference.zip") as archive:
                    archive.extractall(self.output)
                def increase_area(rows):
                    for row in rows:
                        if row["cell_id"] == str(cell):
                            row["A"] = float(row["A"]) * 1.01
                self.edit(f"snapshot_{frame}.txt", increase_area)
                with self.assertRaisesRegex(ValueError, f"cell {cell}: area"):
                    compare_frozen(self.output)

    def test_small_area(self):
        self.check_area(2, 17, 2)

    def test_large_area(self):
        self.check_area(13, 3, 2)

    def test_within_area_limit(self):
        self.check_area(13, 3, 0.99)

    def test_nonfinite(self):
        self.edit("snapshot_2.txt", lambda rows: rows[0].update(vs_y="NaN"))
        with self.assertRaisesRegex(ValueError, "threshold 0.5: cell 1: nonfinite"):
            compare_frozen(self.output)

    def test_missing_cell(self):
        def remove_cell(rows):
            rows[:] = [row for row in rows if row["cell_id"] != "1"]
        self.edit("snapshot_1.txt", remove_cell)
        with self.assertRaisesRegex(ValueError, r"differing IDs \[1\]"):
            compare_frozen(self.output)

    def test_missing_vertex(self):
        self.edit("snapshot_4.txt", lambda rows: rows.pop(0))
        with self.assertRaisesRegex(ValueError, "threshold 1.5: cell 1: vertex count"):
            compare_frozen(self.output)

    def test_missing_threshold(self):
        self.edit("frame_index.csv", lambda rows: rows.pop(3))
        with self.assertRaisesRegex(ValueError, r"missing or unexpected thresholds: \[1.5\]"):
            compare_frozen(self.output)

    def test_duplicate_threshold(self):
        self.edit("frame_index.csv", lambda rows: rows[3].update(scheduled_time="1"))
        with self.assertRaisesRegex(ValueError, "threshold 1.0: invalid or duplicate"):
            compare_frozen(self.output)


if __name__ == "__main__":
    unittest.main()
