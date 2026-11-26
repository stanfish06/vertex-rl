"""Compare a solver output directory with the frozen serial rectangle trajectory."""
import argparse
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from zipfile import Path as ZipPath, ZipFile

HERE = Path(__file__).resolve().parent
DURATION = 6
CADENCE = 0.5
THRESHOLDS = [i * CADENCE for i in range(int(DURATION / CADENCE) + 1)]
TIME = 1e-6
POSITION = 1e-5
AREA_ABSOLUTE = 5e-6
AREA_RELATIVE = 1e-5


def finite(text):
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"nonfinite value {text}")
    return value


def read_index(root):
    index = {}
    ids = set()
    with (root / "frame_index.csv").open("r") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["frame_id", "scheduled_time", "actual_time"]:
            raise ValueError("invalid frame index header")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError("invalid frame index column count")
            t = finite(row["scheduled_time"])
            frame = int(row["frame_id"])
            if frame <= 0 or frame in ids or t in index:
                raise ValueError(f"threshold {t}: invalid or duplicate frame {frame}")
            ids.add(frame)
            index[t] = frame, finite(row["actual_time"])
    if set(index) != set(THRESHOLDS):
        raise ValueError(f"missing or unexpected thresholds: {sorted(set(index) ^ set(THRESHOLDS))}")
    return index


@dataclass
class Cell:
    area: float
    ring: list = field(default_factory=list)


def read_snapshot(root, frame):
    cells = {}
    columns = "cell_id,v_x,v_y,vs_x,vs_y,c_x,c_y,A,A0,dead".split(",")
    with (root / f"snapshot_{frame}.txt").open("r") as stream:
        reader = csv.reader(stream)
        if next(reader, []) != columns:
            raise ValueError("invalid snapshot header")
        for row in reader:
            if len(row) != len(columns):
                raise ValueError("invalid snapshot column count")
            identity = int(row[0])
            try:
                x, y, _, _, _, _, area, _, dead = map(finite, row[1:])
            except ValueError as error:
                raise ValueError(f"cell {identity}: {error}") from error
            if identity < 1 or dead not in (0, 1):
                raise ValueError(f"invalid cell ID or dead flag for cell {identity}")
            if identity not in cells:
                cells[identity] = Cell(area)
            if cells[identity].area != area:
                raise ValueError(f"cell {identity}: inconsistent reported area")
            cells[identity].ring.append((x, y))
    return cells


def signed_double_area(ring):
    return sum(x * v - y * u for (x, y), (u, v) in zip(ring, ring[1:] + ring[:1]))


def compare_cell(ref, other):
    if len(ref.ring) != len(other.ring):
        raise ValueError(f"vertex count {len(other.ring)}, reference {len(ref.ring)}")
    area_error = abs(ref.area - other.area)
    area_limit = AREA_ABSOLUTE + AREA_RELATIVE * abs(ref.area)
    if area_error > area_limit:
        raise ValueError(f"area {other.area}, reference {ref.area}, allowed difference {area_limit}")
    if len(ref.ring) < 3 or signed_double_area(ref.ring) * signed_double_area(other.ring) <= 0:
        raise ValueError("degenerate polygon or reversed winding")
    # One forward cyclic shift must align the whole ring.
    position_error = min(
        max(math.dist(p, q) for p, q in zip(ref.ring, other.ring[k:] + other.ring[:k]))
        for k in range(len(ref.ring))
    )
    if position_error > POSITION:
        raise ValueError(f"vertex displacement {position_error}, limit {POSITION}")
    return area_error, position_error


def compare(reference, candidate):
    ref_index, index = read_index(reference), read_index(candidate)
    maxima = {"actual_time": 0., "area": 0., "vertex_displacement": 0.}
    for t in THRESHOLDS:
        ref_frame, ref_actual = ref_index[t]
        frame, actual = index[t]
        drift = abs(actual - ref_actual)
        if drift > TIME:
            raise ValueError(f"threshold {t}: actual time {actual}, reference {ref_actual}, limit {TIME}")
        maxima["actual_time"] = max(drift, maxima["actual_time"])
        try:
            expected = read_snapshot(reference, ref_frame)
            observed = read_snapshot(candidate, frame)
            if expected.keys() != observed.keys():
                raise ValueError(f"cell count {len(observed)}, reference {len(expected)}, "
                                 f"differing IDs {sorted(expected.keys() ^ observed.keys())}")
            for identity, cell in expected.items():
                try:
                    area, displacement = compare_cell(cell, observed[identity])
                except ValueError as error:
                    raise ValueError(f"cell {identity}: {error}") from error
                maxima["area"] = max(maxima["area"], area)
                maxima["vertex_displacement"] = max(maxima["vertex_displacement"], displacement)
        except (OSError, KeyError, ValueError) as error:
            raise ValueError(f"threshold {t}: {error}") from error
    return maxima


def compare_frozen(output):
    with ZipFile(HERE / "reference.zip") as archive:
        return compare(ZipPath(archive), output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(compare_frozen(args.output), indent=2))
    except (ValueError, OSError) as error:
        parser.exit(1, f"{args.output}: {error}\n")
