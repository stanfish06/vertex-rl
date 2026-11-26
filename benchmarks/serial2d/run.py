"""Launch one fresh single-thread rectangle simulation in a new directory."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from compare import CADENCE, DURATION, HERE, THRESHOLDS, compare_frozen, read_index, read_snapshot


def read_placement(log):
    placements = {}
    zones = {}
    attempts = 0
    divisions = 0
    with log.open() as stream:
        for line in stream:
            if "ERROR:" in line:
                raise ValueError(f"solver reported an error; see {log}")
            if line.startswith("Point "):
                attempts += 1
            elif line.startswith("placement,"):
                _, attempt, cell, x, y = line.strip().split(",")
                placements[int(cell)] = [attempt, cell, x, y]
            elif line.startswith("initial_zone,"):
                _, cell, zone = line.strip().split(",")
                zones[int(cell)] = zone
            divisions += line.count("start divisionProliferative zone division")
    if attempts != 20 or placements.keys() != zones.keys():
        raise ValueError("incomplete fresh placement log")
    return [row + [zones[cell]] for cell, row in placements.items()], divisions


def run(binary, destination, seed):
    binary = binary.resolve(strict=True)
    destination.mkdir(parents=True)
    shutil.copytree(HERE / "initial", destination / "initial")
    output = destination / "output"
    output.mkdir()
    argv = [str(binary), "--reference"]
    if seed is not None:
        argv.extend(["--seed", str(seed)])
    commands = f"set number_of_threads 1\nset outFreq {CADENCE}\nset tMAX {DURATION}\nrun\n"
    environment = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_THREAD_LIMIT="1")
    log = destination / "run.log"
    with log.open("w") as stream:
        start = time.perf_counter()
        subprocess.run(argv, input=commands, text=True, env=environment, cwd=destination,
                       stdout=stream, stderr=subprocess.STDOUT, check=True)
        elapsed = time.perf_counter() - start
    placements, divisions = read_placement(log)
    placement_path = destination / "placement.csv"
    with placement_path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["attempt", "cell_id", "x", "y", "zone"])
        writer.writerows(placements)
    effective_seed = 646 if seed is None else seed
    errors = None
    if effective_seed == 646:
        if placement_path.read_bytes() != (HERE / "placement.csv").read_bytes():
            raise ValueError("seed 646 initial placement or zones differ from reference")
        errors = compare_frozen(output)
    counts = []
    index = read_index(output)
    for threshold in THRESHOLDS:
        frame, _ = index[threshold]
        counts.append(len(read_snapshot(output, frame)))
    if counts[0] != len(placements) or counts[-1] - counts[0] != divisions:
        raise ValueError("division log and sampled cell counts disagree")
    report = dict(seed=effective_seed, explicit_seed=seed is not None, threads=1,
                  duration=DURATION, cadence=CADENCE, launch_to_exit_seconds=elapsed,
                  accepted_centers=len(placements), initial_zone_counts={
                      str(i): sum(row[-1] == str(i) for row in placements) for i in range(4)
                  }, divisions=divisions, sampled_cell_counts=counts)
    if errors is not None:
        report["maximum_errors"] = errors
    print(json.dumps(report, indent=2))
    (destination / "run.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("destination", type=Path, help="a new run directory")
    parser.add_argument("--seed", type=int, help="omit to exercise the built-in default of 646")
    args = parser.parse_args()
    if args.seed is not None and not 0 <= args.seed <= 2**32 - 1:
        parser.error("seed must be an unsigned 32-bit integer")
    try:
        run(args.binary, args.destination, args.seed)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{args.destination}: {error}\n")
