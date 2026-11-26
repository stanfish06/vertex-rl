"""Time alternating dev/new CLI pairs against the frozen 2D rectangle."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import select
import shutil
import statistics
import subprocess
import sys
import time

from compare import CADENCE, DURATION, HERE, FrameMismatch, compare_frozen
from run import read_placement

REPETITIONS = 5
BUILDS = ("dev", "new")


def write_json(path, obj):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def schedule(budgets):
    for threads in budgets:
        for repetition in range(1, REPETITIONS + 1):
            for build in BUILDS if repetition % 2 else reversed(BUILDS):
                yield threads, repetition, build


def metrics(budgets, runs):
    medians = {}
    for threads in budgets:
        for build in BUILDS:
            batch = [r for r in runs if r["threads"] == threads and r["build"] == build]
            if len(batch) == REPETITIONS and all(r["valid"] for r in batch):
                medians[threads, build] = statistics.median(r["launch_to_exit_seconds"] for r in batch)
    results = {}
    for threads in budgets:
        dev, new = (medians.get((threads, build)) for build in BUILDS)
        serial = medians.get((1, "new"))
        results[str(threads)] = {
            "dev_median_seconds": dev, "new_median_seconds": new,
            "dev_over_new": dev / new if dev and new else None,
            "new_serial_over_new": serial / new if serial and new else None,
        }
    return results


def validate(record, directory):
    try:
        record["maximum_errors"] = compare_frozen(directory / "output")
    except ValueError as error:
        if isinstance(error, FrameMismatch):
            record["first_failing_threshold"] = error.threshold
        record["errors"].append(str(error))

    try:
        placements, divisions = read_placement(directory / "run.log")
        with (HERE / "placement.csv").open() as stream:
            reference = list(csv.reader(stream))[1:]
        if placements != reference:
            record["first_failing_threshold"] = 0.0
            record["errors"].append("initial placement or zones differ from the seed-646 reference")
        record.update(accepted_centers=len(placements), divisions=divisions)
        if divisions != 22:
            record["errors"].append(f"solver logged {divisions} divisions, reference 22")
    except ValueError as error:
        record["errors"].append(str(error).replace(str(directory / "run.log"), "run.log.gz"))


def wait_for_exit(process, timeout):
    # Linux pidfds wake at exit. Popen.wait(timeout) can round up by 50 ms.
    try:
        pidfd = os.pidfd_open(process.pid)
    except (AttributeError, OSError):
        return process.wait(timeout=timeout)
    try:
        readable, _, _ = select.select([pidfd], [], [], timeout)
        if not readable:
            raise subprocess.TimeoutExpired(process.args, timeout)
        return process.wait()
    finally:
        os.close(pidfd)


def launch(directory, binary, build, threads, repetition, timeout):
    directory.mkdir(parents=True)
    shutil.copytree(HERE / "initial", directory / "initial")
    (directory / "output").mkdir()
    commands = directory / "commands"
    commands.write_text(f"set number_of_threads {threads}\nset outFreq {CADENCE}\nset tMAX {DURATION}\nrun\n")
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads), OMP_DYNAMIC="FALSE")
    argv = [str(binary), "--reference", "--seed", "646"]
    record = dict(build=build, threads=threads, repetition=repetition, directory=directory.name,
                  valid=False, exit_code=None, launch_to_exit_seconds=None,
                  first_failing_threshold=None, errors=[])
    log = directory / "run.log"
    with commands.open() as stdin, log.open("wb") as stdout:
        start = time.perf_counter()
        try:
            with subprocess.Popen(argv, env=env, cwd=directory, stdin=stdin,
                                  stdout=stdout, stderr=subprocess.STDOUT) as proc:
                try:
                    wait_for_exit(proc, timeout)
                except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                    proc.kill()
                    proc.wait()
                    if isinstance(error, KeyboardInterrupt):
                        record["interrupted"] = True
                        record["errors"].append("benchmark interrupted during this launch")
                    else:
                        record["errors"].append(f"launch timed out after {timeout} seconds")
                record["launch_to_exit_seconds"] = time.perf_counter() - start
                record["exit_code"] = proc.returncode
        except OSError as error:
            record["errors"].append(f"launch failed: {error}")
    if record["exit_code"] not in (None, 0):
        record["errors"].append(f"solver exit code {record['exit_code']}")
    validate(record, directory)
    record["valid"] = not record["errors"]
    write_json(directory / "run.json", record)
    # Keep stdout capture and its flushes identical for both builds inside the timer.
    with log.open("rb") as source, gzip.open(directory / "run.log.gz", "wb", compresslevel=1) as compressed:
        shutil.copyfileobj(source, compressed)
    log.unlink()
    return record


def paired(binaries, destination, full_cpus, timeout):
    budgets = list(dict.fromkeys((1, 2, 4, full_cpus)))
    destination.mkdir(parents=True)
    report = dict(machine=platform.node(), platform=platform.platform(),
                  cpu_affinity=sorted(os.sched_getaffinity(0)), full_cpus=full_cpus,
                  budgets=budgets, seed=646, duration=DURATION, cadence=CADENCE,
                  timeout_seconds=timeout, repetitions=REPETITIONS,
                  binaries={k: {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                            for k, p in binaries.items()},
                  runs=[], complete=False, active_run=None)
    summary = destination / "summary.json"
    for threads, repetition, build in schedule(budgets):
        name = f"threads-{threads}-pair-{repetition}-{build}"
        report["active_run"] = name
        write_json(summary, report)
        print(f"running {name}", flush=True)
        record = launch(destination / name, binaries[build], build, threads, repetition, timeout)
        report["runs"].append(record)
        report["metrics"] = metrics(budgets, report["runs"])
        report["active_run"] = None
        write_json(summary, report)
        print(json.dumps(record), flush=True)
        if record.get("interrupted"):
            return 130
    report["complete"] = True
    write_json(summary, report)
    print(json.dumps(report["metrics"], indent=2), flush=True)
    return 0 if all(r["valid"] for r in report["runs"]) else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("destination", type=Path, help="new results directory")
    parser.add_argument("--full-cpus", type=int, default=len(os.sched_getaffinity(0)),
                        help="CPU count for the full-machine configuration (default: caller's affinity size)")
    parser.add_argument("--timeout", type=float, default=120,
                        help="wall-clock timeout in seconds per process (default: 120)")
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be finite and positive")
    if args.full_cpus < 1 or args.full_cpus > len(os.sched_getaffinity(0)):
        parser.error("--full-cpus must be a positive count within this process's CPU affinity")
    try:
        binaries = {key: path.resolve(strict=True)
                    for key, path in zip(BUILDS, (args.baseline, args.candidate))}
        return paired(binaries, args.destination, args.full_cpus, args.timeout)
    except KeyboardInterrupt:
        print(f"benchmark interrupted; inspect {args.destination}/summary.json", file=sys.stderr)
        return 130
    except OSError as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    sys.exit(main())
