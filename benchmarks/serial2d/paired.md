# Paired 2D standalone benchmark

Use the instrumented dev CLI from commit `c91a7b1` and the candidate standalone CLI. Build both with the same compiler, optimization flags, OpenMP runtime and array capacity. The [serial fixture](README.md) gives the build flags and frozen tolerances.

From the repository root on Linux, provide the two executables and a new directory on the same machine.

```sh
python3 benchmarks/serial2d/paired.py /tmp/dev-cli /tmp/new-cli /tmp/paired-results --timeout 120
```

The harness uses thread counts 1, 2, 4 and the size of its CPU affinity mask, counting coincident budgets once. On an SMT machine, pass `--full-cpus N` to identify the desired physical-core budget. Both builds receive the same OpenMP environment and CLI thread setting at each budget. Each uses seed 646, duration 6, cadence 0.5 and a fresh copy of the rectangle. The timeout is a positive wall-clock duration in seconds for each launch. Choose it before the comparison.

Five pairs run sequentially at each thread count. Dev starts pairs 1, 3 and 5. The candidate starts pairs 2 and 4. The external monotonic timer covers process launch through exit, including placement, integration, legacy writes, frame-index writes and shutdown. Input preparation, frame comparison and log compression take place outside the timer.

Each run has a directory with its inputs, stdin commands, geometry, frame index, combined solver log in `run.log.gz` and `run.json`. The JSON contains its build, thread count, repetition, elapsed seconds, exit code, validation errors and first failing scheduled threshold when known. A timeout retains the available frames. Even if the index ends early, the comparator checks earlier indexed frames for the first difference. A failure after the last valid frame can have a null failing threshold.

The top-level `summary.json` records the hostname, CPU affinity, binary paths and SHA-256 digests, workload, launch order and individual results. The harness replaces this file after each launch and records the active run before launching it. A stopped comparison leaves `complete` false. An interrupt during a solver run kills and reaps that process and validates its available output.

Each median is the third of the five sorted launch-to-exit times for one build at one thread count. A missing or failed run gives that build a null median at that budget. The same-budget ratio divides the dev median by the candidate median. The candidate scaling ratio divides its single-thread median by its median at the indicated thread count. A ratio appears only when both medians exist. Exit status is 0 when all runs qualify, 1 for a completed comparison with failures or a harness error, and 130 on keyboard interruption.
