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

## November 26, 2026 results

The workload ran on `lima-devflow-task-dda41e3bf92e`, an Ubuntu 26.04 LTS VM with four AMD Ryzen 7 8845HS guest cores and one hardware thread per core. Both the instrumented dev baseline (`c91a7b1`) and candidate (`53e40b8`) used g++ 15.2.0, C++17, `-O3 -fopenmp -DVERTEX_ARRAY_CAPACITY=1024` and system Readline 8. Timing used the command above with a 120-second timeout and `OMP_MAX_ACTIVE_LEVELS=1`. The fate trace was disabled during timing.

| Threads | Dev passes | Candidate passes | Dev median (s) | Candidate median (s) | Dev / candidate | Candidate scaling |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 5/5 | 5/5 | 9.291 | 1.547 | 6.00× | 1.00× |
| 2 | 0/5 | 5/5 | — | 1.216 | — | 1.27× |
| 4 (full VM) | 0/5 | 5/5 | — | 1.006 | — | 1.54× |

Every candidate frame reproduced the reference sample time, cell IDs, reported areas and vertex positions exactly. All parallel dev runs exceeded the time tolerance. The first failure occurred at threshold 0.5 for two-thread repetition 4 and at threshold 0 for the other nine runs. Each of those processes exited successfully. The missing medians and ratios reflect these trajectory failures. Four threads reduced the candidate's elapsed time by 35% relative to one thread on this rectangle. [Results for all 30 launches](paired-results.json) · [Ordered fate replay](fate-trace.md)
