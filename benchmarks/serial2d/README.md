# Serial 2D reference

This rectangle exercises fresh placement, adaptive integration, cell division and legacy output on the dev solver. Use it to check a candidate's trajectory before measuring speed.

## Reproduce the run

From the repository root, use a C++17 compiler with OpenMP, GNU Readline development files and Python 3.9 or later. Give each launch a new directory.

```sh
g++ -std=c++17 -O3 -fopenmp -DVERTEX_ARRAY_CAPACITY=1024 model/main.c -lreadline -o /tmp/vertex-reference
python3 benchmarks/serial2d/run.py /tmp/vertex-reference /tmp/vertex-default
python3 benchmarks/serial2d/run.py /tmp/vertex-reference /tmp/vertex-seed646 --seed 646
python3 benchmarks/serial2d/compare.py /tmp/vertex-seed646/output
python3 -m unittest discover -s benchmarks/serial2d -v
```

The driver copies the five inputs, creates the output directory, sets all solver thread counts to one and supplies the CLI commands through stdin. Both seed-646 modes verify the accepted placements and compare every scheduled frame with the checked-in snapshots. Each run retains its inputs, legacy viewer files, frame index, solver log, placement CSV and summary. The timer in the summary covers process launch through exit, including initialization, writes and shutdown. Input preparation and comparison occur outside the timer.

## Workload and build

The rectangle covers `[0,28] × [0,20]`. Three interfaces at y = 15, 10 and 5 have normals pointing toward increasing y. The first row of the interface coordinate file contains their anchors. The second row gives their other endpoints for the archived readers. All mechanical and fate parameters use dev defaults. The run stops when the clock reaches 6 and samples at intervals of 0.5. The threshold-zero sample captures the first accepted motion.

The build starts at `dev` commit `d141ee1ec236bef1f494cb715565a9e05544ef04`. See commit `1bf8ad26d269ddc894da453ec673ee4a0c48519e` for the isolated repair of the shoulder output streams. On this libc the original writer's two closed stream pointers produced transposed X/Y shoulder files. The repair writes through the two open streams. The separate seed change gives `model/main.c` an unsigned `--seed` argument with a default of 646. In `model/_initial_structure.h`, the original `srand(646)` becomes `srand(placement_seed)` at the same site before the 20 candidate draws. Subsequent fate draws continue from that stream. The optional `--reference` flag records placements, zones and the frame index during a fresh run.

The qualification machine was `lima-devflow-task-dda41e3bf92e`, an x86-64 Ubuntu 26.04 LTS VM on an AMD Ryzen 7 8845HS. The compiler was Ubuntu g++ 15.2.0-16ubuntu1 with glibc 2.43. The build used the flags above, Readline 8.3 development headers with an explicit include path, and `-l:libreadline.so.8`. Dev eagerly initializes three square arrays at a default capacity of 50,000, requiring more than 30 GB of memory. This fixture sets capacity to 1,024 at compile time. An independent run with 2,048 produced identical sampled times, cell IDs, counts, areas and vertices. A capacity overflow terminates the solver with an error. For later comparisons, compile both builds with the same capacity.

## Recorded trajectory

The 646 run accepted 17 of 20 candidates. The initial counts were 4 resting, 7 proliferative, 2 hypertrophic and 4 ossification cells. [placement.csv](placement.csv) records the accepted draw number, ID, center and zone before motion. The solver logged 22 divisions. The final frame contains 39 polygons, including cells with a death flag. All 13 snapshots and their index are in [reference.zip](reference.zip).

| Frame | Scheduled time | Actual time | Cells |
| ---: | ---: | ---: | ---: |
| 1 | 0 | 0.0106460173 | 17 |
| 2 | 0.5 | 0.501726508 | 17 |
| 3 | 1 | 1.00039756 | 17 |
| 4 | 1.5 | 1.50000024 | 19 |
| 5 | 2 | 2.00021887 | 22 |
| 6 | 2.5 | 2.50016713 | 24 |
| 7 | 3 | 3.00008392 | 24 |
| 8 | 3.5 | 3.50023437 | 26 |
| 9 | 4 | 4.00021029 | 31 |
| 10 | 4.5 | 4.50016785 | 32 |
| 11 | 5 | 5.00009298 | 36 |
| 12 | 5.5 | 5.50020933 | 36 |
| 13 | 6 | 6.00008535 | 39 |

A default-seed launch and two independent explicit-646 launches agreed exactly at every sample. Each launch's 104 legacy files matched a separately built dev executable byte for byte with only the capacity adjustment and shoulder-stream repair. The archived MATLAB live script and app's `read_single_time` functions both read all 13 frames in GNU Octave 11.1, using `containers.Map` for MATLAB's `dictionary`. The recovered vertex counts, polygon and shoulder coordinates, centroids and death flags matched the snapshots. This check exercised reader functions in Octave. Native MATLAB and Wolfram GUI rendering remain unverified.

## Frozen comparison limits

| Quantity | Limit | Unit |
| --- | ---: | --- |
| Actual time at the same threshold | 0.000001 | simulation time |
| Euclidean displacement per vertex | 0.00001 | simulation length |
| Absolute area allowance | 0.000005 | simulation area |
| Relative area allowance | 0.00001 | dimensionless |

Across repeated serial runs, the largest time difference, vertex displacement and area difference were all zero. The time allowance is 0.0002% of one sampling interval. Adjacent vertex spacings ranged from 0.01046715 to 1.585207, so the displacement allowance is below 0.1% of the shortest sampled edge. Reported areas ranged from 0.26605988 (frame 2, cell 17) to 3.46811676 (frame 13, cell 3). The area limit is the absolute allowance plus the relative allowance times the reference area, about 0.00000766 and 0.00003968 at those extremes.

The comparator matches frames by scheduled threshold and uses each run's frame ID to find its snapshot. It checks actual time first, then checks finite fields, cell IDs, vertex counts, reported areas and vertex displacement with one forward cyclic shift per polygon. Calibration tests reject a four-tolerance shift of the actual time before attempting to open the affected snapshot. They also reject a 0.001-unit vertex movement and a 1% area change to either size extreme. Boundary tests exercise the combined area allowance, Euclidean displacement, ring rotation and reversal, and missing or duplicate data.
