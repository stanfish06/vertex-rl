# Training qualification results

On November 26, 2026, a feedback policy trained on three seeds reproduced the known growth-rate target on five unseen initial geometries. The [protocol](protocol.md) fixed the tissue, search budgets, seeds and scoring rules before evaluation.

| Held-out seed | Trained policy error | Default error | Random-search error |
| ---: | ---: | ---: | ---: |
| 211 | 0.0000047505 | 0.046542 | 0.0016117 |
| 223 | 0.0000048480 | 0.046857 | 0.0016121 |
| 227 | 0.0000050114 | 0.048294 | 0.0016706 |
| 229 | 0.0000051098 | 0.047869 | 0.0016418 |
| 233 | 0.0000049109 | 0.048887 | 0.0017169 |
| **Equal-weight mean** | **0.0000049261** | **0.047690** | **0.0016506** |

Policy search and random search each ran 192 training episodes. All 15 held-out episodes reached the time limit and supplied eight scored intervals each. The reference growth rate was 0.73. The learned policy chose rates from 0.729438 to 0.730068 across held-out decisions. Random search chose a fixed rate of 0.711621. The configured default was 0.15. These scores apply to recovery of one growth rate on a pair of separated 3D cells.

## Four-worker throughput

| Active environments | Solver threads per worker | Simulated time | Elapsed seconds | Simulated time per second |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 4 | 0.222980 | 17.94 |
| 4 | 1 | 4 | 0.054126 | 73.90 |

All four workers replied within the common 30-second readiness deadline in each batch. Timing covered the first reset through the last step return. Each of the four seeds gave exactly equal observations at reset and at all eight decisions, and all episodes reached the time limit. The concurrent-to-sequential rate ratio was **4.12×** for this single short paired measurement.

The host was `lima-devflow-task-dda41e3bf92e`, an Ubuntu VM with four available cores on a Ryzen 7 8845HS. The run used g++ 15.2.0, Python 3.14.4, NumPy 2.5.5, Gymnasium 1.4.0 and one OpenBLAS thread. The [qualification archive](qualification.tar.gz) contains the machine records and native executable hash, reference trajectories, action recordings, scores for every candidate, frozen selections, interval-level held-out returns and both timed replays. In the captured learning output, `truncated` indicates reaching the episode time limit.
