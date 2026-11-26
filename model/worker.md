# Persistent 2D solver

A worker owns one seeded simulation. It accepts a JSON object per line on stdin and flushes one JSON reply per line to stdout. Diagnostics go to stderr. Each process reads the five domain files in its working directory's `initial/` and writes requested frames to its own `output/`. The solver uses one thread per worker.

## Build and run

Run these commands from the repository root on Linux with g++, OpenMP and Python 3.9 or later. The shared library uses the same source as the standalone CLI. The example uses the serial benchmark's allocation capacity and starts in a new directory.

```sh
g++ -std=c++17 -O3 -fopenmp -fPIC -shared -DVERTEX_ARRAY_CAPACITY=1024 \
    model/worker_api.cpp -o /tmp/libvertex2d.so
worker="$(pwd)/model/worker.py"
run_dir="$(mktemp -d)"
cp -r benchmarks/serial2d/initial "$run_dir/"
(cd "$run_dir" && python3 "$worker" /tmp/libvertex2d.so <<'REQUESTS'
{"op":"ping"}
{"op":"reset","seed":646,"time_limit":0.15}
{"op":"step","deadline":0.1}
{"op":"step","deadline":0.2}
{"op":"export","frame_id":1}
{"op":"close"}
REQUESTS
)
```

This sequence returns `ready`, the initial tissue at time 0, the tissue at time 0.1, and the truncated episode at time 0.15. The export writes the eight legacy viewer files for frame 1. Using a fresh directory per episode preserves earlier exports when reusing frame numbers.

## Requests

| Operation | Fields | Reply |
| --- | --- | --- |
| `ping` | none | `ready` |
| `reset` | unsigned 32-bit `seed`, positive finite `time_limit`, optional `parameters` object | Initial tissue and status |
| `step` | absolute simulation time `deadline`, later than the current time | Tissue and status at the earlier of the deadline, time limit or physical ending |
| `measure` | none | Tissue and status at the last commit |
| `export` | positive integer `frame_id` | `ok`, frame ID and actual time |
| `close` | none | `closed`, then process exit |

Reset starts with the compiled mechanical defaults. It accepts the CLI's floating-point mechanical and zone settings, for example `"parameters":{"kA":1.1,"proliferative_CV":1.8}`. Each explicit CV setting also updates the corresponding contact radius through the same setter used by the CLI. Omitted values return to their compiled defaults, including the independent default contact radii. Each reset releases the previous workspace, reads the geometry again, clears the clock and topology/caches, and seeds fresh placement.

A tissue reply contains `time`, `cells` and `measurements`. Each polygon has a solver-slot `id`, integer `zone`, `dead` and `growing` flags, its vertex-average `center`, actual `area`, `target_area` and ordered `vertices`. The live-cell measurements are `cell_count`, `total_area` and `zone_counts` in the five dev zone slots. Polygons flagged dead remain available for inspection and legacy export.

An active episode returns `ok`. An empty tissue or the loss of its last live cell returns `terminated`. Reaching the time limit returns `truncated`. Malformed requests return `error` with a diagnostic. A mechanical or topology failure returns `solver_error`, its diagnostic and the last accepted time and tissue, including earlier accepted substeps of the same request. The worker validates geometry after topology changes and commits the complete observation by value. On failure it releases the trial workspace. Measurement and export use the committed copy. A new `reset` starts the next episode.

The CLI retains dev's adaptive float clock and samples after motion, before topology updates. The worker caps the displacement interval before moving vertices, accumulates decision time in double precision and measures post-topology geometry. Both call the same `Solver2D::reset`, `Solver2D::advance`, `measure_tissue` and `export_tissue` routines.

## Verification

```sh
python3 -m unittest discover -s tests -v
```

The integration tests compile a fixture-sized library and communicate over the real pipes. They check the displacement produced by a shortened step, non-binary decision deadlines, truncation, empty-tissue termination, legacy export, independent seeded reset, parameter restoration, and a failed interval on a 16 × 12 rectangle. An independent worker advances to the reported last successful time and must reproduce every returned cell and tissue measurement. The standalone serial replay and comparator commands are in the [reference benchmark](../benchmarks/serial2d/README.md).
