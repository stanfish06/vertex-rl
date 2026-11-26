# Persistent 2D solver

A worker owns one seeded simulation. It accepts a JSON object per line on stdin and flushes one JSON reply per line to stdout. Diagnostics go to stderr. Each process reads the five domain files in its working directory's `initial/` and writes requested frames to its own `output/`. The solver uses one thread per worker by default. Supply a positive integer with `--threads N` at process launch to use the same force-thread budget across resets. Built-in fate draws and topology mutations run in serial order.

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

## Training intervals in both dimensions

A training reset sets `mode` to `tissue` or `cell`. A 3D reset also supplies its `experiment`. An optional `zones` array specifies overrides by zone ID. Growth, stop, division and death permissions are named `allow_growth`, `allow_stop`, `allow_division` and `allow_death`. The native reports expose a lifetime `id` and `age` in simulation time for each live cell. The Python cell policy receives its own geometry-derived state and zone. Tissue scoring can use the complete geometry.

A tissue step accepts absolute assignments in the declared control order, for example

```json
{"op":"step","deadline":0.1,"controls":[{"zone":1,"parameter":"growth_rate","value":0.5},{"zone":1,"parameter":"target_size","value":2}]}
```

The four selectable parameters are `growth_rate`, `target_size`, `division_rate` and `death_rate`. The worker checks the whole batch before assignment. Python checks researcher-configured bounds before dispatch. `target_size` is the zone’s growth ceiling in area or volume units. An action below an already attained target stops further growth until the ceiling rises. Growth adds target volume at the configured rate in 3D. In 2D its rate tapers linearly to zero at the ceiling. Repeating the action preserves accumulated growth. Built-in division and death hazards use the accepted simulated duration. Zone eligibility governs these built-in events.

A cell step accepts exactly one choice per live agent. In 2D the division orientation is the angle of the cut line in radians. In 3D it is a nonzero plane normal. For example, a 3D worker with a single living agent accepts

```json
{"op":"step","deadline":0.1,"actions":[{"id":1,"action":"divide","orientation":[1,0,0]}]}
```

The other choices are `grow`, `stop_growing` and `die`. An eligible physical cut receives one attempt at the decision boundary. Daughters inherit the growth flag and split the parent's target in proportion to their measured size. A successful cut retires the parent ID and allocates two new IDs. A rejected division keeps the parent active. Death retires the agent immediately. If all agents die at the boundary, the episode ends at that time. Built-in fates are disabled while the cell policy is attached.

For a 2D training experiment the ordered zone rules support `all`, `layer` with origin and axis, `radial` with center, `interface` with oriented polyline points and `positive_side`, and `surface` with closed polygon points and `inside`. Layers and radial bands accept inclusive lower and upper limits. A polyline's positive side lies to the left of its nearest directed segment. A generic polygon wall uses point-in-polygon membership and closest-segment boundary projection. The legacy material tables provide properties for zones 0 through 4. Additional zone IDs use the base material and their own fate permissions and controls.

Early policy divisions can create daughters smaller than the widths of contact layers. For training-mode planar contacts, each vertex interacts with the nearest feature of the other polygon. Interior repulsion requires containment. The worker computes contact reactions and polygon areas in local coordinates and checks simple positive-area polygons before committing a training step. Planar cut faces use a small initial clearance capped at 0.001 simulation length units. Standalone runs and workers without a training mode use the reference contact and fate calculations.
