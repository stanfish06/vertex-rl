# 3D vertex simulation

The 3D solver evolves separate closed triangle meshes inside a fixed triangulated wall. The CLI writes a time series for ParaView. The persistent worker exposes the same mechanics, division and export routines at requested simulation times.

## Build and run

Install a C++17 compiler, OpenMP, make and the json-c development package. On Ubuntu the latter is `libjson-c-dev`. The combined 2D/3D CLI also needs GNU Readline development files.

From the repository root:

```sh
make -C model vertex3d
model/vertex3d run examples/3d/cube.json /tmp/vertex-cube --threads 4
```

Open `/tmp/vertex-cube/tissue.pvd` in ParaView. The example starts with a cell of volume 1 and surface area 6, bisects it at time 0.05, and advances the daughters through time 0.1. At the division frame, threshold the cell-data `cell_id` array to select either daughter. IDs 2 and 3 have disjoint point indices, volumes of 0.5 and zone IDs 10 and 20. Each output frame is VTK XML PolyData. Its polygons carry `cell_id`, `zone_id` and `volume`.

For a standalone CLI that selects either dimension, build `make -C model vertex` and run `model/vertex --3d run examples/3d/cube.json /tmp/vertex-cube`. Invoking `model/vertex` directly opens the existing 2D CLI. Set `CPPFLAGS=-DVERTEX_ARRAY_CAPACITY=1024` when building a fixture-sized 2D CLI or worker.

## Configure the experiment

[The example](../../examples/3d/cube.json) gives a complete JSON configuration. Coordinates and surface errors use simulation length units. Each cell's optional positive `target_volume` defaults to its initial volume. `growing` defaults to true. Supply an empty `cells` array for an empty tissue. The seed controls the built-in division draws and orientations. Initial meshes are explicitly specified.

`domain` and each cell's `mesh` accept these shapes:

| Shape | Fields |
| --- | --- |
| `box` | `min`, `max` 3-vectors |
| `sphere` | `center`, `radius`, `max_surface_error` |
| `ellipsoid` | `center`, `radii`, `max_surface_error` |
| `cylinder` | `center`, `radius`, `height`, `max_surface_error`; axis parallel to z |
| `mesh` | `vertices` as 3-vectors and `faces` as triples of zero-based indices |
| `import` | `path` to an outward-wound OBJ or STL surface |

The facet interior and edge deviation determine curved-wall resolution. Generators reject curved shapes requiring more than two million facets. Wall imports accept closed, connected, outward-wound OBJ and ASCII or binary STL files. OBJ faces may be simple planar concave polygons. For file-backed experiments, relative import paths resolve against the configuration's directory. Inline worker experiments resolve paths against the process working directory. At reset every cell must fit wholly inside the wall and be disjoint from other initial meshes.

Zones form an ordered list, each with an integer `id` and a `rule`. Classification uses the volume centroid and chooses the first matching rule, including on a shared interface. The default covers the domain as zone 0. Supply rules that cover every location a centroid can reach.

| Rule `kind` | Fields |
| --- | --- |
| `all` | matches every point |
| `layer` | `origin`, `axis`, optional inclusive `lower` and `upper` signed distances |
| `radial` | `center`, optional inclusive `lower` and `upper` radii |
| `surface` | closed `mesh` in any shape format above, and `inside` boolean |
| `interface` | `path` to an oriented OBJ/STL interface and `positive_side` boolean |

A zone's `growth_rate` adds target volume per unit accepted time when the cell is growing and `allow_growth` is true. An eligible cell with actual volume at or above `division_volume` divides in an accepted interval of duration Δt with probability `1 − exp(−division_rate × Δt)` when `allow_division` is true. Both rates default to zero. Built-in divisions draw a uniformly oriented plane through the volume centroid. Each daughter gets a new ID, the parent's growth flag and its share of target volume in proportion to measured volume. Subsequent growth uses its own zone's permissions. The standalone `divisions` list can request a plane at a specified time and emits a frame after the cut. The cap supports sections bounded by simple loops, including concave loops. A section with a hole, a plane coinciding with an existing face, or a split producing disconnected or nonmanifold daughters rejects the cut and preserves the parent. A rejected scheduled cut stops the standalone run with a diagnostic.

The optional `mechanics` object sets `volume_stiffness` (default 1), `surface_tension` (0.05), contact `pressure` at zero gap (1), `contact_range` (0.1), broad-phase grid `bin_width` (1), `max_displacement` (0.05), `max_step` (0.01) and `min_step` (1e-8). Force equals the negative gradient of the volume-elastic and surface energies plus repulsive cell contact. Overdamped motion uses unit mobility. The solver shortens a crossing or invalid mesh trial until it can commit motion or reaches the minimum step. A final interval shorter than the minimum can land exactly on a deadline.

## Persistent worker

```sh
model/vertex3d worker /tmp/vertex-frames --threads 2 <<'REQUESTS'
{"op":"ping"}
{"op":"reset","seed":646,"time_limit":1,"experiment":"examples/3d/cube.json"}
{"op":"divide","cell_id":1,"normal":[1,0,0]}
{"op":"export","frame_id":1}
{"op":"step","deadline":0.025}
{"op":"close"}
REQUESTS
```

The common launcher also supports `python3 model/worker.py model/vertex3d --dimension 3 --threads 2`. It uses `output/` in the worker's current directory. Each request and reply occupies one JSON line. Use a separate output directory for each run whose frames you want to keep.

A reset accepts either an experiment path or its JSON object. The request's seed and limit govern that episode. The worker waits for explicit `divide` and `export` requests for interventions and frames. CLI sampling and event fields in experiment files apply to standalone runs. The reply includes the committed time, per-cell geometry and target volumes, tissue volume, cell and zone counts, and cumulative counts of accepted steps and rejected trials preceding accepted steps. Zone counts follow configuration order. A requested cut that a zone or the geometry rejects returns `division_accepted: false` with the parent intact.

Step deadlines are absolute, capped at the episode limit. Empty tissue returns `terminated`, and the time limit returns `truncated`. A failed trial returns `solver_error` with the committed time and geometry. Measurement and export remain available. Reset explicitly to begin the next episode. An invalid request returns `error`.

The training reset modes add bounded zone controls or a single action per live cell at each decision boundary. [Training worker protocol](../worker.md#training-intervals-in-both-dimensions).

## Verification

```sh
python3 -m unittest discover -s tests -v
```

Native and integration tests cover closed-mesh mechanics, plane cuts, centroid classification, curved and imported walls, contact and swept crossings. Integration checks read the daughter fields through the VTK collection, compare tissue states at one, two and four threads, and reproduce the last accepted geometry after a wall-crossing failure. The generated box and an imported copy of its faces must give identical accept/retry counts and final state.
