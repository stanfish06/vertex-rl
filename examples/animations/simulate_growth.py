"""Record a regional growth and division run as compressed JSON.

The square is split into four quadrants by its two midlines and the cube into
eight octants by its three midplanes. Each region has its own target size and
base growth rate. Growth slows as a cell gets larger: the solver divides the
regional rate by 1 + size / growth_half_size. A cell whose measured size passes
divide_at of its regional target divides across its longest axis. Daughters
inherit the growth flag and keep growing. A compressed cell, whose target
exceeds its measured size by more than the inhibit margin, pauses growth until
it has room again. This contact inhibition lets the tissue fill the domain
while capping every cell at the same elastic pressure, so large cells cannot
crush small neighbors. Initial 3D cells are placed at seeded random positions;
2D cells use the solver's seeded random placement. Each JSON sample is one
decision interval of `sample` time units.

Build the workers from the repository root; the 2D capacity must hold every
vertex of the filled square (about 24 per cell):
    make -C model vertex3d
    make -C model libvertex2d.so CPPFLAGS=-DVERTEX_ARRAY_CAPACITY=5000
Then, with vertex_rl installed, the published runs are:
    python examples/animations/simulate_growth.py 2 /tmp/quadrants.json.gz --limit 20
    python examples/animations/simulate_growth.py 3 /tmp/octants.json.gz --limit 4
"""
import argparse
import gzip
import json
import signal
import sys
from pathlib import Path

import numpy as np

from vertex_rl import CellEnv, Experiment, SolverError

SETTINGS = {
    2: dict(
        side=10.0,
        targets=[0.8, 1.3, 1.9, 2.8],
        rates=[1.2, 1.6, 2.0, 2.4],
        divide_at=0.8,
        inhibit=0.12,
        sample=0.05,
        parameters=dict(kA=10, e_g0=0.3, lS=0.1, growth_half_size=1.0, max_edges=24),
    ),
    3: dict(
        side=3.0,
        targets=[0.25, 0.35, 0.45, 0.55, 0.3, 0.5, 0.7, 0.9],
        rates=[0.7, 0.8, 0.9, 1.0, 0.75, 0.85, 0.95, 1.05],
        divide_at=0.8,
        inhibit=0.01,
        sample=0.05,
        initial_cells=14,
        initial_radius=0.25,
        # A level-2 octahedral sphere (66 vertices) for seeds and rebuilt daughters.
        sphere_error=0.07,
        threads=1,
        # Mesh forces act per vertex, so stiffness is high for 66-vertex cells.
        # The barrier contact and wall reaction keep compressed cells apart.
        mechanics=dict(volume_stiffness=3000, surface_tension=1.0, pressure=30, contact_range=0.05,
                       contact_barrier=True, smoothing=2, bin_width=0.5, max_step=0.005,
                       max_displacement=0.005, growth_half_size=0.4, remesh_level=2),
    ),
}
SEEDS = {2: 17, 3: 7}


def random_centers(rng, count, side, radius, gap):
    """Rejection-sample non-overlapping sphere centers clear of the walls."""
    centers = []
    for _ in range(100000):
        p = rng.uniform(radius + gap, side - radius - gap, 3)
        if all(np.linalg.norm(p - q) > 2*radius + gap for q in centers):
            centers.append(p)
            if len(centers) == count:
                return centers
    raise RuntimeError("could not place the requested initial cells")


def experiment(dimension, settings, seed, limit):
    side = settings["side"]
    zones, cells = [], []
    for region in range(2**dimension):
        lo = [((region >> axis) & 1) * side / 2 for axis in range(dimension)]
        hi = [p + side / 2 for p in lo]
        if dimension == 2:
            x, y = lo
            right, top = hi
            rule = dict(kind="surface", inside=True, points=[[x, y], [right, y], [right, top], [x, top]])
        else:
            rule = dict(kind="surface", inside=True, mesh=dict(shape="box", min=lo, max=hi))
        zones.append(dict(id=10 + region, rule=rule, allow_growth=True, allow_division=True,
                          allow_death=False, growth_rate=settings["rates"][region],
                          target_size=settings["targets"][region]))
    shape = "rectangle" if dimension == 2 else "box"
    geometry = dict(domain=dict(shape=shape, min=[0]*dimension, max=[side]*dimension))
    parameters = {}
    if dimension == 3:
        rng = np.random.default_rng(seed)
        radius = settings["initial_radius"]
        error = radius * settings["sphere_error"]
        for center in random_centers(rng, settings["initial_cells"], side, radius, 0.1):
            cells.append(dict(mesh=dict(shape="sphere", center=center.tolist(), radius=radius,
                                        max_surface_error=error)))
        geometry.update(cells=cells, mechanics=settings["mechanics"])
    else:
        parameters = settings["parameters"]
    binary = "vertex3d" if dimension == 3 else "libvertex2d.so"
    model = Path(__file__).resolve().parents[2] / "model"
    return Experiment(dimension=dimension, geometry=geometry, zones=zones,
                      zone_measures={z["id"]: (side/2)**dimension for z in zones},
                      decision_interval=settings["sample"], time_limit=limit, parameters=parameters,
                      binary=model/binary, threads=settings.get("threads", 1))


def long_axis(vertices):
    points = np.asarray(vertices, dtype=float)
    points = points - points.mean(axis=0)
    _, vectors = np.linalg.eigh(points.T @ points)
    return vectors[:, -1]


def interval_actions(cells, settings, dimension):
    actions = {}
    for c in cells:
        actions[c["id"]] = "stop_growing" if c["target"] > c["size"] + settings["inhibit"] else "grow"
        if c["size"] > settings["divide_at"] * settings["targets"][c["zone"]]:
            axis = long_axis(c["vertices"])
            if dimension == 2:
                # The 2D angle is the direction of the cut line.
                orientation = float(np.arctan2(axis[1], axis[0]) + np.pi/2)
            else:
                orientation = axis.tolist()
            actions[c["id"]] = dict(fate="divide", orientation=orientation)
    return actions


def record(dimension, destination, limit, seed, overrides):
    settings = {**SETTINGS[dimension], **overrides}
    env = CellEnv(experiment(dimension, settings, seed, limit))
    frames = []
    size_key = "area" if dimension == 2 else "volume"
    ending = dict(reason="time_limit")
    try:
        env.reset(seed=seed)
        while True:
            cells = []
            for c in env._state["cells"]:
                cell = dict(id=c["id"], zone=c["zone"]-10, center=c["center"], vertices=c["vertices"],
                            size=c[size_key], target=c["target_"+size_key], growing=c["growing"])
                if dimension == 3:
                    cell["faces"] = c["faces"]
                cells.append(cell)
            frames.append(dict(time=env.time, cells=cells))
            if not env.agents:
                ending = dict(reason=env._state["status"])
                break
            if env.time >= limit - 1e-9:
                break
            if len(frames) % 20 == 0:
                print(f"t={env.time:.2f} cells={len(cells)}", file=sys.stderr, flush=True)
            try:
                env.step(interval_actions(cells, settings, dimension))
            except SolverError as error:
                ending = dict(reason="solver_error", time=error.time, message=str(error))
                break
            except KeyboardInterrupt:
                ending = dict(reason="interrupted")
                break
    finally:
        env.close()
    with gzip.open(destination, "wt") as stream:
        json.dump(dict(d=dimension, side=settings["side"], seed=seed, settings=settings,
                       frames=frames, ending=ending), stream, separators=(",", ":"), allow_nan=False)
    print(f"{len(frames)} samples, {len(cells)} cells, t={frames[-1]['time']:.3g}, {ending}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dimension", type=int, choices=(2, 3))
    parser.add_argument("destination", type=Path)
    parser.add_argument("--limit", type=float, default=12.)
    parser.add_argument("--seed", type=int, help="placement and solver seed (default 17 in 2D, 7 in 3D)")
    parser.add_argument("--settings", type=json.loads, default={},
                        help="JSON object overriding entries of SETTINGS for this dimension")
    args = parser.parse_args()
    # Keep the frames recorded so far when a time limit terminates the run.
    signal.signal(signal.SIGTERM, signal.default_int_handler)
    seed = SEEDS[args.dimension] if args.seed is None else args.seed
    record(args.dimension, args.destination, args.limit, seed, args.settings)
