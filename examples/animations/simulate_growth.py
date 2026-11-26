"""Record native cell geometry, target sizes and the end reason as compressed JSON.

With vertex_rl installed and both workers built, run from the repository root:
    python examples/animations/simulate_growth.py 3 /tmp/octants.json.gz
"""
import argparse
import gzip
import json
from pathlib import Path

from vertex_rl import CellEnv, Experiment, SolverError


TARGETS = (1.35, 2.25, 3.50, 5.20, 1.75, 2.85, 4.25, 5.80)
RATES = (0.60, 1.95, 4.20, 8.40, 1.20, 2.85, 6.15, 9.60)
STEP = 0.1


def experiment(dimension, limit):
    side = 16 if dimension == 2 else 10
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
            for offset in [(1.55, 1.65, 2.7), (3.5, 3.45, 1.85)]:
                center = [a + b for a, b in zip(lo, offset)]
                cells.append(dict(mesh=dict(shape="sphere", center=center, radius=.52, max_surface_error=.10)))
        zones.append(dict(id=10 + region, rule=rule, allow_growth=True, allow_division=True,
                          allow_death=False, growth_rate=RATES[region], target_size=TARGETS[region]))
    shape = "rectangle" if dimension == 2 else "box"
    geometry = dict(domain=dict(shape=shape, min=[0]*dimension, max=[side]*dimension))
    parameters = {}
    if dimension == 3:
        geometry.update(cells=cells, mechanics=dict(volume_stiffness=3, surface_tension=.05,
                        pressure=3, max_step=.025, max_displacement=.015, contact_range=.07))
    else:
        parameters = dict(kA=2., resting_T=.12)
    binary = "vertex3d" if dimension == 3 else "libvertex2d.so"
    model = Path(__file__).resolve().parents[2] / "model"
    return side, Experiment(dimension=dimension, geometry=geometry, zones=zones,
                            zone_measures={z["id"]: (side/2)**dimension for z in zones},
                            decision_interval=STEP, time_limit=limit, parameters=parameters,
                            binary=model/binary)


def interval_actions(cells, frame, dimension):
    # Birth IDs are deliberately treated exactly like other live IDs.
    actions = {c["id"]: "grow" for c in cells}
    if frame >= 10 and frame % 4 == 0:
        for region in range(2**dimension):
            mature = [c for c in cells if c["zone"] == region and c["size"] > .76*TARGETS[region]]
            if mature:
                parent = max(mature, key=lambda c: c["size"])
                axis = (frame//4 + region) % dimension
                orientation = axis*1.5707963268 if dimension == 2 else [int(k == axis) for k in range(3)]
                actions[parent["id"]] = dict(fate="divide", orientation=orientation)
    return actions


def record(dimension, destination, limit):
    side, config = experiment(dimension, limit)
    seed = 53 if dimension == 2 else 646
    env = CellEnv(config)
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
            if env.time >= limit:
                break
            try:
                env.step(interval_actions(cells, len(frames)-1, dimension))
            except SolverError as error:
                ending = dict(reason="solver_error", time=error.time, message=str(error))
                break
    finally:
        env.close()
    count = 2**dimension
    with gzip.open(destination, "wt") as stream:
        json.dump(dict(d=dimension, side=side, seed=seed, targets=TARGETS[:count], rates=RATES[:count],
                       frames=frames, ending=ending), stream, separators=(",", ":"), allow_nan=False)
    print(f"{len(frames)} samples, {len(cells)} cells, t={frames[-1]['time']:.3g}, {ending['reason']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dimension", type=int, choices=(2, 3))
    parser.add_argument("destination", type=Path)
    parser.add_argument("--limit", type=float, default=120.)
    args = parser.parse_args()
    record(args.dimension, args.destination, args.limit)
