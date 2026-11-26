"""Experiment configuration and geometry preparation, outside the worker's step loop."""
import copy
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np

MODEL = Path(__file__).resolve().parents[1] / "model"
PARAMETERS = {"growth_rate", "target_size", "division_rate", "death_rate"}
PERMISSIONS = {"allow_growth", "allow_stop", "allow_division", "allow_death"}


def number(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite number, received {value!r}")
    return float(value)


def zone_id(value):
    if type(value) is not int or not 0 <= value < 2**31:
        raise ValueError(f"zone ID must be a nonnegative integer, received {value!r}")
    return value


@dataclass(frozen=True)
class Control:
    zone: int
    parameter: str
    lower: float
    upper: float

    @property
    def identifier(self):
        return f"zone_{self.zone}.{self.parameter}"

    def __post_init__(self):
        zone_id(self.zone)
        if self.parameter not in PARAMETERS:
            raise ValueError(f"unknown decision parameter {self.parameter!r}")
        lo, hi = number(self.lower, "lower"), number(self.upper, "upper")
        if lo < 0 or hi < lo or (self.parameter == "target_size" and lo == 0):
            raise ValueError(f"invalid {self.identifier} bounds [{lo}, {hi}]")


def absolute_paths(obj, base):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key == "path":
                obj[key] = str((base / value).resolve())
            else:
                absolute_paths(value, base)
    elif isinstance(obj, list):
        for value in obj:
            absolute_paths(value, base)


def read_curve(spec):
    if "path" in spec:
        with Path(spec["path"]).open() as f:
            return np.array([[float(x) for x in row] for row in csv.reader(f) if row])
    return np.asarray(spec["points"], dtype=float)


def polygon(spec):
    shape = spec["shape"]
    if shape == "rectangle":
        lo, hi = np.asarray(spec["min"], float), np.asarray(spec["max"], float)
        if lo.shape != (2,) or hi.shape != (2,) or not np.all(hi > lo):
            raise ValueError("rectangle requires increasing minimum and maximum coordinates")
        points = np.array([lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]])
    elif shape in ("circle", "ellipse"):
        center = np.asarray(spec["center"], float)
        radii = np.asarray([spec["radius"]]*2 if shape == "circle" else spec["radii"], float)
        segments = spec.get("segments", 64)
        if (center.shape != (2,) or radii.shape != (2,) or not np.all(radii > 0)
                or type(segments) is not int or segments < 3):
            raise ValueError("invalid curved wall center, radii or segment count")
        theta = np.arange(segments) * (2 * math.pi / segments)
        points = center + radii * np.column_stack((np.cos(theta), np.sin(theta)))
    elif shape == "polygon":
        points = read_curve(spec)
    else:
        raise ValueError(f"unknown 2D wall shape {shape!r}")
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("wall needs finite 2D points")
    if len(points) > 1 and np.array_equal(points[0], points[-1]):
        points = points[:-1]
    if len(points) < 3 or len(np.unique(points, axis=0)) != len(points):
        raise ValueError("wall needs at least three distinct points without repeated vertices")
    area = signed_area(points)
    if area == 0:
        raise ValueError("wall has zero area")
    # Reject crossing/touching nonadjacent sides before native placement.
    for i, (a, b) in enumerate(zip(points, np.roll(points, -1, axis=0))):
        for j in range(i + 2, len(points) - (i == 0)):
            c, d = points[j], points[(j+1) % len(points)]
            if np.any(np.maximum(np.minimum(a,b),np.minimum(c,d)) > np.minimum(np.maximum(a,b),np.maximum(c,d))):
                continue
            def turn(u,v,w):
                x, y = v-u
                return x*(w[1]-u[1])-y*(w[0]-u[0])
            if turn(a,b,c)*turn(a,b,d)<=0 and turn(c,d,a)*turn(c,d,b)<=0:
                raise ValueError("wall polygon intersects itself")
    return points if area > 0 else points[::-1]


def signed_area(points):
    q = points - points[0]
    return float(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1))/2)


class Experiment:
    def __init__(self, *, dimension, geometry, decision_interval, time_limit, zone_measures=None,
                 zones=None, controls=None, binary=None, threads=1, parameters=None):
        if dimension not in (2, 3) or type(dimension) is not int:
            raise ValueError("dimension must be 2 or 3")
        if type(threads) is not int or not 1 <= threads < 2**31:
            raise ValueError("threads must be a positive native integer")
        self.dimension, self.threads = dimension, threads
        self.decision_interval = number(decision_interval, "decision_interval")
        self.time_limit = number(time_limit, "time_limit")
        if min(self.decision_interval, self.time_limit) <= 0:
            raise ValueError("decision_interval and time_limit must be positive")
        if dimension == 2 and self.time_limit > float(np.finfo(np.float32).max):
            raise ValueError("time_limit exceeds the 2D solver clock range")
        self.binary = Path(binary or MODEL / ("vertex3d" if dimension == 3 else "libvertex2d.so")).resolve()
        self.parameters = dict(parameters or {})
        self.initial = self.points = None
        base = Path.cwd()
        if isinstance(geometry, (str, Path)):
            path = Path(geometry).resolve()
            if dimension == 2 and path.is_dir():
                self.initial = path / "initial" if (path / "initial").is_dir() else path
                spec = {}
                base = self.initial.parent
            else:
                spec = json.loads(path.read_text())
                base = path.parent
        else:
            spec = copy.deepcopy(geometry)
        if dimension == 2 and "shape" in spec:
            spec = {"domain":spec}
        absolute_paths(spec, base)
        if dimension == 2 and self.initial is None:
            self.points = polygon(spec["domain"])
        fallback = [{"id": i} for i in range(4)] if self.initial else [{"id":0,"rule":{"kind":"all"}}]
        source_zones = spec.get("zones") or fallback
        source_by_id = {z["id"]: z for z in source_zones}
        self.zones = copy.deepcopy(source_zones if zones is None else [
            {**source_by_id.get(z["id"], {}), **z} for z in zones
        ])
        if self.initial is None and len(self.zones) == 1:
            self.zones[0].setdefault("rule", {"kind":"all"})
        absolute_paths(self.zones, base)
        ids = [zone_id(z["id"]) for z in self.zones]
        if not ids or len(ids) != len(set(ids)):
            raise ValueError("zones must have unique IDs")
        for z in self.zones:
            supported = PARAMETERS | PERMISSIONS | {"id","rule"}
            if dimension == 3:
                supported.add("division_volume")
            if z.keys() - supported:
                raise ValueError(f"unknown zone settings {sorted(z.keys()-supported)}")
            for field in PARAMETERS & z.keys():
                if number(z[field], field) < 0 or (field == "target_size" and z[field] == 0):
                    raise ValueError(f"invalid zone {z['id']} {field}")
            for field in PERMISSIONS & z.keys():
                if type(z[field]) is not bool:
                    raise ValueError(f"{field} must be boolean")
            if dimension == 2 and "rule" in z and "path" in z["rule"]:
                z["rule"]["points"] = read_curve(z["rule"]).tolist()
                del z["rule"]["path"]
        if any("rule" in z for z in self.zones) and not all("rule" in z for z in self.zones):
            raise ValueError("specify a spatial rule for every zone or use legacy 2D classification")
        if dimension == 3:
            if self.parameters:
                raise ValueError("3D mechanical parameters belong in geometry['mechanics']")
            # Spatial rules and standalone mechanics live in the experiment. Training
            # parameters and permissions travel separately in the reset request.
            spec["zones"] = [{k:v for k,v in z.items() if k in {
                "id","rule","allow_growth","allow_division","growth_rate","division_volume","division_rate"
            }} for z in self.zones]
            if len(ids)==1:
                spec["zones"][0].setdefault("rule", {"kind":"all"})
        self.geometry = spec
        if zone_measures is None:
            domain = spec.get("domain", {})
            entire_domain = len(ids)==1 and self.zones[0].get("rule",{}).get("kind")=="all"
            if entire_domain and self.points is not None:
                zone_measures = {ids[0]: signed_area(self.points)}
            elif entire_domain and domain.get("shape")=="box":
                zone_measures = {ids[0]: np.prod(np.subtract(domain["max"],domain["min"]))}
            else:
                raise ValueError("provide zone_measures in simulation area or volume for these regions")
        if set(zone_measures) != set(ids):
            raise ValueError(f"zone_measures must cover zone IDs {ids}")
        self.zone_measures = {i:number(zone_measures[i], f"zone {i} measure") for i in ids}
        if min(self.zone_measures.values()) <= 0:
            raise ValueError("zone measures must be positive")
        self.controls = tuple(controls if controls is not None else
                              [c for i in ids for c in (Control(i,"growth_rate",0,1),Control(i,"target_size",0.01,10))])
        if (any(not isinstance(c,Control) or c.zone not in ids for c in self.controls)
                or len({c.identifier for c in self.controls}) != len(self.controls)):
            raise ValueError("controls must have distinct identifiers in the configured zones")
        if dimension == 2 and any(c.upper > float(np.finfo(np.float32).max) for c in self.controls):
            raise ValueError("control bounds exceed the 2D solver range")

    def stage(self, directory):
        if self.dimension == 3:
            (directory/"experiment.json").write_text(json.dumps(self.geometry, allow_nan=False))
            return [str(self.binary),"worker","--threads",str(self.threads)]
        folder = directory / "initial"
        folder.mkdir()
        if self.initial:
            for name in ("domain_coords","domain_v_ns","domain_v_ts","interface_coords","interface_v_ns"):
                shutil.copyfile(self.initial/(name+".txt"), folder/(name+".txt"))
        else:
            p = self.points
            tangent = np.roll(p,-1,axis=0)-p
            tangent /= np.linalg.norm(tangent,axis=1)[:,None]
            normals = np.column_stack((-tangent[:,1],tangent[:,0]))
            for name, values in (("domain_coords",p),("domain_v_ts",tangent),("domain_v_ns",normals)):
                np.savetxt(folder/(name+".txt"), values, delimiter=",", fmt="%.17g")
            # The legacy loader reads three interface rows before training rules take effect.
            lo,hi = p.min(axis=0),p.max(axis=0)
            y = lo[1]+np.array([.75,.5,.25])*(hi[1]-lo[1])
            anchors = [[x for v in y for x in (end,v)] for end in (lo[0],hi[0])]
            np.savetxt(folder/"interface_coords.txt", anchors, delimiter=",")
            np.savetxt(folder/"interface_v_ns.txt", [[0,1]]*3, delimiter=",")
        return [sys.executable,str(MODEL/"worker.py"),str(self.binary),"--threads",str(self.threads)]

    def reset_payload(self, seed, mode):
        fields = PARAMETERS | PERMISSIONS | {"id"}
        if self.dimension == 2:
            fields.add("rule")
        request = {"op":"reset","seed":seed,"mode":mode,"time_limit":self.time_limit,
                   "zones":[{k:v for k,v in z.items() if k in fields} for z in self.zones]}
        if self.dimension == 3:
            request["experiment"] = "experiment.json"
        else:
            request["parameters"] = self.parameters
        return request
