"""Validation and native layout for the 2D training worker."""
import ctypes as C
import math


def real(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def identifier(value):
    if type(value) is not int or not 0 <= value < 2**31:
        raise ValueError("invalid ID")
    return value


def fields(obj, allowed):
    if not isinstance(obj, dict) or obj.keys() - set(allowed):
        raise ValueError(f"expected fields from {allowed}")


def boolean(value):
    if type(value) is not bool:
        raise ValueError("expected boolean")
    return value


def point(value):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("expected a 2D point")
    return [real(x, "coordinate") for x in value]


class ZoneSpec(C.Structure):
    _fields_ = ([(name, C.c_int) for name in ("id", "rule", "flags", "positive", "count")]
                + [(name, C.c_double) for name in ("growth", "target", "division", "death", "lower", "upper")]
                + [("points", C.POINTER(C.c_double))])


PARAMETERS = ("growth_rate", "target_size", "division_rate", "death_rate")
FATES = ("grow", "stop_growing", "divide", "die")
PERMISSIONS = ("allow_growth", "allow_stop", "allow_division", "allow_death")


def parameter_value(name, value):
    value = real(value, name)
    if value < 0 or (name == "target_size" and value == 0) or value > 3.4028234663852886e38:
        raise ValueError(f"invalid {name}")
    return value


def rule_config(obj):
    kind = obj["kind"]
    points = []
    lower, upper, positive = -math.inf, math.inf, True
    if kind == "all":
        fields(obj, ["kind"])
    elif kind in ("layer", "radial"):
        if kind == "layer":
            fields(obj, ["kind", "origin", "axis", "lower", "upper"])
            axis = point(obj["axis"])
            length = math.hypot(*axis)
            if not math.isfinite(length) or not length:
                raise ValueError("invalid layer axis")
            points = [point(obj["origin"]), [x / length for x in axis]]
        else:
            fields(obj, ["kind", "center", "lower", "upper"])
            points = [point(obj["center"])]
            lower = 0
        if "lower" in obj:
            lower = real(obj["lower"], "lower")
        if "upper" in obj:
            upper = real(obj["upper"], "upper")
        if lower > upper or (kind == "radial" and lower < 0):
            raise ValueError("invalid zone limits")
    elif kind in ("interface", "surface"):
        flag = "positive_side" if kind == "interface" else "inside"
        fields(obj, ["kind", "points", flag])
        points = [point(p) for p in obj["points"]]
        positive = boolean(obj.get(flag, True))
        if len(points) < (2 if kind == "interface" else 3) or any(a == b for a, b in zip(points, points[1:])):
            raise ValueError("invalid zone polyline")
    else:
        raise ValueError(f"unknown zone rule {kind}")
    index = ("legacy", "all", "layer", "radial", "interface", "surface").index(kind)
    return index, points, lower, upper, positive


def zones_config(zones):
    if not isinstance(zones, list):
        raise ValueError("zones must be a list")
    selected = {}
    for z in zones:
        fields(z, ["id", "rule", *PARAMETERS, *PERMISSIONS])
        key = identifier(z["id"])
        if key in selected:
            raise ValueError("duplicate zone ID")
        selected[key] = z
    custom = any("rule" in z for z in zones)
    if custom and not all("rule" in z for z in zones):
        raise ValueError("all configured zones require rules")
    if not custom and set(selected) - set(range(5)):
        raise ValueError("legacy zones use IDs 0 through 4")
    ids = list(selected) if custom else list(range(5))
    specs = (ZoneSpec * len(ids))()
    for i, key in enumerate(ids):
        z = selected.get(key, {})
        permissions = (True, True, custom or key == 1, custom or key == 3)
        flags = sum(1 << j for j, (name, fallback) in enumerate(zip(PERMISSIONS, permissions))
                    if boolean(z.get(name, fallback)))
        params = [parameter_value(name, z[name]) if name in z else default
                  for name, default in zip(PARAMETERS, (5, 1.5 if custom else 0, 0.8, 0.1))]
        rule, points, lower, upper, positive = rule_config(z["rule"]) if custom else (0, [], 0, 0, True)
        coords = (C.c_double * (2 * len(points)))(*(x for p in points for x in p))
        specs[i] = ZoneSpec(key, rule, flags, positive, len(points), *params, lower, upper, coords)
    return ids, specs


def controls(items, ids):
    if not isinstance(items, list):
        raise ValueError("controls must be a list")
    result, seen = [], set()
    for item in items:
        fields(item, ["zone", "parameter", "value"])
        zone, name = identifier(item["zone"]), item["parameter"]
        if zone not in ids or name not in PARAMETERS or (zone, name) in seen:
            raise ValueError(f"unknown or duplicate control {zone}.{name}")
        seen.add((zone, name))
        result.append((zone, PARAMETERS.index(name), parameter_value(name, item["value"])))
    return result


def actions(items, live):
    if not isinstance(items, list) or len(items) != len(live):
        raise ValueError("expected one action per live cell")
    result = {}
    for item in items:
        name, key = item["action"], identifier(item["id"])
        fields(item, ["id", "action", "orientation"] if name == "divide" else ["id", "action"])
        if key not in live or key in result or name not in FATES:
            raise ValueError("unknown or duplicate action/agent ID")
        angle = real(item["orientation"], "angle") if name == "divide" else 0
        result[key] = (FATES.index(name), angle)
    return [(key, *result[key]) for key in sorted(result)]
