"""Occupied-shape measurements, held reference samples and time-weighted scores."""
from bisect import bisect_right
from copy import deepcopy
from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class ShapeMoments:
    size: float
    centroid: np.ndarray
    covariance: np.ndarray
    elongation: float
    axis: np.ndarray
    orientation: np.ndarray


def shape_moments(cell):
    """Integrate over a polygon interior or a consistently wound closed triangle mesh."""
    vertices = np.asarray(cell["vertices"], dtype=float)
    if (vertices.ndim != 2 or vertices.shape[1] not in (2, 3) or
            len(vertices) < vertices.shape[1] + 1 or not np.isfinite(vertices).all()):
        raise ValueError("invalid cell vertices")
    dimension = vertices.shape[1]
    origin = vertices[0]
    offsets = vertices - origin
    if dimension == 2:
        simplices = np.stack((offsets, np.roll(offsets, -1, axis=0)), axis=1)
        a, b = simplices.transpose(1, 0, 2)
        weights = (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]) / 2
    else:
        faces = np.asarray(cell["faces"])
        if (faces.ndim != 2 or faces.shape[1] != 3 or not len(faces) or faces.dtype.kind not in 'iu'
                or faces.min() < 0 or faces.max() >= len(vertices)):
            raise ValueError("invalid cell faces")
        simplices = offsets[faces]
        a, b, c = simplices.transpose(1, 0, 2)
        weights = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6
    signed_size = weights.sum()
    if not np.isfinite(signed_size) or signed_size == 0:
        raise ValueError("cell has zero or nonfinite size")
    weights = weights / signed_size
    sums = simplices.sum(axis=1)
    centroid = weights @ sums / (dimension + 1)
    # In a uniform d-simplex, E[b_i*b_j] is (1 + (i==j))/((d+1)(d+2))
    # for barycentric coordinates b. Signed simplices also cover concave cells.
    second = np.einsum("n,ni,nj->ij", weights, sums, sums)
    second += np.einsum("n,nki,nkj->ij", weights, simplices, simplices)
    covariance = second / ((dimension + 1) * (dimension + 2)) - np.outer(centroid, centroid)
    values, vectors = np.linalg.eigh(covariance)
    if values[-2] <= 0 or not np.isfinite(values).all():
        raise ValueError("cell has singular shape covariance")
    elongation = math.sqrt(values[-1] / values[-2])
    # Account for round-off at the inclusive 5% tie boundary.
    if (values[-1] - values[-2]) / values[-1] <= 0.05 + 32 * np.finfo(float).eps:
        axis = np.zeros(dimension)
        orientation = np.eye(dimension) / dimension
    else:
        axis = vectors[:, -1]
        if axis[np.argmax(np.abs(axis))] < 0:
            axis = -axis
        orientation = np.outer(axis, axis)
    return ShapeMoments(float(abs(signed_size)), origin + centroid, covariance,
                        elongation, axis, orientation)


def _mean_tensor(tensors):
    return np.mean(tensors, axis=0).tolist() if tensors else None


def measure_tissue(cells, zone_measures, dimension):
    """Return JSON-serializable measurements in configured zone order."""
    if dimension not in (2, 3) or not zone_measures:
        raise ValueError("a tissue needs dimension 2 or 3 and at least one zone")
    records = {zone: ([], []) for zone in zone_measures}
    if not all(math.isfinite(m) and m > 0 for m in zone_measures.values()):
        raise ValueError("zone areas or volumes must be positive and finite")
    cells = list(cells)
    shapes = [shape_moments(cell) for cell in cells]
    centers = np.array([s.centroid for s in shapes])
    for i, (cell, shape) in enumerate(zip(cells, shapes)):
        if len(shape.centroid) != dimension or cell['zone'] not in records:
            raise ValueError("cell has wrong dimension or unknown zone")
        orientations, arrangements = records[cell['zone']]
        orientations.append(shape.orientation)
        delta = centers - centers[i]
        squared = np.einsum("ij,ij->i", delta, delta)
        distinct = squared > 0
        if np.any(distinct):
            shortest = np.min(squared[distinct])
            # Zero displacement has no direction. Among distinct centers give
            # all nearest neighbors one equal share of the focal cell's vote.
            tied = np.isclose(squared, shortest, rtol=32*np.finfo(float).eps, atol=0)
            directions = delta[tied] / np.sqrt(squared[tied, None])
            arrangements.append(directions.T @ directions / len(directions))
    zones = []
    for zone, (orientations, arrangements) in records.items():
        count = len(orientations)
        zones.append(dict(id=zone, count=count, density=count / zone_measures[zone],
                          orientation=_mean_tensor(orientations), arrangement=_mean_tensor(arrangements)))
    sizes = [s.size for s in shapes]
    return dict(count=len(cells), size_percentiles=np.percentile(sizes, [10, 50, 90]).tolist()
                if sizes else None, zones=zones)


def _tensor_error(observed, target):
    if observed is None or target is None:
        return float((observed is None) != (target is None))
    a, b = np.asarray(observed), np.asarray(target)
    if a.shape != b.shape:
        raise ValueError("tensor dimensions differ")
    return float(np.linalg.norm(a - b)) / math.sqrt(2)


def _zone_table(sample, expected):
    zones = {zone["id"]: zone for zone in sample["zones"]}
    if len(zones) != len(sample["zones"]) or zones.keys() != expected.keys():
        raise ValueError("measurement zone IDs must match the configured zones exactly")
    return zones


def _validate_sample(sample, zone_measures):
    count = sample["count"]
    if type(count) is not int or count < 0:
        raise ValueError("reference cell count must be a nonnegative integer")
    if count:
        sizes = np.asarray(sample['size_percentiles'], dtype=float)
        if (sizes.shape != (3,) or not np.isfinite(sizes).all() or np.any(sizes <= 0)
                or np.any(np.diff(sizes) < 0)):
            raise ValueError("reference size percentiles must be positive, finite and ordered")
    counts = []
    for zone in _zone_table(sample, zone_measures).values():
        n = zone['count']
        if type(n) is not int or n < 0 or not math.isfinite(zone['density']) or zone['density'] < 0:
            raise ValueError("invalid reference zone count or density")
        counts.append(n)
        if bool(n) != (zone['orientation'] is not None):
            raise ValueError("reference orientation is required exactly for populated zones")
        for name in ('orientation', 'arrangement'):
            if zone[name] is not None:
                tensor = np.asarray(zone[name], dtype=float)
                if tensor.shape not in ((2,2),(3,3)) or not np.isfinite(tensor).all():
                    raise ValueError(f"invalid reference {name} tensor")
    if sum(counts) != count:
        raise ValueError("reference zone counts must sum to the tissue count")


def error_components(observed, target, zone_measures):
    """Unweighted component errors for the selected reference sample."""
    count, target_count = observed['count'], target['count']
    if count and target_count:
        p = np.asarray(observed['size_percentiles'])
        q = np.asarray(target['size_percentiles'])
        if p.shape != (3,) or q.shape != (3,) or np.any(q <= 0):
            raise ValueError("three positive reference size percentiles are required")
        size_error = float(np.mean(np.abs(p - q) / q))
    else:
        size_error = float(bool(count) != bool(target_count))
    actual, reference = _zone_table(observed, zone_measures), _zone_table(target, zone_measures)
    density, orientation, arrangement = [], [], []
    for zone, measure in zone_measures.items():
        a, b = actual[zone], reference[zone]
        density.append(abs(a['density'] - b['density']) / max(b['density'], 1 / measure))
        orientation.append(_tensor_error(a['orientation'], b['orientation']))
        arrangement.append(_tensor_error(a['arrangement'], b['arrangement']))
    return dict(count=abs(count-target_count) / max(target_count, 1), size=size_error,
                density=float(np.mean(density)), orientation=float(np.mean(orientation)),
                arrangement=float(np.mean(arrangement)))


class TargetScorer:
    """Match each end time to its latest reference row at or before that time."""

    def __init__(self, samples, time_limit, zone_measures, reset_time=0.):
        if not (math.isfinite(reset_time) and math.isfinite(time_limit) and time_limit > reset_time):
            raise ValueError("time limit must be finite and later than reset")
        self.zone_measures = dict(zone_measures)
        if not self.zone_measures or not all(math.isfinite(v) and v > 0 for v in self.zone_measures.values()):
            raise ValueError("zone measures must be positive and finite")
        rows = []
        seen = {}
        for row, sample in enumerate(samples, 1):
            t = float(sample['time'])
            if not math.isfinite(t):
                raise ValueError(f"reference row {row} has nonfinite time {t}")
            if t in seen:
                raise ValueError(f"duplicate timestamp {t} in rows {seen[t]} and {row}")
            seen[t] = row
            measurements = deepcopy(sample['measurements'])
            _validate_sample(measurements, self.zone_measures)
            rows.append((t, measurements))
        if not rows:
            raise ValueError("target reference has no samples")
        rows.sort(key=lambda item: item[0])
        self.times, self.samples = zip(*rows)
        if self.times[0] > reset_time:
            raise ValueError(f"first reference time {self.times[0]} is after reset time {reset_time}")
        if self.times[-1] < time_limit:
            raise ValueError(f"last reference time {self.times[-1]} is before episode limit {time_limit}")

    def at(self, time):
        if not math.isfinite(time) or not self.times[0] <= time <= self.times[-1]:
            raise ValueError(f"time {time} outside reference coverage")
        return self.samples[bisect_right(self.times, time)-1]

    def __call__(self, measurements, time):
        errors = error_components(measurements, self.at(time), self.zone_measures)
        score = -(errors['count'] + errors['size'] + errors['density'] +
                  (errors['orientation'] + errors['arrangement']) / 2) / 4
        if not math.isfinite(score):
            raise ValueError("nonfinite target matching reward")
        return score


def episode_error(intervals):
    """Score one completed episode from (start time, actual end time, reward)."""
    durations, costs = [], []
    for start, end, reward in intervals:
        if not all(map(math.isfinite, (start, end, reward))) or end < start:
            raise ValueError("score needs finite times and reward and a nonnegative duration")
        durations.append(end-start)
        costs.append(-(end-start)*reward)
    duration = math.fsum(durations)
    if not duration:
        raise ValueError("episode needs positive scored duration")
    return math.fsum(costs)/duration
