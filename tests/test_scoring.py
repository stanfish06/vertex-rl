import copy
import math
import unittest

import numpy as np

from vertex_rl.scoring import TargetScorer, episode_error, error_components, measure_tissue, shape_moments


def rectangle(x=0, y=0, w=2., h=1., angle=0., zone=1):
    points = np.array([[-w,-h], [w,-h], [w,h], [-w,h]]) / 2
    c, s = math.cos(angle), math.sin(angle)
    return dict(vertices=points @ [[c,s],[-s,c]] + [x,y], zone=zone)


def cuboid():
    points = np.array([[0,0,0], [2,0,0], [2,1,0], [0,1,0],
                       [0,0,1], [2,0,1], [2,1,1], [0,1,1]], dtype=float)
    faces = [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
             [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]
    return dict(vertices=points, faces=faces, zone=1)


class ScoringTest(unittest.TestCase):
    def test_spatial_moments_and_retessellation(self):
        a = rectangle(x=10, y=20)
        initial = shape_moments(a)
        self.assertAlmostEqual(initial.size, 2)
        np.testing.assert_allclose(initial.centroid, [10,20])
        np.testing.assert_allclose(initial.covariance, np.diag([4/12,1/12]), atol=1e-14)
        v = a['vertices']
        a['vertices'] = np.vstack([v[0]+(v[1]-v[0])*i/10 for i in range(10)] + list(v[1:]))
        np.testing.assert_allclose(shape_moments(a).covariance, initial.covariance, atol=1e-14)
        b = cuboid()
        initial = shape_moments(b)
        self.assertAlmostEqual(initial.size, 2)
        np.testing.assert_allclose(initial.covariance, np.diag([4/12,1/12,1/12]), atol=1e-14)
        self.assertAlmostEqual(initial.elongation, 2)
        face = b['faces'].pop()
        b['vertices'] = np.vstack([b['vertices'], np.mean(b['vertices'][face], axis=0)])
        b['faces'].extend([[i,j,8] for i,j in zip(face, face[1:]+face[:1])])
        np.testing.assert_allclose(shape_moments(b).orientation, initial.orientation, atol=1e-14)
        np.testing.assert_allclose(shape_moments(b).centroid, [1,.5,.5])

    def test_isotropy_threshold_and_reversed_axis(self):
        for scale in (.01, 1, 100):
            for h, isotropic in ((.95, True), (.951, True), (.949, False)):
                shape = shape_moments(rectangle(w=scale, h=math.sqrt(h)*scale))
                self.assertEqual(bool(np.linalg.norm(shape.axis) == 0), isotropic)
                if isotropic:
                    np.testing.assert_array_equal(shape.orientation, np.eye(2)/2)
        np.testing.assert_allclose(shape_moments(rectangle()).orientation,
                                   shape_moments(rectangle(angle=math.pi)).orientation, atol=1e-14)

    def test_rotation_scale_and_empty_zones(self):
        zones = {1:4, 2:100}
        ref = measure_tissue([rectangle(zone=1), rectangle(x=4, zone=2)], zones, 2)
        for angle, expected in ((math.pi/6, .5), (math.pi/2, 1)):
            rotated = [rectangle(zone=1, angle=angle),
                       rectangle(x=4*math.cos(angle), y=4*math.sin(angle), zone=2, angle=angle)]
            errors = error_components(measure_tissue(rotated, zones, 2), ref, zones)
            for field in ('orientation','arrangement'):
                self.assertAlmostEqual(errors[field], expected)
        empty = measure_tissue([], zones, 2)
        err = error_components(empty, ref, zones)
        self.assertEqual(err, dict(count=1, size=1, orientation=1, arrangement=1, density=1))
        self.assertTrue(all(v == 0 for v in error_components(empty, empty, zones).values()))
        lone = measure_tissue([rectangle(zone=1)], zones, 2)
        self.assertIsNone(lone['zones'][0]['arrangement'])
        self.assertEqual(error_components(lone, empty, zones)['density'], .5)
        # Opposite aligned neighbors on either side form a single column.
        line = measure_tissue([rectangle(x=x) for x in (-3,0,3)], {1:12}, 2)
        np.testing.assert_allclose(line['zones'][0]['arrangement'], np.diag([1.,0]), atol=1e-14)

    def test_ties_across_zones_and_equal_zone_weights(self):
        zones = {1:100, 2:30}
        cells = [rectangle(), rectangle(x=3, zone=2), rectangle(y=3, zone=2)]
        for order in (cells, list(reversed(cells))):
            tissue = measure_tissue(order, zones, 2)
            np.testing.assert_allclose(tissue['zones'][0]['arrangement'], np.eye(2)/2, atol=1e-14)
        large = [rectangle(x=3*i) for i in range(10)] + [rectangle(x=33,zone=2,angle=math.pi/2)]
        swapped = [rectangle(x=3*i,angle=math.pi/2) for i in range(10)] + [rectangle(x=33,zone=2)]
        err = error_components(measure_tissue(swapped,zones,2),measure_tissue(large,zones,2),zones)
        self.assertAlmostEqual(err['orientation'],1.)
        self.assertAlmostEqual(err['arrangement'],0.)

    def test_three_percentiles_count_normalization_and_density_floor(self):
        zones = {1:50}
        a = measure_tissue([rectangle(w=i, x=i*5) for i in (1,2,7)], zones, 2)
        np.testing.assert_allclose(a['size_percentiles'], [1.2,2,6])
        b = copy.deepcopy(a)
        b['size_percentiles'] = [.6,2,12]
        self.assertAlmostEqual(error_components(b,a,zones)['size'], .5)
        for target in (0,10,100):
            a['count'], b['count'] = target, target+2
            err = error_components(b,a,zones)
            self.assertAlmostEqual(err['count'], 2/max(1,target))
        a['zones'][0]['density'], b['zones'][0]['density'] = 0., .02
        self.assertAlmostEqual(error_components(b,a,zones)['density'],1.)

    def test_hold_duplicate_coverage_and_group_weights(self):
        zones = {1:5}
        empty = measure_tissue([],zones,2)
        one = measure_tissue([rectangle()],zones,2)
        rows = [dict(time=t,measurements=m) for t,m in [(3.,one),(0.,empty),(2.,empty)]]
        scorer = TargetScorer(rows,3,zones)
        self.assertEqual(scorer(empty,2.7),0)
        self.assertEqual(scorer(one,3.),0)
        self.assertAlmostEqual(scorer(empty,3.),-(1+1+1+.5)/4)
        with self.assertRaisesRegex(ValueError,r'3.0 in rows 1 and 4'):
            TargetScorer(rows+[rows[0]],3,zones)
        with self.assertRaisesRegex(ValueError,r'2.0.*reset time 0'):
            TargetScorer([rows[0],rows[2]],3,zones)
        with self.assertRaisesRegex(ValueError,r'3.0.*limit 4'):
            TargetScorer(rows,4,zones)

    def test_duration_weighting(self):
        self.assertAlmostEqual(episode_error([(0,1,0),(1,2,0),(2,2.5,-1)]),.2)
        with self.assertRaisesRegex(ValueError,'positive scored duration'):
            episode_error([(0,0,0)])
