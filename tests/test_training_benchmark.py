import copy
import math
import unittest
from unittest.mock import Mock, patch

import numpy as np

from benchmarks.training import run


class TrainingBenchmarkTest(unittest.TestCase):
    def test_policy_uses_simulation_units_and_current_observation(self):
        self.assertEqual(run.feedback([0.,1.], 1., 1.),1.)
        self.assertGreater(run.feedback([0.,1.], 2., 1.),1.)
        self.assertEqual(run.feedback([-10000.,1.], 1., 1.),run.CONTROL.lower)
        self.assertEqual(run.feedback([10000.,1.], 1., 1.),run.CONTROL.upper)

    def test_optimizer_obeys_budget(self):
        for budget in (4,7,23,64):
            history = []
            def objective(weights):
                history.append(weights)
                return np.linalg.norm(weights-[.12,.75])
            run.search(objective,budget)
            self.assertEqual(len(history),budget)
            if budget == 64:
                self.assertLess(min(np.linalg.norm(w-[.12,.75]) for w in history),.0005)

    def test_missing_score_invalidates_comparison_and_seeds_count_equally(self):
        seeds = [1,2]
        records = dict(policy=[dict(seed=1,score=1.,duration=2),dict(seed=2,score=0.,duration=8)],
                       default=[dict(seed=s,score=.55) for s in seeds],
                       random=[dict(seed=s,score=.51) for s in seeds])
        summary = run.compare_methods(records,seeds)
        self.assertTrue(summary['passed'])
        self.assertEqual(summary['means']['policy'],.5)
        records['random'][0].update(score=None,failure="solver error")
        del records['default'][1]
        summary = run.compare_methods(records,seeds)
        self.assertIsNone(summary['passed'])
        self.assertIsNone(summary['means'])
        self.assertEqual(summary['missing'],[dict(method='default',seed=2,reason='missing episode'),
                                             dict(method='random',seed=1,reason='solver error')])

    def test_first_replay_mismatch_per_seed(self):
        records = [dict(time=t,observation=[1.,2.],reason='active' if t<1 else 'time_limit')
                   for t in (0.,.5,1.)]
        batch = dict(episodes=[dict(seed=s,records=copy.deepcopy(records)) for s in run.BATCH])
        other = copy.deepcopy(batch)
        episodes = other['episodes']
        episodes[0]['records'][0]['observation'][0] += 2e-5
        episodes[0]['records'][1]['observation'][0] += 2e-5
        episodes[1]['records'][1]['observation'][0] += 1e-12
        episodes[1]['records'][2]['reason'] = 'empty_tissue'
        episodes[2]['records'].pop()
        episodes[3]['records'][2]['time'] = math.nan
        differences = run.replay_mismatches(batch,other)
        self.assertEqual([d['at'] for d in differences],['reset',2,2,2])
        for mismatch,field in zip(differences,['observation','reason','missing interval','time']):
            self.assertIn(field,mismatch['field'])

    def test_startup_failure_keeps_deadline_and_closes_every_environment(self):
        envs = [Mock() for _ in run.BATCH]
        with patch.object(run,'make_env',side_effect=envs), patch.object(
                run,'start_batch', side_effect=TimeoutError('missing worker IDs 2, 3')) as start:
            result = run.batch(None,{},True)
        start.assert_called_once_with(envs,timeout=run.READY_TIMEOUT)
        for env in envs:
            env.close.assert_called_once()
        self.assertIsNone(result['rate'])
        self.assertIn('missing worker IDs 2, 3',result['startup_failure'])
        self.assertEqual(result['ready_timeout_seconds'],run.READY_TIMEOUT)
