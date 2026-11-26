"""Run the predeclared 3D tissue learning or throughput qualification."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np

from vertex_rl import Control, Experiment, TissueEnv
from vertex_rl.scoring import TargetScorer, episode_error
from vertex_rl.transport import start_batch

TRAIN = (101, 113, 127)
HELD = (211, 223, 227, 229, 233)
BATCH = (307, 311, 313, 317)
CANDIDATES = 64
REFERENCE_RATE = .73
DEFAULT_RATE = .15
INTERVAL = .125
LIMIT = 1.
READY_TIMEOUT = 30.
CONTROL = Control(0, "growth_rate", 0., 2.)
MEASURES = {0: 512.}


def write_json(path, obj):
    with path.open("x") as file:
        json.dump(obj, file, indent=2, allow_nan=False)
        file.write("\n")


def geometry(seed, max_step):
    rng = np.random.Generator(np.random.PCG64(seed))
    cells = []
    for x in (-1.5,1.5):
        center = np.array([x,0,0])
        half = np.array([1.3,.85,.7]) * rng.uniform(.85,1.15,3) / 2
        cells.append({"mesh":{"shape":"box","min":(center-half).tolist(),"max":(center+half).tolist()}})
    return {"domain":{"shape":"box","min":[-4]*3,"max":[4]*3}, "cells":cells,
            "mechanics":{"volume_stiffness":2., "surface_tension":.01, "max_step":max_step}}


def make_env(binary, seed, max_step, scorer=None):
    experiment = Experiment(dimension=3, geometry=geometry(seed,max_step), binary=binary,
                            zones=[dict(id=0, growth_rate=DEFAULT_RATE, target_size=10.,
                                        allow_division=False, allow_death=False)],
                            controls=[CONTROL], zone_measures=MEASURES,
                            decision_interval=INTERVAL, time_limit=LIMIT, threads=1)
    return TissueEnv(experiment, scorer=scorer)


def references(binary, seeds):
    result = {}
    for seed in seeds:
        env = make_env(binary, seed, .02)
        try:
            _, info = env.reset(seed=seed)
            rows = [dict(time=info['time'], measurements=env.measurements)]
            while True:
                _, _, terminated, truncated, info = env.step([REFERENCE_RATE])
                rows.append(dict(time=info['time'], measurements=env.measurements))
                if terminated or truncated:
                    break
            result[seed] = dict(geometry=geometry(seed,.02), rows=rows)
        finally:
            env.close()
    return result


def sigmoid_bounded(logit, control):
    probability = math.exp(-float(np.logaddexp(0., -logit)))
    return control.lower + (control.upper-control.lower)*probability


def feedback(weights, current, initial):
    bias, gain = weights
    return sigmoid_bounded(bias + gain*math.log(current / initial), CONTROL)


def scored_episode(env, seed, *, weights=None, rate=None):
    steps = []
    try:
        _, info = env.reset(seed=seed)
        start = info['time']
        initial = env.measurements['size_percentiles'][1]
        while True:
            current = env.measurements['size_percentiles'][1]
            value = float(rate) if weights is None else feedback(weights, current, initial)
            _, reward, terminated, truncated, info = env.step([value])
            end = info['time']
            steps.append(dict(start=start, end=end, action=[value], reward=float(reward)))
            if terminated or truncated:
                score = episode_error((s['start'],s['end'],s['reward']) for s in steps)
                return dict(seed=seed, score=score, end_reason=info['reason'],
                            steps=steps)
            start = end
    except Exception as error:
        # Earlier scored intervals are preserved for diagnosis. A worker or
        # solver failure invalidates the entire method/seed episode.
        return dict(seed=seed, score=None, failure=f"{type(error).__name__}: {error}", steps=steps)


def search(objective, budget):
    """Derivative-free direct policy search on the mean training-episode error."""
    n = 0

    def evaluate(weights):
        nonlocal n
        n += 1
        return objective(weights), weights

    simplex = [evaluate(np.array(x)) for x in [(0.,.75),(.5,.75),(0.,1.25)]]
    while n < budget:
        simplex.sort(key=lambda v:v[0])
        best, second, worst = simplex
        midpoint = (best[1]+second[1])/2
        reflected = evaluate(2*midpoint-worst[1])
        if n == budget:
            break
        if reflected[0] < best[0]:
            expanded = evaluate(midpoint+2*(reflected[1]-midpoint))
            simplex[-1] = min((reflected,expanded),key=lambda v:v[0])
        elif reflected[0] < second[0]:
            simplex[-1] = reflected
        else:
            better = min((reflected,worst),key=lambda v:v[0])
            contracted = evaluate(midpoint+.5*(better[1]-midpoint))
            if contracted[0] < better[0]:
                simplex[-1] = contracted
            else:
                for j in (1,2):
                    if n == budget:
                        break
                    simplex[j] = evaluate((best[1]+simplex[j][1])/2)


def evaluated_candidate(envs, parameters):
    episodes = [scored_episode(env, seed, **parameters) for seed,env in envs.items()]
    valid = all(e['score'] is not None for e in episodes)
    mean = math.fsum(e['score'] for e in episodes)/len(episodes) if valid else None
    return dict(parameters=parameters, scores=[{k:v for k,v in e.items() if k != 'steps'} for e in episodes],
                mean=mean)


def select(binary, dest):
    refs = references(binary, TRAIN)
    write_json(dest/'training-targets.json', refs)
    history = dict(policy=[],random=[])
    with ExitStack() as stack:
        envs = {}
        for seed in TRAIN:
            scorer = TargetScorer(refs[seed]['rows'], LIMIT, MEASURES)
            envs[seed] = make_env(binary,seed,.02,scorer)
            stack.callback(envs[seed].close)
        def objective(weights):
            row = evaluated_candidate(envs,dict(weights=weights.tolist()))
            history['policy'].append(row)
            return row['mean'] if row['mean'] is not None else math.inf
        search(objective, CANDIDATES)
        rng = np.random.Generator(np.random.PCG64(10517))
        for value in rng.uniform(CONTROL.lower, CONTROL.upper, CANDIDATES):
            history['random'].append(evaluated_candidate(envs,dict(rate=float(value))))
    write_json(dest/'selection.json',history)
    chosen = {}
    for method, rows in history.items():
        valid = [r for r in rows if r['mean'] is not None]
        chosen[method] = min(valid,key=lambda r:r['mean'])['parameters'] if valid else None
    frozen = {**chosen, 'default':dict(rate=DEFAULT_RATE)}
    write_json(dest/'frozen.json',frozen)
    return frozen


def compare_methods(records, seeds):
    misses, means = [], {}
    for method in ('policy','default','random'):
        by_seed = {r['seed']:r for r in records[method]}
        for seed in seeds:
            r = by_seed.get(seed, {})
            if r.get('score') is None:
                misses.append(dict(method=method,seed=seed,reason=r.get('failure','missing episode')))
        if all(by_seed.get(s,{}).get('score') is not None for s in seeds):
            means[method] = math.fsum(by_seed[s]['score'] for s in seeds)/len(seeds)
    if misses:
        return dict(complete=False, missing=misses, means=None, passed=None)
    return dict(complete=True, missing=[], means=means,
                passed=means['policy'] < min(means['random'],means['default']))


def evaluate_heldout(binary, frozen, dest):
    records = {method:[] for method in frozen}
    refs = {}
    for seed in HELD:
        try:
            ref = references(binary,[seed])[seed]
            refs[seed] = ref
            scorer = TargetScorer(ref['rows'],LIMIT,MEASURES)
        except Exception as error:
            for episodes in records.values():
                episodes.append(dict(seed=seed,score=None,failure=f"reference failed: {error}"))
            continue
        for method, parameters in frozen.items():
            if parameters is None:
                records[method].append(dict(seed=seed,score=None,failure="no completed training candidate"))
                continue
            env = make_env(binary,seed,.02,scorer)
            try:
                records[method].append(scored_episode(env,seed,**parameters))
            finally:
                env.close()
    write_json(dest/'heldout-targets.json',refs)
    summary = compare_methods(records,HELD)
    write_json(dest/'evaluation.json',dict(episodes=records,**summary))
    return summary


def serializable(obj):
    if isinstance(obj,np.ndarray):
        return obj.tolist()
    if isinstance(obj,dict):
        return {k:serializable(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [serializable(v) for v in obj]
    if isinstance(obj,np.generic):
        return obj.item()
    return obj


def first_difference(a, b, field):
    if isinstance(a,dict) and isinstance(b,dict):
        if a.keys()!=b.keys():
            return f"{field} keys differ: {sorted(a.keys() ^ b.keys())}"
        for key in a:
            difference=first_difference(a[key],b[key],f"{field}.{key}")
            if difference:
                return difference
    elif isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b):
            return f"{field} length {len(a)} != {len(b)}"
        for i,(x,y) in enumerate(zip(a,b)):
            difference=first_difference(x,y,f"{field}[{i}]")
            if difference:
                return difference
    elif type(a) is float and type(b) is float:
        if not (math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-9)):
            return f"{field}: {a} != {b}"
    elif type(a)!=type(b) or a!=b:
        return f"{field}: {a!r} != {b!r}"
    return None


def replay(env, seed, actions):
    records = []
    start = time.perf_counter()
    finished = start
    try:
        observation,info = env.reset(seed=seed)
        finished = time.perf_counter()
        records.append(dict(time=info['time'],observation=observation,measurements=env.measurements))
        for action in actions:
            observation,reward,terminated,truncated,info = env.step(action)
            finished = time.perf_counter()
            reason = info['reason'] if terminated or truncated else 'active'
            records.append(dict(time=info['time'],observation=observation,measurements=env.measurements,
                                reward=reward,reason=reason))
            if terminated or truncated:
                return dict(seed=seed,start=start,finish=finished,records=records)
        raise RuntimeError("recorded actions ended before episode termination or truncation")
    except Exception as error:
        return dict(seed=seed,start=start,finish=time.perf_counter(),records=records,
                    failure=f"{type(error).__name__}: {error}")


def batch(binary, actions, concurrent):
    result = dict(ready_timeout_seconds=READY_TIMEOUT, active_workers=len(BATCH) if concurrent else 1,
                  ready_workers=0, threads_per_worker=1, available_cpus=len(os.sched_getaffinity(0)),
                  episodes=[], seconds=None, simulated_duration=None, rate=None)
    with ExitStack() as stack:
        envs = []
        for seed in BATCH:
            env = make_env(binary,seed,.0005)
            stack.callback(env.close)
            envs.append(env)
        try:
            start_batch(envs,timeout=READY_TIMEOUT)
        except Exception as error:
            result['startup_failure'] = f"{type(error).__name__}: {error}"
            return result
        result['ready_workers'] = len(envs)
        args = [(env,seed,actions[seed]) for env,seed in zip(envs,BATCH)]
        if concurrent:
            with ThreadPoolExecutor(max_workers=len(args)) as pool:
                futures = [pool.submit(replay,*arg) for arg in args]
                episodes = [f.result() for f in futures]
        else:
            episodes = [replay(*arg) for arg in args]
        result['seconds'] = max(e['finish'] for e in episodes)-min(e['start'] for e in episodes)
    for e in episodes:
        records = e['records']
        e['simulated_duration'] = (records[-1]['time']-records[0]['time']) if records else 0.
    result['simulated_duration'] = math.fsum(e['simulated_duration'] for e in episodes)
    if not any('failure' in e for e in episodes):
        result['rate'] = result['simulated_duration']/result['seconds']
    result['episodes'] = serializable(episodes)
    return result


def replay_mismatches(sequential, concurrent):
    mismatches = []
    seq = {r['seed']:r for r in sequential['episodes']}
    con = {r['seed']:r for r in concurrent['episodes']}
    for seed in BATCH:
        a,b=seq.get(seed),con.get(seed)
        if a is None or b is None:
            mismatches.append(dict(seed=seed,at="reset",field="missing batch episode"))
            continue
        for i in range(max(len(a['records']),len(b['records']))):
            if i>=len(a['records']) or i>=len(b['records']):
                issue="missing interval result"
            else:
                issue=first_difference(a['records'][i],b['records'][i],"result")
            if issue:
                mismatches.append(dict(seed=seed,at="reset" if i==0 else i,field=issue))
                break
        else:
            if 'failure' in a or 'failure' in b:
                mismatches.append(dict(seed=seed,at=len(a['records']),field=a.get('failure') or b['failure']))
    return mismatches


def throughput(binary, dest):
    actions = {s:np.random.Generator(np.random.PCG64(50000+s)).uniform(0.,2.,(8,1)).tolist() for s in BATCH}
    write_json(dest/'actions.json',dict(controls=[CONTROL.identifier],decision_interval=INTERVAL,
                                        time_limit=LIMIT,actions=actions))
    sequential = batch(binary,actions,False)
    concurrent = batch(binary,actions,True)
    mismatches = replay_mismatches(sequential,concurrent)
    ratio = None
    if not mismatches and sequential['rate'] and concurrent['rate']:
        ratio = concurrent['rate']/sequential['rate']
    write_json(dest/'throughput.json',dict(sequential=sequential,concurrent=concurrent,
                                           mismatches=mismatches,ratio=ratio))
    return dict(sequential_rate=sequential['rate'],concurrent_rate=concurrent['rate'],
                ratio=ratio,mismatches=mismatches)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['learn','throughput'])
    parser.add_argument('destination',type=Path,help='new directory for raw evidence')
    parser.add_argument('--binary',type=Path,required=True,help='3D worker executable')
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    dest = args.destination
    dest.mkdir(parents=True,exist_ok=False)
    write_json(dest/'machine.json',dict(host=platform.node(),platform=platform.platform(),
                                        cpu_affinity=sorted(os.sched_getaffinity(0)),
                                        executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                                        python=sys.version,numpy=np.__version__))
    if args.mode=='learn':
        frozen = select(binary,dest)
        summary = evaluate_heldout(binary,frozen,dest)
        passed = summary['passed']
    else:
        summary = throughput(binary,dest)
        passed = summary['ratio'] is not None and summary['ratio']>1.
    print(json.dumps(summary,indent=2))
    return 0 if passed else 1


if __name__=='__main__':
    sys.exit(main())
