import yaml, json
from assistant.agent import A1
from benchmark.verifiers.eval_interface import eval_tool_dict
import argparse
import os
from joblib import Parallel, delayed, cpu_count
from tqdm import tqdm

def pmap_multi(pickleable_fn, data, n_jobs=None, verbose=1, desc=None, **kwargs):
    if n_jobs is None:
        n_jobs = cpu_count() - 1

    results = Parallel(n_jobs=n_jobs, verbose=verbose, timeout=None)(
        delayed(pickleable_fn)(*d, **kwargs) for i, d in tqdm(enumerate(data),desc=desc)
    )

    return results


def run_case(task_name, task_yaml, preddir, tmpdir, logdir, evaldir):

    try:
        with open(f'{evaldir}/{task_name}.json', 'r') as f:
            eval_res = json.load(f)
        if eval_res['cost'] > 0.0:
            print(f'{task_name} already evaluated with cost {eval_res["cost"]}, skipping...')
            return eval_res['success']
    except Exception as e:
        print(f'Failed to read existing evaluation result for {task_name}, re-evaluating... Error: {e}')
        return False


parser = argparse.ArgumentParser()
parser.add_argument('--task_yaml_dir', type=str, default='benchmark/tasks')
parser.add_argument('--llm', type=str, default='gemini-3-flash-preview')
parser.add_argument('--pred_res_dir', type=str, default='/root/autodl-tmp/biomni_v20260208_pred_results')
parser.add_argument('--temp_res_dir', type=str, default='/root/autodl-tmp/biomni_v20260208_tmp')
parser.add_argument('--log_dir', type=str, default='/root/autodl-tmp/biomni_v20260208_log')
parser.add_argument('--evaluate_result', type=str, default='/root/autodl-tmp/biomni_v20260208_evaluate_results')
args = parser.parse_args()

preddir = os.path.join(args.pred_res_dir, args.llm)
tmpdir = os.path.join(args.temp_res_dir, args.llm)
logdir = os.path.join(args.log_dir, args.llm)
evaldir = os.path.join(args.evaluate_result, args.llm)

os.makedirs(preddir, exist_ok=True)
os.makedirs(tmpdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)
os.makedirs(evaldir, exist_ok=True)

task_names = [task_name for task_name in eval_tool_dict.keys()]
task_yamls = [os.path.join(args.task_yaml_dir, f'{task_name}.yaml') for task_name in task_names]

res = pmap_multi(run_case,
                zip(task_names,task_yamls),
                preddir=preddir,
                tmpdir=tmpdir,
                logdir=logdir,
                evaldir=evaldir,
                n_jobs=8,
                desc='Evaluating tasks')
# print(res)
print(f'success rate: {sum(res)}/{len(res)}={sum(res)/len(res)}')