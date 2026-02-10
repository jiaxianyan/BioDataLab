import yaml, json
from assistant.agent import A1
from benchmark.verifiers_valid.eval_interface import eval_tool_dict
import argparse
import os
from joblib import Parallel, delayed, cpu_count
from tqdm import tqdm
import pandas as pd

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


def check_valid(task_name, task_yaml, preddir, tmpdir, logdir, evaldir):
    
    with open(task_yaml, 'r') as f:
        task_config = yaml.safe_load(f)

    eval_func = eval_tool_dict[task_config['Verifier']]
    eval_res = eval_func(d=preddir)
    eval_res = bool(eval_res)
    
    return eval_res

parser = argparse.ArgumentParser()
parser.add_argument('--task_yaml_dir', type=str, default='benchmark/tasks')
parser.add_argument('--llm', type=str, default='doubao-seed-1-8')
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

res_success = pmap_multi(run_case,
                zip(task_names,task_yamls),
                preddir=preddir,
                tmpdir=tmpdir,
                logdir=logdir,
                evaldir=evaldir,
                n_jobs=8,
                desc='Evaluating tasks')

res_valid = pmap_multi(check_valid,
                zip(task_names,task_yamls),
                preddir=preddir,
                tmpdir=tmpdir,
                logdir=logdir,
                evaldir=evaldir,
                n_jobs=8,
                desc='Evaluating tasks')


# --- 核心统计部分 ---
    
# 1. 加载映射关系
df_meta = pd.read_csv('kdd_scripts/benchmark.csv')
# 建立 Task Name -> task_type 的映射字典
nametotype = dict(zip(df_meta['Task Name'], df_meta['task_type']))

# 2. 汇总当前运行结果
results_df = pd.DataFrame({
    'task_name': task_names,
    'is_success': res_success,
    'is_valid': res_valid
})

# 3. 映射类型并计算统计值
results_df['task_type'] = results_df['task_name'].map(nametotype).fillna('unknown')

# 分组统计：计算平均值（即率）和数量
summary = results_df.groupby('task_type').agg({
    'is_success': ['mean', 'count'],
    'is_valid': ['mean']
})

# 整理列名方便阅读
summary.columns = ['Success Rate', 'Count', 'Valid Rate']

print("\n" + "="*50)
print(f"OVERALL RESULTS ({args.llm})")
print("="*50)
print(summary)
print("-" * 50)
print(f"Total Success Rate: {results_df['is_success'].mean():.2%}")
print(f"Total Valid Rate:   {results_df['is_valid'].mean():.2%}")
print("="*50)

# 可选：保存详细结果到 CSV
# summary.to_csv(f'eval_summary_{args.llm}.csv')