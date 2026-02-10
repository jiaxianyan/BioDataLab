import yaml, json
from assistant.agent import A1
from benchmark.verifiers.eval_interface import eval_tool_dict
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--task_yaml', type=str)
parser.add_argument('--llm', type=str, default='gpt-5.1')
parser.add_argument('--pred_res_dir', type=str, default='./pred_results_new')
parser.add_argument('--temp_res_dir', type=str, default='./tmp_new')
parser.add_argument('--log_dir', type=str, default='./log')
parser.add_argument('--evaluate_result', type=str, default='./evaluate_results')
args = parser.parse_args()

preddir = os.path.join(args.pred_res_dir, args.llm)
tmpdir = os.path.join(args.temp_res_dir, args.llm)
logdir = os.path.join(args.log_dir, args.llm)
evaldir = os.path.join(args.evaluate_result, args.llm)

os.makedirs(preddir, exist_ok=True)
os.makedirs(tmpdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)
os.makedirs(evaldir, exist_ok=True)

log_file_path = os.path.join(logdir, os.path.basename(args.task_yaml).replace('.yaml','.json'))
agent = A1(path='./operation_env', llm=args.llm, debug=True)

with open(args.task_yaml, 'r') as f:
    task_config = yaml.safe_load(f)

print(task_config)

print('============================== <Task description> ==============================')
new_task_desc = task_config['Task'].replace('<r></r>', preddir)
new_task_desc = new_task_desc.replace('<t></t>', tmpdir)

print(new_task_desc)
print('============================== </Task description> ==============================')
trajectory, u = agent.go(new_task_desc, log_dir=logdir, task_id=os.path.basename(args.task_yaml).replace('.yaml',''))

with open(log_file_path, 'w') as f:
        json.dump(trajectory, f, indent=4)

eval_func = eval_tool_dict[task_config['Verifier']]
eval_res = eval_func(d=preddir)

with open(f'{evaldir}/{os.path.basename(args.task_yaml).replace(".yaml","_eval.json")}', 'w') as f:
    json.dump({'success': eval_res,
               'cost': u}, 
               f, 
               indent=4)

print(f'Evaluation result: {eval_res}')
# if eval_res:
#     print(f'✅ task {t_yaml} successes!')
# else:
#     print(f'❌ task {t_yaml} fails!')

