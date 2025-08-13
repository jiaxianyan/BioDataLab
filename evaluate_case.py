import os, yaml
from assistant.agent import A1
from data.biodatalab_data.benchmark.eval_tools.eval import eval_1

agent = A1(path='./data', llm='gpt-4.1', debug=False)
task_dir = 'data/biodatalab_data/benchmark/tasks'
benchmark_task_files = os.listdir(task_dir)
task = '1.yaml'
# for task in benchmark_task_files:

with open(os.path.join(task_dir, task), 'r') as f:
    task_config = yaml.safe_load(f)
print(task_config)
agent.go(task_config['Task_description'])
eval_res = eval_1(task_config['Output_data'], task_config['Ref_data'])
if eval_res:
    print('task 1 successes!')
else:
    print('task 1 fails!')
