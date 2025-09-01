import os, yaml
from assistant.agent import A1
from data.biodatalab_data.benchmark.eval_tools.eval_interface import eval_tool_dict

benchmark_dir = 'data/biodatalab_data/benchmark/tasks'
benchmark_datasets = os.listdir(benchmark_dir)
tasks_yamls = []
for dataset in benchmark_datasets:
    tasks_yamls.extend([os.path.join(benchmark_dir, dataset, f) for f in os.listdir(os.path.join(benchmark_dir, dataset)) if f.endswith('.yaml')])

success_num = 0
for t_yaml in tasks_yamls:
    agent = A1(path='./operation_env', llm='gpt-4.1', debug=False)
    with open(t_yaml, 'r') as f:
        task_config = yaml.safe_load(f)

    print(task_config)
    agent.go(task_config['Task_description'])

    eval_func = eval_tool_dict[task_config['Evaluate_tool']]

    eval_res = eval_func(task_config['Output_data'], task_config['Ref_data'])

    if eval_res:
        print(f'task {t_yaml} successes!')
        success_num += 1
    else:
        print(f'task {t_yaml} fails!')

print(f'Success Rate: {success_num/len(tasks_yamls)}')