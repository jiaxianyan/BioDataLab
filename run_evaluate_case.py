import os, yaml
from assistant.agent import A1
from data.biodatalab_data.benchmark.eval_tools.eval_interface import eval_tool_dict

agent = A1(path='./operation_env', llm='gpt-4.1', debug=False)
t_yaml = 'data/biodatalab_data/benchmark/tasks/ASMdb/3.yaml'
with open(t_yaml, 'r') as f:
    task_config = yaml.safe_load(f)

print(task_config)

if len(task_config['Input_data']) == 0:
    input_data_des = ''
else:
    input_files = ','.join(task_config['Input_data'])
    input_data_des = f'The input files to this task are: {input_files}\n'

print('============================== Task description ==============================')
print(input_data_des + task_config['Task_description'])
agent.go(input_data_des + task_config['Task_description'])

eval_func = eval_tool_dict[task_config['Evaluate_tool']]

eval_res = eval_func(task_config['Output_data'], task_config['Ref_data'])

if eval_res:
    print(f'task {t_yaml} successes!')
else:
    print(f'task {t_yaml} fails!')

