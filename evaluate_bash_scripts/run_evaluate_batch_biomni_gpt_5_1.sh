LLM_NAME="gpt-5.1"
LOG_DIR="/root/autodl-tmp/biomni_v20260208_bash_log/${LLM_NAME}"

mkdir -p "${LOG_DIR}"

python - <<'PY'
from benchmark.verifiers.eval_interface import eval_tool_dict
for k in eval_tool_dict.keys():
    print(k)
PY

python - <<'PY' | parallel --jobs 16 --bar \
"python run_evaluate_batch_biomni_for_bash.py \
  --specific_task {} \
  --llm ${LLM_NAME} \
  > ${LOG_DIR}/{}.log 2>&1"
from benchmark.verifiers.eval_interface import eval_tool_dict
for k in eval_tool_dict.keys():
    print(k)
PY

python run_only_evaluate.py --llm ${LLM_NAME}
