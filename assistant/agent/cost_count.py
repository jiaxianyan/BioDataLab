from dataclasses import dataclass, field
from typing import Dict, Optional, Any
from langchain_core.callbacks import BaseCallbackHandler

@dataclass
class RunUsage:
    input_tokens: int = 0
    output_tokens: int = 0
    total_tokens: int = 0
    cost: float = 0.0
    missing_calls: int = 0  # 记录拿不到usage的次数
    calls: int = 0
    
def extract_usage(ai_msg):
    usage = getattr(ai_msg, "usage_metadata", None)
    if isinstance(usage, dict) and (("total_tokens" in usage) or ("input_tokens" in usage)):
        return {
            "input_tokens": usage.get("input_tokens", 0),
            "output_tokens": usage.get("output_tokens", 0),
            "total_tokens": usage.get("total_tokens", 0),
        }

    resp_meta = getattr(ai_msg, "response_metadata", None)
    if isinstance(resp_meta, dict):
        um = resp_meta.get("usageMetadata") or resp_meta.get("usage_metadata")
        if isinstance(um, dict):
            return {
                "input_tokens": um.get("promptTokenCount", 0),
                "output_tokens": um.get("candidatesTokenCount", 0),
                "total_tokens": um.get("totalTokenCount", 0),
            }
    return None

class RunUsageCallback(BaseCallbackHandler):
    def __init__(self, agent):
        self.agent = agent  # 直接拿 agent 引用，写入 agent.run_usage

    def on_llm_end(self, response, **kwargs):
        u = self.agent.run_usage
        u.calls += 1

        ai_msg = None
        try:
            ai_msg = response.generations[0][0].message
        except Exception:
            pass

        usage = extract_usage(ai_msg) if ai_msg is not None else None

        # 兜底：有些实现把 usage 放 llm_output
        if usage is None:
            llm_output = getattr(response, "llm_output", None) or {}
            um = llm_output.get("usage") or llm_output.get("usageMetadata") or llm_output.get("usage_metadata")
            if isinstance(um, dict):
                usage = {
                    "input_tokens": um.get("promptTokenCount", um.get("input_tokens", 0)),
                    "output_tokens": um.get("candidatesTokenCount", um.get("output_tokens", 0)),
                    "total_tokens": um.get("totalTokenCount", um.get("total_tokens", 0)),
                }

        if usage is None:
            u.missing_calls += 1
            return

        u.input_tokens += int(usage["input_tokens"] or 0)
        u.output_tokens += int(usage["output_tokens"] or 0)
        u.total_tokens += int(usage["total_tokens"] or 0)
        # u.cost += calc_cost(...)  # 需要你按渠道填价格表