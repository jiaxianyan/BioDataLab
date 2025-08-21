# -*- coding: utf-8 -*-
# file: tri_stage_tool_builder_yaml.py
# Run: python tri_stage_tool_builder_yaml.py
# Reqs: pip install --upgrade openai pyyaml
# Env : export OPENAI_API_KEY=sk-...

import json
import os
import re
import sys
import textwrap
from dataclasses import dataclass
from typing import Any, Dict, Optional
from pprint import pformat
import ast

try:
    import yaml  # PyYAML
except Exception as e:
    raise RuntimeError("Please `pip install pyyaml` to load the YAML config.") from e

try:
    from openai import OpenAI
except Exception as e:
    raise RuntimeError(
        "Please `pip install --upgrade openai`. "
        "Docs: https://platform.openai.com/docs"
    ) from e


API_SECRET_KEY = "sk-di6EQiNDxyGzzoFa93389fA755B14eAfAf31B863F2E6C369"
BASE_URL = "https://aihubmix.com/v1"
os.environ["OPENAI_API_KEY"] = API_SECRET_KEY
os.environ["OPENAI_BASE_URL"] = BASE_URL
os.environ["GOOGLE_API_KEY"] = API_SECRET_KEY
os.environ["GOOGLE_BASE_URL"] = BASE_URL

# =========================
# Prompts
# =========================

SYSTEM_STEP1 = """You are an expert Python toolsmith for bioinformatics.
Task: Given a TOOL REQUIREMENT (natural language), produce ONE production-ready Python function.

STRICT OUTPUT: a SINGLE JSON object with keys:
- function_name (snake_case)
- function_code (string)  # exact Python code, no markdown fences

CODING RULES:
- function_code MUST start with 'def {function_name}(' at top-level.
- Use only standard library unless the requirement explicitly needs more.
- Provide a complete docstring (Args, Returns, Raises, Examples).
- Robust error handling and deterministic behavior.
- If external CLIs/programs are needed (e.g., samtools, annovar, bcftools):
  - DO NOT rely on system PATH.
  - The function MUST accept their executable paths or root dir as arguments (e.g., annovar_dir, samtools_path)
    OR define explicit path params with defaults as None, then validate them.
  - Add preflight checks with clear error messages when paths are missing/invalid.
- Add type hints for all parameters and return type.
- Avoid network calls unless explicitly required by the TOOL REQUIREMENT.
- No extra commentary; ONLY the JSON object.
"""

USER_STEP1_TMPL = """TOOL REQUIREMENT:
{task}

OUTPUT FORMAT:
{{
  "function_name": "snake_case_name",
  "function_code": "def snake_case_name(...):\\n    ..."
}}
"""

# =========================
# Step 2 Prompts (REVISED)
# =========================

SYSTEM_STEP2 = """You are an expert registry writer.
Given a TOOL REQUIREMENT and a PRE-DEFINED function_name, produce a STRICT JSON object.

STRICT OUTPUT: a SINGLE JSON object with EXACT keys:
- description: string
- name: string               # MUST equal the provided function_name
- required_parameters: [     # list of objects with required fields
    {
      "name": "str",
      "type": "str",
      "description": "str"
    }
  ]
- optional_parameters: [     # list of objects with required fields
    {
      "name": "str",
      "type": "str",
      "default": <any>,      # may be str/number/bool/null
      "description": "str"
    }
  ]

HARD RULES:
- 'name' MUST exactly equal the provided function_name.
- ONLY use the requirement to infer params/behavior; do not require reading the source code.
- Output MUST be valid JSON (UTF-8), no comments, no trailing commas, no extra text.
"""

USER_STEP2_TMPL = """TOOL REQUIREMENT:
{task}

BIND function_name TO:
{name}

Please output ONLY the JSON object with fields:
description, name, required_parameters, optional_parameters.
"""


SYSTEM_STEP3 = """You are a meticulous environment and dependency analyst.
Given a Python tool function (exact source code), list all environment needs:

- Python library dependencies (pip/conda names, minimal versions if applicable)
- External programs/CLIs (names, minimal versions, typical install methods)
- OS-level prerequisites if relevant
- Suggested environment setup (venv/conda), example install commands
- Notes about paths, permissions, large-file handling, performance tips

OUTPUT FORMAT (Markdown):
# {function_name}
<clear, actionable dependency & environment notes>

Be explicit and complete. No need to be overly formal—clarity first.
"""

USER_STEP3_TMPL = """PYTHON TOOL FUNCTION:
{tool_code}

Please produce the Markdown section described above.
"""


# =========================
# Data structures
# =========================

@dataclass
class Step1Result:
    function_name: str
    function_code: str


@dataclass
class Step2Result:
    description: Dict[str, Any]


# =========================
# LLM helpers
# =========================

def _client(api_key: Optional[str] = None) -> OpenAI:
    if api_key:
        return OpenAI(api_key=api_key)
    return OpenAI()

def _strip_code_fences_if_any(s: str) -> str:
    s2 = s.strip()
    if s2.startswith("```"):
        s2 = re.sub(r"^```[a-zA-Z0-9_+-]*\s*", "", s2)
        s2 = re.sub(r"\s*```$", "", s2)
    return s2

def llm_step1(task: str, model: str, client: OpenAI) -> Step1Result:
    """
    Step 1: 根据任务描述生成 Python 工具函数。
    允许在函数前出现 import / from ... import 等顶层语句，但必须包含一个顶层
    def {function_name}(...): 的函数定义。
    """
    messages = [
        {"role": "system", "content": SYSTEM_STEP1},
        {"role": "user", "content": USER_STEP1_TMPL.format(task=task)},
    ]
    resp = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=0,
        response_format={"type": "json_object"},
    )
    content = resp.choices[0].message.content
    # print(content)  # Debug: print the raw response content

    # 解析 JSON（容忍意外的代码围栏）
    obj = json.loads(_strip_code_fences_if_any(content))
    for key in ("function_name", "function_code"):
        if key not in obj:
            raise ValueError(f"LLM(step1) missing key: {key}")

    fn = obj["function_name"]
    raw_code = obj["function_code"]

    # 校验函数名 snake_case
    if not re.match(r"^[a-z_][a-z0-9_]*$", fn):
        raise ValueError("function_name must be snake_case.")

    # 用 AST 解析，确保代码语法正确，且存在一个顶层函数定义与 fn 同名
    try:
        mod = ast.parse(raw_code, filename="<function_code>", mode="exec")
    except SyntaxError as e:
        raise ValueError(f"function_code is not valid Python: {e}") from e

    found_top_level_def = False
    for node in mod.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == fn:
            found_top_level_def = True
            break

    if not found_top_level_def:
        raise ValueError(f"function_code must define a top-level function named `{fn}`")

    # 允许前置 import/常量/注释等；只需把代码整理成以换行结尾
    code = raw_code.rstrip() + "\n"

    return Step1Result(function_name=fn, function_code=code)


def _validate_step2_payload(obj: Dict[str, Any], bind_name: str) -> None:
    # Required top-level keys
    required_keys = {"description", "name", "required_parameters", "optional_parameters"}
    missing = required_keys - set(obj.keys())
    if missing:
        raise ValueError(f"LLM(step2) missing keys: {sorted(missing)}")

    # name must match
    if obj.get("name") != bind_name:
        raise ValueError("Description.name must equal the provided function_name.")

    # description: non-empty string
    if not isinstance(obj.get("description"), str) or not obj["description"].strip():
        raise ValueError("description must be a non-empty string.")

    # required_parameters: list of dicts with {name,type,description}
    rp = obj.get("required_parameters")
    if not isinstance(rp, list):
        raise ValueError("required_parameters must be a list.")
    for i, it in enumerate(rp):
        if not isinstance(it, dict):
            raise ValueError(f"required_parameters[{i}] must be an object.")
        for k in ("name", "type", "description"):
            if k not in it:
                raise ValueError(f"required_parameters[{i}] missing key: {k}")
            if not isinstance(it[k], str) or not it[k].strip():
                raise ValueError(f"required_parameters[{i}].{k} must be a non-empty string.")

    # optional_parameters: list of dicts with {name,type,default,description}
    op = obj.get("optional_parameters")
    if not isinstance(op, list):
        raise ValueError("optional_parameters must be a list.")
    for i, it in enumerate(op):
        if not isinstance(it, dict):
            raise ValueError(f"optional_parameters[{i}] must be an object.")
        for k in ("name", "type", "default", "description"):
            if k not in it:
                raise ValueError(f"optional_parameters[{i}] missing key: {k}")
        # name/type/description require non-empty strings; default can be any JSON type
        for k in ("name", "type", "description"):
            if not isinstance(it[k], str) or not it[k].strip():
                raise ValueError(f"optional_parameters[{i}].{k} must be a non-empty string.")


def llm_step2(task: str, bind_name: str, model: str, client: OpenAI) -> Step2Result:
    messages = [
        {"role": "system", "content": SYSTEM_STEP2},
        {"role": "user", "content": USER_STEP2_TMPL.format(task=task, name=bind_name)},
    ]
    resp = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=0,
        response_format={"type": "json_object"},
    )
    content = resp.choices[0].message.content
    obj = json.loads(_strip_code_fences_if_any(content))
    _validate_step2_payload(obj, bind_name)
    return Step2Result(description=obj)


def llm_step3(tool_code: str, function_name: str, model: str, client: OpenAI) -> str:
    sys_prompt = SYSTEM_STEP3.format(function_name=function_name)
    messages = [
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": USER_STEP3_TMPL.format(tool_code=tool_code)},
    ]
    resp = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=0,
    )
    md = resp.choices[0].message.content.strip()
    # Ensure header
    first = md.splitlines()[0].strip() if md else ""
    header = f"# {function_name}"
    if first != header:
        md = header + "\n\n" + md
    return md


# =========================
# File ops
# =========================

def append_description_to_json_array(json_path: str, new_obj: Dict[str, Any]) -> None:
    """
    Append (or upsert) a tool description object into a JSON file containing an array.

    Guarantees:
    - If file is missing, or exists but is 0-byte / whitespace-only, treat as [].
    - If file contains a JSON array, append; if an item with the same "name" exists, replace it.
    - Atomic write via temp file + os.replace, so partial writes won't corrupt the file.
    - Raises ValueError if the existing file is valid JSON but not a list.

    Expected object shape (validated elsewhere):
      { "name": str, "description": str, "required_parameters": [...], "optional_parameters": [...] }
    """
    if not isinstance(new_obj, dict) or not new_obj.get("name"):
        raise ValueError("New description object must be a dict and contain a non-empty 'name' field.")

    # Ensure parent dir
    os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)

    arr: list = []
    if os.path.exists(json_path):
        # Handle empty/whitespace-only file
        if os.path.getsize(json_path) > 0:
            with open(json_path, "r", encoding="utf-8") as f:
                raw = f.read()
            if raw.strip():
                try:
                    existing = json.loads(raw)
                except json.JSONDecodeError as e:
                    raise ValueError(f"Invalid JSON in {json_path}: {e}")
                if existing is None:
                    arr = []
                elif isinstance(existing, list):
                    arr = existing
                else:
                    raise ValueError(f"{json_path} must contain a JSON array, got {type(existing).__name__}.")
        # else: zero-byte → keep arr = []

    # Upsert by name
    name = new_obj["name"]
    for i, item in enumerate(arr):
        if isinstance(item, dict) and item.get("name") == name:
            arr[i] = new_obj
            break
    else:
        arr.append(new_obj)

    # Atomic write
    tmp = json_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(arr, f, ensure_ascii=False, indent=2)
        f.write("\n")
    os.replace(tmp, json_path)


def append_function_to_py(py_path: str, function_code: str) -> None:
    if not os.path.exists(py_path):
        raise FileNotFoundError(f"Python file not found: {py_path}")
    with open(py_path, "r", encoding="utf-8") as f:
        text = f.read()
    if not text.endswith("\n"):
        text += "\n"
    if not text.endswith("\n\n"):
        text += "\n"
    text += function_code
    if not text.endswith("\n"):
        text += "\n"
    with open(py_path, "w", encoding="utf-8") as f:
        f.write(text)

def save_function_source(out_path: str, function_code: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(function_code)

def append_markdown(md_path: str, section_md: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(md_path)), exist_ok=True)
    exists = os.path.exists(md_path)
    with open(md_path, "a", encoding="utf-8") as f:
        if exists and os.path.getsize(md_path) > 0:
            f.write("\n\n")
        f.write(section_md.strip() + "\n")

def append_description_to_python_list(py_path: str, var_name: str, obj: Dict[str, Any]) -> None:
    if not os.path.exists(py_path):
        raise FileNotFoundError(f"Python file not found: {py_path}")
    with open(py_path, "r", encoding="utf-8") as f:
        code = f.read()

    m = re.search(rf"(^\s*{re.escape(var_name)}\s*=\s*\[)", code, flags=re.MULTILINE)
    if not m:
        raise ValueError(f"Variable `{var_name}` list assignment not found in {py_path} (expect '{var_name} = [').")
    start_list_idx = m.end()

    # find matching ]
    bracket = 1
    i = start_list_idx
    end_list_idx = None
    while i < len(code) and bracket > 0:
        ch = code[i]
        if ch == "[":
            bracket += 1
        elif ch == "]":
            bracket -= 1
            if bracket == 0:
                end_list_idx = i
                break
        i += 1
    if end_list_idx is None:
        raise ValueError("Unmatched '[' ... ']' in registry list.")

    assign_line_start = code.rfind("\n", 0, m.start()) + 1
    line = code[assign_line_start:m.start()]
    base_indent = re.match(r"^\s*", line).group(0)
    elem_indent = base_indent + "    "

    list_content = code[start_list_idx:end_list_idx].strip()
    py_literal = pformat(obj, width=100, compact=False)

    if not list_content:
        insertion = "\n" + textwrap.indent(py_literal, elem_indent) + "\n" + base_indent
    else:
        prev = code[:end_list_idx].rstrip()
        if prev.endswith(","):
            insertion_core = textwrap.indent(py_literal, elem_indent)
        else:
            insertion_core = "," + "\n" + textwrap.indent(py_literal, elem_indent)
        insertion = "\n" + insertion_core + "\n" + base_indent

    new_code = code[:start_list_idx] + insertion + code[end_list_idx:]
    with open(py_path, "w", encoding="utf-8") as f:
        f.write(new_code)


# =========================
# Pipeline runner (one-shot)
# =========================

def load_config(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"YAML config not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Top-level YAML must be a mapping.")
    return cfg

def run_pipeline(cfg: Dict[str, Any]) -> None:
    task = cfg.get("task")
    if not task or not isinstance(task, str):
        raise ValueError("Config must have a string field: task")

    model = cfg.get("model", {}) or {}
    m1 = model.get("step1", "gpt-5")
    m2 = model.get("step2", m1)
    m3 = model.get("step3", m1)

    out = cfg.get("output", {}) or {}
    py_file = out.get("py_file")
    out_tool_src = out.get("out_tool_src")
    registry_json = out.get("registry_json")   # <—— 新字段
    deps_md = out.get("deps_md")
    if not all([py_file, registry_json, deps_md]):
        raise ValueError("output.py_file, output.registry_json, output.deps_md are required.")

    api_key = cfg.get("openai_api_key") or os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not set (env var or `openai_api_key` in YAML).")
    client = _client(api_key)

    # --- Step 1 ---
    s1 = llm_step1(task=task, model=m1, client=client)
    append_function_to_py(py_file, s1.function_code)
    if out_tool_src:
        save_function_source(out_tool_src, s1.function_code)
    print(f"[Step1] Appended `{s1.function_name}` to {py_file}")

    # --- Step 2 (write to JSON file) ---
    s2 = llm_step2(task=task, bind_name=s1.function_name, model=m2, client=client)
    append_description_to_json_array(registry_json, s2.description)
    print(f"[Step2] Appended/updated description for `{s2.description.get('name')}` in {registry_json}")

    # --- Step 3 ---
    md = llm_step3(tool_code=s1.function_code, function_name=s1.function_name, model=m3, client=client)
    append_markdown(deps_md, md)
    print(f"[Step3] Appended dependency notes for `{s1.function_name}` to {deps_md}")


def main():
    # 不使用 argparse；从环境变量或默认文件读取配置
    cfg_path = os.getenv("TRI_STAGE_CONFIG", "tool_builder.yml")
    try:
        cfg = load_config(cfg_path)
        run_pipeline(cfg)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
