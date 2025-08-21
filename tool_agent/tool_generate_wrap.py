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
# Prompts (REVISED to wrap external tools)
# =========================

SYSTEM_STEP1 = """You are a senior Python engineer specializing in wrapping external Linux CLI tools for agent use.

TASK: Given a TOOL REQUIREMENT (natural language), generate ONE production-ready Python function that **wraps one external command-line tool as a single high-level interface** (i.e., one function that covers the tool's core subcommands/modes). This is NOT a fine-grained task tool—it's a general wrapper for an external program.

STRICT OUTPUT: a SINGLE JSON object with keys:
- function_name (snake_case)
- function_code (string)  # exact Python code, no markdown fences

CODING HARD RULES:
- The generated code must be pure Python using only the standard library.
- Linux-only assumptions are fine. Do not write Windows-specific code.
- The function MUST start with 'def {function_name}(' at top-level.
- The function MUST accept the **full path to the external executable** (e.g., '/path/to/trimmomatic' or '/path/to/tool/bin/tool') as a required parameter. DO NOT rely on system PATH.
- The function MUST perform **preflight validations**:
  - ensure the executable exists and is executable;
  - ensure required input files/dirs exist; ensure output directory is creatable/writable if applicable.
- One function covers the tool's core modes/subcommands via a parameter (e.g., mode='PE'|'SE' etc.). Also include an optional `extra_args: Optional[list[str]] = None` to pass through uncommon flags verbatim.
- Build the command using a clear list of strings. Use `subprocess.run` with `check=False`, capture stdout/stderr.
- The function should **return execution results only**, e.g., a dict like:
  {
    "returncode": int,
    "stdout": str,
    "stderr": str,
    "outputs": {...}   # optional, such as output file paths if determinable
  }
  Do NOT include version queries, dry-run, resource/timeout knobs, or extra logging unless absolutely required by the requirement.
- Provide a complete docstring with Args, Returns, Raises, and Examples (example calls must show the executable full path is passed explicitly).
- Add type hints for all parameters and the return type.
- Avoid any network calls.
- Be deterministic and robust with clear error messages.
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

SYSTEM_STEP2 = """You are an expert tool-registry writer.

Given a TOOL REQUIREMENT and a PRE-DEFINED function_name (the Python wrapper for an external CLI), produce a STRICT JSON object with EXACT keys:
- description: string
- name: string               # MUST equal the provided function_name
- required_parameters: [     # list of objects with required fields
    {{
      "name": "str",
      "type": "str",
      "description": "str"
    }}
  ]
- optional_parameters: [     # list of objects with required fields
    {{
      "name": "str",
      "type": "str",
      "default": <any>,      # may be str/number/bool/null
      "description": "str"
    }}
  ]

HARD RULES:
- 'name' MUST exactly equal the provided function_name.
- description MUST explicitly state **which external tool is wrapped** AND **what it does** (the concrete purpose/use).
- Output MUST be valid JSON (UTF-8), no comments, no trailing commas, no extra text.
- Only use the requirement to infer params/behavior; do not require reading the source code.
"""

USER_STEP2_TMPL = """TOOL REQUIREMENT:
{task}

BIND function_name TO:
{name}

Please output ONLY the JSON object with fields:
description, name, required_parameters, optional_parameters.
"""

# =========================
# Step 3 Prompts (REVISED to include bash installer)
# =========================

SYSTEM_STEP3 = """You are a meticulous environment and dependency analyst.

GOAL: Produce a Markdown section for the Python wrapper function **{function_name}** that:
1) Summarizes the external tool's purpose and I/O at a high level (from the given code and task context).
2) Provides a **ready-to-run Linux bash installation script** that:
   - downloads the external tool from an official or stable source URL;
   - installs it **into a subdirectory of the current working directory** (e.g., ./<tool>/<version>/bin/...);
   - ensures the **executable file resides inside that subdirectory**;
   - prints the final executable full path;
   - uses only POSIX shell commands (curl/wget, tar/unzip), no package managers.
   - Do NOT use sudo and do NOT modify system PATH; all files stay under the current directory.
3) Shows a minimal Python usage example calling the generated function with the **executable's full path** (no reliance on PATH).
4) Lists Python dependencies (usually standard library only). Keep it brief and practical.

IMPORTANT:
- Linux-only.
- Keep the script simple and robust; include basic checks and exit on error (set -euo pipefail).
- Do NOT include resource or timeout tuning.
- Return results are simple (returncode/stdout/stderr/outputs), consistent with the wrapper contract.

OUTPUT FORMAT (Markdown):
# {function_name}
<clear, actionable, concise notes plus the install script and example>
"""

USER_STEP3_TMPL = """PYTHON TOOL FUNCTION (EXACT SOURCE CODE):
{tool_code}

TASK CONTEXT (for purpose/IO summary):
{task}

Please produce the Markdown as specified.
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
    # Use official env-based initialization; base URL can be set via OPENAI_BASE_URL if needed
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
    Step 1: 生成“外部工具封装”的 Python 函数（单函数覆盖核心子命令），必须以可执行文件全路径为必参。
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
    obj = json.loads(_strip_code_fences_if_any(content))
    for key in ("function_name", "function_code"):
        if key not in obj:
            raise ValueError(f"LLM(step1) missing key: {key}")

    fn = obj["function_name"]
    raw_code = obj["function_code"]

    # Validate function_name snake_case
    if not re.match(r"^[a-z_][a-z0-9_]*$", fn):
        raise ValueError("function_name must be snake_case.")

    # Validate code parses and contains a top-level function of that name
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

    code = raw_code.rstrip() + "\n"
    return Step1Result(function_name=fn, function_code=code)

def _validate_step2_payload(obj: Dict[str, Any], bind_name: str) -> None:
    required_keys = {"description", "name", "required_parameters", "optional_parameters"}
    missing = required_keys - set(obj.keys())
    if missing:
        raise ValueError(f"LLM(step2) missing keys: {sorted(missing)}")

    if obj.get("name") != bind_name:
        raise ValueError("Description.name must equal the provided function_name.")

    if not isinstance(obj.get("description"), str) or not obj["description"].strip():
        raise ValueError("description must be a non-empty string.")

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

    op = obj.get("optional_parameters")
    if not isinstance(op, list):
        raise ValueError("optional_parameters must be a list.")
    for i, it in enumerate(op):
        if not isinstance(it, dict):
            raise ValueError(f"optional_parameters[{i}] must be an object.")
        for k in ("name", "type", "default", "description"):
            if k not in it:
                raise ValueError(f"optional_parameters[{i}] missing key: {k}")
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

def llm_step3(tool_code: str, function_name: str, task: str, model: str, client: OpenAI) -> str:
    sys_prompt = SYSTEM_STEP3.format(function_name=function_name)
    messages = [
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": USER_STEP3_TMPL.format(tool_code=tool_code, task=task)},
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
    - Atomic write via temp file + os.replace.
    - Raises ValueError if the existing file is valid JSON but not a list.
    """
    if not isinstance(new_obj, dict) or not new_obj.get("name"):
        raise ValueError("New description object must be a dict and contain a non-empty 'name' field.")

    os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)

    arr: list = []
    if os.path.exists(json_path):
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

    name = new_obj["name"]
    for i, item in enumerate(arr):
        if isinstance(item, dict) and item.get("name") == name:
            arr[i] = new_obj
            break
    else:
        arr.append(new_obj)

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
# Config & Pipeline
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
    registry_json = out.get("registry_json")
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

    # --- Step 2 ---
    s2 = llm_step2(task=task, bind_name=s1.function_name, model=m2, client=client)
    append_description_to_json_array(registry_json, s2.description)
    print(f"[Step2] Appended/updated description for `{s2.description.get('name')}` in {registry_json}")

    # --- Step 3 ---
    md = llm_step3(tool_code=s1.function_code, function_name=s1.function_name, task=task, model=m3, client=client)
    append_markdown(deps_md, md)
    print(f"[Step3] Appended dependency notes for `{s1.function_name}` to {deps_md}")

def main():
    cfg_path = os.getenv("TRI_STAGE_CONFIG", "tool_builder.yml")
    try:
        cfg = load_config(cfg_path)
        run_pipeline(cfg)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
