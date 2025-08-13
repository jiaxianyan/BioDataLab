import importlib

def read_biodatalab_eval_module2api():
    fields = [
        "eval",
    ]

    module2api = {}
    for field in fields:
        module_name = f"data.biodatalab_data.benchmark.eval_tools.{field}"
        module = importlib.import_module(module_name)
        module2api[f"assistant.tool.{field}"] = module.description
    return module2api
