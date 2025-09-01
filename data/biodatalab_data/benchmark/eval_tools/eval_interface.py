# import importlib

# def read_biodatalab_eval_module2api():
#     fields = [
#         "eval",
#     ]

#     module2api = {}
#     for field in fields:
#         module_name = f"data.biodatalab_data.benchmark.eval_tools.{field}"
#         module = importlib.import_module(module_name)
#         module2api[f"assistant.tool.{field}"] = module.description
#     return module2api

from .common_eval import eval_list
from .AlphaFold_Protein_Structure_Database.eval import eval_AlphaFold_Protein_Structure_Database_fasta_file
from .ASMdb.eval import eval_ASMdb_fastq

eval_tool_dict = {'eval_list':eval_list,
                  'eval_AlphaFold_Protein_Structure_Database_fasta_file':eval_AlphaFold_Protein_Structure_Database_fasta_file,
                  'eval_ASMdb_fastq':eval_ASMdb_fastq}