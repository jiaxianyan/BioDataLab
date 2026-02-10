![logo](visualization/logo.jpg)
# BioDataLab
This is the offical implement of "BioDataLab: Benchmarking LLM Agents on Real-World Biological Database Curation for Data-Driven Scientific Discovery".


![logo](visualization/intro_demo.jpg)

## Quick Start

### Install the envrironment
```bash
conda create -f environment.yml
```

### Basic Usage of BioDataLab

```
# evaluate on one task
conda activate biomni_e1
python3 run_evaluate_case_biomni.py --task_yaml=benchmark/tasks/fusionneoantigen_annotate_2.yaml
```

```
# run bash scirpt on all tasks
conda activate biomni_e1
bash evaluate_bash_scripts\run_evaluate_batch_biomni_gemini-3-flash-preview.sh
```