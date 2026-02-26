<!-- ![logo](visualization/logo.jpg) -->
# BioDataLab
This is the offical implement of "BioDataLab: Benchmarking LLM Agents on Real-World Biological Database Curation for Data-Driven Scientific Discovery".


![logo](visualization/intro_demo.jpg)

## Quick Start

### Install the envrironment
```bash
conda create -f environment.yml
```
### LLMs API 
Setting the LLM API_KEY and BASE_URL in `assistant\llm.py`.
```python
API_KEY = ""
BASE_URL = ""
```

### Basic Usage of BioDataLab
Evaluate on one task:
```bash
conda activate biomni_e1
python3 run_evaluate_case_biomni.py --task_yaml=benchmark/tasks/fusionneoantigen_annotate_2.yaml
```
Run bash scirpt on all tasks
```bash
conda activate biomni_e1
bash evaluate_bash_scripts\run_evaluate_batch_biomni_gemini-3-flash-preview.sh
```

## Evaluation Results

| Model | Open-World SR | Open-World VR | Structured Info SR | Structured Info VR | Functional Feature SR | Functional Feature VR | Data Refinement SR | Data Refinement VR | Overall SR | Overall VR | Cost |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Proprietary Models** | | | | | | | | | | | |
| GPT-5.2 | 4.5 | 86.4 | 0.0 | 87.5 | 34.5 | 82.8 | 28.0 | 84.0 | 18.0 | 85.0 | 0.22 |
| GPT-5.1 | 0.0 | 68.2 | 0.0 | 41.7 | 10.3 | 44.8 | 12.0 | 32.0 | 6.0 | 46.0 | 0.17 |
| Gemini-3.0-Pro | 18.2 | 95.5 | 41.7 | 100.0 | 51.7 | 96.6 | 44.0 | 100.0 | 40.0 | 98.0 | 0.19 |
| Gemini-3.0-Flash | 18.2 | 72.7 | 16.7 | 75.0 | 44.8 | 89.6 | 44.0 | 80.0 | 32.0 | 80.0 | 0.04 |
| Claude-4.5-Sonnet | 18.2 | 95.5 | 33.3 | 91.7 | 41.4 | 96.6 | 40.0 | 96.0 | 34.0 | 95.0 | 1.82 |
| Claude-4.5-Haiku | 4.5 | 100.0 | 33.3 | 100.0 | 37.9 | 100.0 | 44.0 | 100.0 | 31.0 | 100.0 | 0.69 |
| Seed-1.8 | 13.6 | 100.0 | 20.8 | 100.0 | 31.0 | 96.6 | 28.0 | 96.0 | 24.0 | 98.0 | 0.02 |
| **Open Source Models** | | | | | | | | | | | |
| Kimi-K2.5 | 13.6 | 40.9 | 25.0 | 20.8 | 48.3 | 55.2 | 48.0 | 44.0 | 35.0 | 41.0 | 0.18 |
| DeepSeek-V3.2 | 4.5 | 100.0 | 16.7 | 100.0 | 51.7 | 100.0 | 44.0 | 100.0 | 31.0 | 100.0 | 0.43 |
| Qwen-3-Max | 9.1 | 95.5 | 4.2 | 100.0 | 34.5 | 93.1 | 24.0 | 96.0 | 19.0 | 96.0 | 0.08 |
| GLM-4.7 | 13.6 | 100.0 | 25.0 | 100.0 | 44.8 | 100.0 | 44.0 | 100.0 | 33.0 | 100.0 | 0.23 |