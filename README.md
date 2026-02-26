<!-- ![logo](visualization/logo.jpg) -->
# BioDataLab
This is the offical implement of "BioDataLab: Benchmarking LLM Agents on Real-World Biological Database Curation for Data-Driven Scientific Discovery". If you encounter any issues, please reach out to jiaxianyan@mail.ustc.edu.cn.

## Introduction

![logo](visualization/intro_demo.jpg)

### Benchmark Overview
We introduce BioDataLab, a rigorous benchmark comprising 100 tasks meticulously derived from 57 high-impact database publications, covering 9 biological domains and 7 data modalities. BioDataLab evaluates the capability of autonomous agents to transform raw, heterogeneous biological resources into structured, analysis-ready databases.
Tasks are classfied by their primary intention into 4 types: open-world data retrieval, structured data extraction, functional feature annotation, and data refinement and integration. 
The following are 4 reprensentative examples of each type:
![logo](visualization/overview.jpg)

### Benchmark Access
A overall summary of all 100 tasks are listed in [BioDataLab](BioDataLab.csv).

First, download the necessary input from [huggingface]() and unzip it as `./benchmark/dataset`.

Then, the whole benchmark directory structure should look like this:

Each task 
```
|-- BioDataLab/
|---- BioDataLab.csv
|---- benchmark/
|------ datasets/           # subdirectory for input data
|-------- ...
|------ verifiers/          # subdirectory for succsess rate verfiers 
|-------- ...
|------ verifiers_valid/    # subdirectory for valid rate verfiers 
|-------- ...
|------ tasks/              # subdirectory for detailed task description yaml files 
|-------- ...
|------ gold_programs/      # subdirectory for programs used to generate groundtruth results if applicable 
|-------- ...
|------ gold_results/       # subdirectory for groundtruth
|-------- ...

```


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

If you want to evaluate on one task, for example [fusionneoantigen_annotate_2](benchmark\tasks\fusionneoantigen_annotate_2.yaml), you can run:
```bash
conda activate biomni_e1
python3 run_evaluate_case_biomni.py --task_yaml=benchmark/tasks/fusionneoantigen_annotate_2.yaml
```

If you want to evaluate all tasks, we also provided the scripts in the [directory](evaluate_bash_scripts). For example, if you want to evaluate the gemini-3-flash-preview model:
```bash
conda activate biomni_e1
bash evaluate_bash_scripts\run_evaluate_batch_biomni_gemini-3-flash-preview.sh
```

## Evaluation Results
![logo](visualization/result.jpg)

## Contributors
**Student Contributors**: Jiaxian Yan, Xi Fang, Chenmin Wu, Jintao Zhu, Yuhang Yang, [Zaixi Zhang](https://zaixizhang.github.io/), Meijing Fang, and Chenxi Du

**Supervisors**: [Qi Liu](http://staff.ustc.edu.cn/~qiliuql/), [Kai Zhang](http://home.ustc.edu.cn/~sa517494/)

**Affiliation**: State Key Laboratory of Cognitive Intelligence, USTC; Peking University; Princeton University; Zhejiang University, Tsinghua University

## Contact
We welcome all forms of feedback! Please raise an issue for bugs, questions, or suggestions. This helps our team address common problems efficiently and builds a more productive community. If you encounter any issues, please reach out to jiaxianyan@mail.ustc.edu.cn.

## License
This project is licensed under the terms of the MIT license. See [LICENSE](LICENSE) for additional details.

## Citation
If you find our work helpful, please kindly cite:
```
@article {Yan2025.04.22.648951,
	author = {Yan, Jiaxian and Fang, Xi and Zhu, Jintao and Wu, Chenmin and Yang, Yuhang and Fang, Meijing and Du, Chenxi and Zhang, Kai and Zhang, Zaixi and Liu, Qi},
	title = {Benchmarking LLM Agents on Real-World Biological Database Curation for Data-Driven Scientific Discovery},
	year = {2026},
	journal = {underview}
}
```
