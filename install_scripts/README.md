# Biomni Environment Setup

This directory contains scripts and configuration files to set up a comprehensive bioinformatics environment with various tools and packages.

1. Clone the repository:
   ```bash
   git clone https://github.com/snap-stanford/Biomni.git
   cd Biomni/biomni_env
   ```

2. Setting up the environment:
- (a) run the following script:

```bash
conda env create -f environment.yml
bash new_software_v005.sh
bash setup_biodatalab_tools.sh
```

Note: we have only tested this setup.sh script with Ubuntu 22.04, 64 bit.

3. Lastly, to activate the biomni environment:
```bash
conda activate biomni_e1
```
