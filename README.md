# &#19918; <img align="right" width="300" src="https://github.com/sunyatin/qna/blob/main/archives/model.png">

This repository contains all data[^1] and scripts to reproduce the results presented in the study "*Questioning Neanderthal admixture*" (doi: https://doi.org/10.1101/2023.04.05.535686). It also contains the `demes`-formatted encoding of the twenty accepted scenarios from the introduced structured model. Scripts can further be re-used to perform novel simulations or inference.

At the root of the repository are five directories. All relevant information can be found in the README files of each respective directory.

| Folder         | Description                                 |
|----------------|----------------------------------------------|
| **archives**       | Stores simulated genetic data[^1], `demes`-formatted histories, observed statistics, empirical genetic maps.   |
| **bin**            | Stores external softwares.          |
| **param_files**    | Stores parameter files for all types of analyses.     |
| **scripts**        | Stores the *python3* scripts.    |
| **scripts_slurm**  | Stores *bash* scripts for analyses run on SLURM. |


At the root of the repository are also five files:

| Scripts         | Description                                 |
|----------------|----------------------------------------------|
| **general.sh** |  ***Bash* pipeline to simulate genetic data[^1] (incl. aDNA) and compute genetic summary statistics.** Note that if you want to re-run all the pipeline from scratch, you'll need to perform the run selection analysis (within `analyses.R`) after simulating the *n* data from the  parameter prior distributions (corresponds to the first entry of the `general.sh`  script). |
| **further.sh** | ***Bash* pipeline to perform additional simulations and analyses:** for robustness assessment (cf. SupMat). |
| **analyses.R** | ***R* script for all statistical analyses:** run selection, model comparison, and all figure plotting. |

| Files         | Description                                 |
|----------------|----------------------------------------------|
| **bo_1k5f.est** | Specifies the structured model and the prior distribution of parameters. |
| **published.est** | Specifies the prior distribution of the mutation rate (per generation) to simulate *published* models. |

# :large_blue_diamond: Setup

> Note that scripts run on **UNIX**/**Linux**. If your OS is Windows, please use a virtual machine.

Download the repository archive. For easier installation (to avoid changing absolute paths), create a directory in your HOME folder (`~`) and decompress the archive content inside at:

> `~/gitdir/qna`

Mandatory programs:
- **python 3.7+** (check that *python3.7-dev* is also installed)
- **R 3.6+**
- **openjdk 11.0+**
- **gsl**
- **openblas**

## Conda environment

Create a *conda* environment. Here, we create an environment called "*qna*":

```bash
# Create the qna environment
conda create --name qna

# Activate
conda activate qna

# Add channel
conda config --add channels conda-forge

# Install gsl and openblas:
conda install -c conda-forge gsl
conda install -c conda-forge openblas
```

## Dependencies

### **python3.7+**
- numpy *[1.21.6]*
- scipy *[1.7.3]*
- pandas *[1.3.5]*
- seaborn *[0.11.2]*
- msprime *[1.1.0]*
- demes *[0.2.2]*
- demesdraw *[0.3.0]*
- stdpopsim *[0.1.3b1]*
- scikit-allel *[1.3.5]*

Commands:
```python
python3.7 -m pip install numpy==1.21.6
python3.7 -m pip install scipy==1.7.3
python3.7 -m pip install pandas==1.3.5
python3.7 -m pip install seaborn==0.11.2
python3.7 -m pip install msprime==1.1.0
python3.7 -m pip install demes==0.2.2
python3.7 -m pip install demesdraw==0.3.0
python3.7 -m pip install stdpopsim==0.1.3
python3.7 -m pip install scikit-allel==1.3.5

```

Other dependencies:
- argparse
- copy
- decimal
- gzip
- logging
- matplotlib
- os
- random
- re
- shutil
- subprocess
- sys
- time
- yaml

### R
Within the *conda* environment, install within *R* using `install.packages(name_of_library)`:
- cowplot
- corrplot
- dplyr
- ggbeeswarm
- ggdist
- ggbump
- ggmap
- ggplot2
- ggsci
- ggridges
- minpack.lm
- paletteer
- plotrix
- reshape2
- scales
- scico
- viridis

To install all at once, in a *R* session:
```R
pkg <- c("cowplot", "corrplot", "dplyr", "ggbeeswarm", "ggdist", "ggbump", "ggmap", "ggplot2", "ggsci", "ggridges", "minpack.lm", "paletteer", "plotrix", "reshape2", "scales", "scico", "viridis")
pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
cat("Packages that will be installed: "); print(pkg)
if(length(pkg)) install.packages(pkg)
```

If installation fails for some packages, try using *conda*, e.g.:
```bash
conda install r-minpack.lm
conda install r-stringi
```

## External programs

For legal considerations, we do not include external programs in the repository. Users are asked to download them separately and store the executable in the `~/gitdir/qna/bin` directory:

- `sprime.jar` (*S'* analysis) : download from [here](https://faculty.washington.edu/browning/sprime.jar).
- `psmc` software (*PSMC* analysis) : download from [here](https://github.com/lh3/psmc).
- `caller-static-2015` (*CRF* analysis) : to be asked from the [authors](https://doi.org/10.1038/nature12961).
- executable of the single-sample `neanderthal_dating` (to be stored within `bin/Neanderthal_dating/bin/`) : download from [here](https://github.com/priyamoorjani/Neanderthal_dating).

## Troubleshooting

If you get error messages related to missing libraries, consider changing the LD_LIBRARY_PATH variable. For instance, if you get a message stating that "libgsl.so.25" is not found, try:

```bash
# Locate the library path
sudo find / -name "libgsl.so.25"
# Usually, it will be in your HOME folder, under conda/envs/qna/lib or miniconda3/envs/qna/lib

# Replace the LD_LIBRARY_PATH variable with the path returned by the previous command
# This will have to be done each time you open a new terminal.
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/USERNAME/miniconda3/envs/qna/lib
export LD_LIBRARY_PATH
```

# :large_blue_diamond: Usage

If you want to simulate data for other runs that the ones selected for the study, you'll have to change the values of these variables in the `further.sh` file: `ACCEPTED_RUNS_LIST_FILE` & `ACCEPTED_RUNS_PAR_DIR` (in the first lines)

# :large_blue_diamond: Recipes

## Re-running statistical analyses using the original data

Assuming you are in `~/gitdir/qna`, if you want rerun the statistical analyses using the genetic data that were simulated for the original paper, just copy the `archives/Final.Blake` directory to `~/gitdir/qna`.

```bash
cd ~/gitdir/qna
cp -r archives/Final.Blake .
```
## Converting a `demes`-formatted history to a `ms` command

# :large_blue_diamond: Contact

remi (dot) tournebize (at) univ-tlse3 (dot) fr


[^1]: Genetic data are in EIGENSTRAT format with a modification: `0` encodes the **ancestral homozygous** genotype (instead of the *derived* homozygous genotype in pure EIGENSTRAT). Subsequently, `2` encodes the derived homozygous genotype. All genotypes are mono- or bi-allelic at most.
