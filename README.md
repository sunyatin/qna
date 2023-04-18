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
| **general.sh** |  ***Bash* pipeline to simulate genetic data[^1] (incl. aDNA) and compute genetic summary statistics.** Note that if you want to re-run all the pipeline from scratch, you'll need to perform the *run selection* analysis (within `analyses.R`) after simulating the *n* data from the  parameter prior distributions (corresponds to the first entry of the `general.sh`  script). |
| **further.sh** | ***Bash* pipeline to perform additional simulations and analyses:** for robustness assessment (cf. SupMat). |
| **analyses.R** | ***R* script for all statistical analyses:** run selection, model comparison, and all figure plotting. |

| Files         | Description                                 |
|----------------|----------------------------------------------|
| **bo_1k5f.est** | Specifies the structured model and the prior distribution of parameters. |
| **published.est** | Specifies the prior distribution of the mutation rate (per generation) to simulate *published* models. |

# :large_blue_diamond: Setup

> Note that the pipelines can only run on **UNIX**/**Linux**. If you use Windows, you could use a virtual machine.

Download the repository archive. For easier installation (to avoid changing absolute paths), create a directory in your HOME folder: `~/gitdir/qna` and uncompress the archive content inside.

If you want to simulate data for other runs that the ones selected for the study, you'll have to change the values of these variables in the `further.sh` file: `ACCEPTED_RUNS_LIST_FILE` & `ACCEPTED_RUNS_PAR_DIR`

You will need to have:
- python 3.7+ (check that python3.7-dev is also installed, `sudo apt-get install python3.7-dev`)
- R 3.6+ (`sudo apt-get install r-base=3.6.4`)
- openjdk 11.0+
- gcc

## Dependencies



wget https://repo.anaconda.com/archive/Anaconda3-2023.03-Linux-x86_64.sh
cd XX
bash Anaconda3-2023.03-Linux-x86_64.sh

conda create --name qna
conda activate qna
conda install -c conda-forge r-base=4.1.2


conda config --add channels conda-forge
conda install msprime>=1.1.0



For **python3.7+** (install using `python3.7 -m pip3 install name_of_module` or alternatively with `conda`):
- cython
- numpy
- pybind11
- pythran
- scipy
- pandas

Then:
- scikit-allel *[>= 1.3.5]* `python3.7 -m pip3 install scikit-allel>=1.3.5`
- argparse
- copy
- decimal
- demes *[>= 0.2.2]* `python3.7 -m pip3 install demes==0.2.2`
- gzip
- logging
- msprime *[>=1.1.0]* `python3.7 -m pip3 install msprime>=1.1.0`
- os
- re
- shutil
- stdpopsim *[>= 0.1.3b1]* `python3.7 -m pip3 install stdpopsim>=0.1.3`
- subprocess
- sys
- time
- yaml

For **R** (install using `install.packages(name_of_library)`):
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

To install all at once, in a R session:

```R
pkg <- c("cowplot", "corrplot", "dplyr", "ggbeeswarm", "ggdist", "ggbump", "ggmap", "ggplot2", "ggsci", "ggridges", "minpack.lm", "paletteer", "plotrix", "reshape2", "scales", "scico", "viridis")
pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(pkg)) install.packages(pkg)
```


## External programs

For legal considerations, we do not include external programs in the repository. Users are asked to download them separately and store the executable in the **bin** directory:

- `sprime.jar` (*S'* analysis) : download from [here](https://faculty.washington.edu/browning/sprime.jar).
- `psmc` software (*PSMC* analysis) : download from [here](https://github.com/lh3/psmc).
- `caller-static-2015` (*CRF* analysis) : to be asked from the [authors](https://doi.org/10.1038/nature12961).
- executable of the single-sample `neanderthal_dating` (to be stored within `bin/Neanderthal_dating/bin/`) : download from [here](https://github.com/priyamoorjani/Neanderthal_dating).

# :large_blue_diamond: Recipes

## Re-running statistical analyses using original data

Assuming you are in `~/gitdir/qna`, if you want rerun the statistical analyses using the genetic data that were simulated for the original paper, just copy the "*Final.Blake*" directory to the current directory.

```bash
cp -r archives/Final.Blake .
```
## Converting a `demes`-formatted history to a `ms` command

# :large_blue_diamond: Contact

remi (dot) tournebize (at) univ-tlse3 (dot) fr


[^1]: Genetic data are in EIGENSTRAT format with a modification: `0` encodes the **ancestral homozygous** genotype (instead of the *derived* homozygous genotype in pure EIGENSTRAT). Subsequently, `2` encodes the derived homozygous genotype. All genotypes are mono- or bi-allelic at most.
