# &#19918; <img align="right" width="300" src="https://github.com/sunyatin/qna/blob/main/archives/model.png">

This repository contains all scripts and files to reproduce the results presented in the study "*Questioning Neanderthal admixture*" (doi: https://doi.org/10.1101/2023.04.05.535686). This repository also contains the [`demes`](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html)-formatted encoding of the twenty accepted scenarios from the introduced structured model. Scripts can also be used to perform novel simulations or inference. The data[^1] are contained in another repository named [qna_data](https://github.com/sunyatin/qna_data).

## Layout

```
archives/			# contains the accepted demographic models, the observed statistics and the list of accepted runs
bin/				# external binaries
param_files/			# parameter files for various types of analyses
scripts/			# all custom python3 scripts
scripts_slurm			# bash scripts for computations run on SLURM

./
	general.sh		# general pipeline to simulate genetic data and calculate statistics
	further.sh		# further simulations/analyses for robustness assessment (cf. Supplementary Materials)
	analyses.R		# R script for final statistical analyses (run selection, model comparison) and plotting

	bo_1k5f.est		# structured model with prior parameter distributions
	published.est		# prior distribution on mutation rates for published model simulations
```

Note that all relevant information can be found in the README files of each respective directory. The `archives/` directory contains the [`demes`](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html)-formatted histories, observed statistics.

> The `general.sh` file contains the full pipeline to simulate genetic data (incl. aDNA) and to compute genetic summary statistics. If you want to re-run the pipeline from scratch, you must perform the run selection analysis (within `analyses.R`) after simulating the *n* data from the  parameter prior distributions (see the first entry of the `general.sh` script).


# Setup

> All scripts were run on **Linux**. If your OS is Windows, I would advice using a virtual machine with Linux.

Create a directory `gitdir/qna` in your HOME folder and git clone the repository inside:

```bash
# Create a local directory in your HOME folder
mkdir -p ~/gitdir/qna

# Change current folder
cd ~/gitdir/qna

# Git clone the repository
git clone git@github.com:sunyatin/qna.git .
```

Then, create a `conda` or [`miniconda`](https://docs.conda.io/en/latest/miniconda.html) environment. Here, we create an environment called "*qna*":

```bash
# Create the environment
conda create --name qna

# Activate
conda activate qna

# Add channel
conda config --add channels conda-forge

# Install gsl and openblas:
conda install -c conda-forge gsl
conda install -c conda-forge openblas
```

## Requirements

You will need the following: **python 3.7+** (check that *python3.7-dev* is also installed), **R 3.6+**, **openjdk 11.0+**, **gsl**, **openblas**.

To install required libraries for python3, within the `conda` environment: `python3[.7] -m pip install requirements_python.txt`

To install required libraries for R, within the `conda` environment and a R session:
```R
pkg <- read.table("requirements_R.txt", header=F, stringsAsFactors=F)[,1]
pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
cat("Packages that will be installed: "); print(pkg)
if(length(pkg)) install.packages(pkg)
```
> If installation fails for some packages (esp. *minpack.lm* or *stringi* in R), try using `conda` directly, e.g.: `conda install NAME_OF_THE_PACKAGE`

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

# Usage

If you want to simulate data for other runs that the ones selected for the study, you'll have to change the values of these variables in the `further.sh` file: `ACCEPTED_RUNS_LIST_FILE` & `ACCEPTED_RUNS_PAR_DIR` (in the first lines)


## Reproducing analyses and figures from Tournebize & Chikhi (2023)

To reproduce the analyses and plotting, first download the simulated data from the [qna_data](https://github.com/sunyatin/qna_data) repository.

```bash
# Create a new local directory in your HOME folder
mkdir -p ~/gitdir/qna_data

# Set the current directory
cd ~/gitdir/qna_data

# Clone the qna_data repository inside
git clone git@github.com:sunyatin/qna_data.git .

# Then copy the content to the qna directory
cp -R Final.Blake ~/gitdir/qna
cp -R genetic_maps ~/gitdir/archives

# Then run the sections of your interest in a R session
# Note that you have to be in the conda environment before running R
Rscript analyses.R
```

## Converting a `demes`-formatted history to a `ms` command

# Contact

remi (dot) tournebize (at) univ-tlse3 (dot) fr


[^1]: Genetic data are in [EIGENSTRAT](https://reich.hms.harvard.edu/software/InputFileFormats) format with a modification: `0` encodes the **ancestral homozygous** genotype (instead of the *derived* homozygous genotype in pure EIGENSTRAT). Subsequently, `2` encodes the derived homozygous genotype. All genotypes are mono- or bi-allelic at most.
