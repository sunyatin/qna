# &#19918; <img align="right" width="300" src="https://github.com/sunyatin/qna/blob/main/archives/model.png">

This repository contains all scripts and files to reproduce the results presented in the study "*Questioning Neanderthal admixture*" (doi: https://doi.org/10.1101/2023.04.05.535686). This repository also contains the [`demes`](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html)-formatted encoding of the twenty accepted scenarios from the introduced structured model. Scripts can also be used to perform novel simulations or inference. The data[^1] are contained in another repository named [qna_data](https://github.com/sunyatin/qna_data).

The genetic data[^1] (modified EIGENSTRAT format) for the twenty accepted runs of our structured model can be found in the qna_data repository, [here](https://github.com/sunyatin/qna_data/tree/main/Final.Blake/data/20x30Mbp).

## Layout

```
archives/			# contains the accepted demographic models, the observed statistics and the list of accepted runs
	accepted_runs/		# demographic models of the twenty accepted runs
		model_plots/	# png and svg plots of the accepted models
		par_yaml/	# YAML and parameter files of the accepted models
	obs/			# observed/empirical values of the genetic summary statistics
		estimation/	# contains the script and data for estimation of CRF and S' values

bin/				# external binaries

param_files/			# parameter files for various types of analyses

scripts/			# all custom python3 scripts

scripts_slurm/			# bash scripts for computations run on SLURM

./
	general.sh		# general pipeline to simulate genetic data and calculate statistics
	further.sh		# further simulations/analyses for robustness assessment (cf. Supplementary Materials)
	analyses.R		# R script for final statistical analyses (run selection, model comparison) and plotting

	bo_1k5f.est		# structured model with prior parameter distributions
	published.est		# prior distribution on mutation rates for published model simulations

	requirements_python.txt	# library requirements for python3
	requirements_R.txt	# library requirements for R
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
```

## Requirements

You will need the following: **python 3.7+** (check that `python3.7-dev` is also installed), **R 3.6.3**, **openjdk 11.0+**, **gsl**, **openblas**:

```bash
conda install -c conda-forge python=3.7 r-base=3.6.3 openjdk=11.0 gsl openblas
```

To install required libraries for python3, within the `conda` environment:

```bash
python3[.7] -m pip install -r requirements_python.txt
```

To install required libraries for R, within the `conda` environment, first run:
```bash
conda install -c conda-forge r-lattice r-mass
```
Then, **in a R session**:
```R
pkg <- read.table("requirements_R.txt", header=F, stringsAsFactors=F)[,1]
pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
cat("Packages that will be installed: "); print(pkg)
if(length(pkg)) install.packages(pkg)
```
> If installation fails for some packages (esp. *minpack.lm* or *stringi* in R), try using `conda` directly: `conda install -c conda-forge r-{NAME_OF_THE_PACKAGE}`. For instance: `conda install -c conda-forge r-stringi`

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

## General important notes

For all calls to the `simulate.py` script, you need to make sure that the alias `python3` in your Linux console refers to the version 3.7+:

```bash
# Check version of python3
python3 --version

# If the version is lower than 3.7, change the alias in the current session:

# First locate the path to python3.7
which python3.7
#=> Would usually return:
# /usr/bin/python3.7

# Change the alias with what the previous command returned
alias python3="/usr/bin/python3.7"

# Note that if you want to permanently change the alias, open the ~/.bashrc file
nano ~/.bashrc

# Add the alias
alias python3="/usr/bin/python3.7"

# Save the file and reload it
source ~/.bashrc

# Check that the version is the proper one
python3 --version
#=> It should now be 3.7+
```

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

## Running the model from scratch

If you want to re-run the pipeline from scratch, or re-simulate or re-analyze data for other runs than the ones selected in Tournebize & Chikhi (2023), the values of two variables should be changed in `general.sh` and `further.sh`:
- `ACCEPTED_RUNS_LIST_FILE` path to a file containing the ID of the selected runs, on separate lines
- `ACCEPTED_RUNS_PAR_DIR` path to a directory containing the parameter files (`*.par`) of the selected runs

```bash
# Simulate a prior distribution of genetic data
Run the first section in general.sh

# Select the runs
Run the first section in analyses.R

# Change the values of the variables in general.sh
ACCEPTED_RUNS_LIST_FILE=your_file_containing_the_runs_selected_in_the_previous_R_analysis
ACCEPTED_RUNS_PAR_DIR=the_directory_of_your_prior_simulations

# Run additional simulations (e.g. with larger genome size)
Run relevant sections in general.sh
```

## Simulate and calculate summary statistics for a particular model

The following command will simulate genetic data under the run "1014230". It will simulate 10 chromosomes of 30 Mbp each, sample 50 individuals in YRI and CEU, 1 Neanderthal at time 50 kya, and calculate all genetic summary statistics. The output files will be written into the repository `test/` with prefix `output_1014230.*`.

```bash
python3.7 scripts/simulate.py \
--model qian_demo \
--write_geno1 \
--bin_dir bin \
-i archives/accepted_runs/par_yaml/1014230.par \
-o test/output_1014230 \
-g 10 30 \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s scripts/sumstats.py param_files/stats_MAC_diplo.spar \
--stats stats SE sprime psmc ld crf \
--ceu CEU --yri YRI --vindija Vindija --altai Altai --psmc_pops CEU YRI \
--draw param_files/qian_png.dpar
```

## Converting between `demes` or `msprime` models and `ms` commands

- ***ms => demes***
```python
import demes
# No is the reference effective size
Model = demes.from_ms(ms_command, N0=No)
```

- ***ms => msprime***
```python
import demes, msprime
# No is the reference effective size
Model = msprime.Demography.from_demes(demes.from_ms(ms_command, N0=No))
```

- ***msprime => ms***
```python
import demes, msprime
# msprimeDemography is the demographic model in the `msprime` format
# No is the reference effective size
ms_command = demes.to_ms(msprimeDemography.to_demes(), N0=No)
print(ms_command)
```

- ***demes => ms***
```python
import demes
# No is the reference effective size
Model = demes.load(path_to_the_yaml_file)
ms_command = demes.to_ms(Model, N0=No)
print(ms_command)
```

- ***demes => msprime***
```python
import demes, msprime
Model = msprime.Demography.from_demes(demes.load(path_to_the_yaml_file))
```

# Notes on the observed data

Observed values for the genetic summary statistics are contained in the [`archives/obs`](https://github.com/sunyatin/qna/tree/main/archives/obs) directory. Note that for the S' and the CRF analyses, we used the analyzed data from the original studies, that we further summarized using a custom script which is located in [`archives/obs/estimation/estimation.py`](https://github.com/sunyatin/qna/blob/main/archives/obs/estimation/estimation.py).

# Contact

remi (dot) tournebize (at) univ-tlse3 (dot) fr


[^1]: Genetic data are in [EIGENSTRAT](https://reich.hms.harvard.edu/software/InputFileFormats) format with a modification: `0` encodes the **ancestral homozygous** genotype (instead of the *derived* homozygous genotype in pure EIGENSTRAT). Subsequently, `2` encodes the derived homozygous genotype. All genotypes are mono- or bi-allelic at most.
