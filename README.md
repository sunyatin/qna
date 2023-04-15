<center>
䷎
</center>

<p style='text-align: justify;'>

This repository contains all data and scripts to reproduce the results presented in the study "*Questioning Neanderthal admixture*" (doi: https://doi.org/10.1101/2023.04.05.535686). It also contains the `demes`-formatted encoding of the twenty accepted scenarios from the introduced structured model. Scripts can further be re-used to perform novel simulations or inference.

</p>

At the root of the repository are five directories:

| Folder         | Description                                 |
|----------------|----------------------------------------------|
| **archives**       | Stores simulated genetic data, `demes`-formatted histories, observed statistics, empirical genetic maps.   |
| **bin**            | Stores external softwares.          |
| **param_files**    | Stores parameter files for all types of analyses.     |
| **scripts**        | Stores the *python3* scripts.    |
| **scripts_slurm**  | Stores *bash* scripts for analyses run on SLURM. |


All relevant information can be found in the README files of each respective directory.

At the root of the repository are also five files:

| Scripts         | Description                                 |
|----------------|----------------------------------------------|
| `general.sh` |  ***Bash* pipeline to simulate genetic data and compute genetic summary statistics.** Note that if you want to re-run all the pipeline from scratch, you'll need to perform the *run selection* analysis (within `analyses.R`) after simulating the *n* data from the  parameter prior distributions (corresponds to the first entry of the `general.sh`  script). |
| `further.sh` | ***Bash* pipeline to perform additional simulations and analyses:** for robustness assessment. |
| `analyses.R` | ***R* script for all statistical analyses:** run selection, model comparison, and all figure plotting. |

| Files         | Description                                 |
|----------------|----------------------------------------------|
| `bo_1k5f.est` | Specifies the structured model and the prior distribution of parameters. |
| `published.est` | Specifies the prior distribution of the mutation rate (per generation) to simulate *published* models. |



# General notes

<p style='text-align: justify;'>

Genetic data are in EIGENSTRAT format with a modification: `0` encodes the **ancestral homozygous** genotype (instead of the *derived* homozygous genotype in pure EIGENSTRAT). Subsequently, `2` encodes the derived homozygous genotype. All genotypes are mono- or bi-allelic at most.

</p>

# Recipes
