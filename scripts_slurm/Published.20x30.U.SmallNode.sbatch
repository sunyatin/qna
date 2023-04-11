#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks 19
#SBATCH --ntasks-per-node=19
#SBATCH --mem=185000

###SBATCH -J 20x30H
###SBATCH --time=01:00:00
###SBATCH --mail-user=toto@mail.com
###SBATCH --mail-type=END
###SBATCH --gres=gpu:4

export OMP_NUM_THREADS=1
MODEL=$1

module load cuda/9.1.85.3
module purge
module load conda/4.9.2 java
conda activate rtournebize

date
echo ""
echo $SLURM_NTASKS
parallel -j $SLURM_NTASKS --delay .2 --joblog logs/${MODEL}.log bash scripts_slurm/Published.20x30.U.sh ::: ${MODEL} ::: {1..50} ::: 0

echo ""
date

#___
