#________________Notes__________________#
#D> Date
#T> Topic
#G> Genome size or properties
#R> Computational resources
#N> Special note
#_______________________________________#
export SCRIPT_SIMULATE=scripts/simulate.py
export SCRIPT_SUMSTATS=scripts/sumstats.py
export BIN_DIR=~/gitdir/qna/bin
export SCRIPT_FURTHER=scripts/further.py


mkdir -p Final.Blake/further
ACCEPTED_RUNS_LIST_FILE=archives/list_accepted_runs.txt
ACCEPTED_RUNS_PAR_DIR=archives/accepted_runs/par_yaml

#===========================================================================#
#T> Simulate 1M parameter combinations with statistical pre-screening
#G> G=10x7Mbp
#G> recomb uniform
#===========================================================================#

CMD="$SCRIPT_SIMULATE \
--model qian_demo \
--write_geno1 \
--bin_dir $BIN_DIR \
--remove_dataset \
-i bo_1k5f.est \
-o Final.Blake/data/1M/\{} \
-g 10 7 \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--seed \{} \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar --stats stats SE sprime crf psmc ld \
--ceu CEU --yri YRI --vindija Vindija --altai Altai --psmc_pops CEU \
--conditions_1 pi.CEU:0.5:1.5 pi.YRI:0.5:1.8 \
--conditions_2 pi.CEU:0.5:1.5 pi.YRI:0.5:1.8 Fst:0.05:0.35 AFS.CEU.1:0.18:0.35 AFS.YRI.1:0.28:0.42 SPRIME.match.mean:0.65:0.90 SPRIME.length.mean:100000:600000 DCFS.1:0.2:0.4"
seq 1 1000000 | parallel --jobs 20 python3.7 $CMD


#===========================================================================#
#D> 04.01.2023
#T> Simulate larger genome sizes for the 20 accepted runs
#G> 20x30Mbp
#R> ocean
#===========================================================================#

#N> seed=42
#G> recomb uniform
CMD="$SCRIPT_SIMULATE \
--model qian_demo \
--write_geno1 \
--bin_dir $BIN_DIR \
-i $ACCEPTED_RUNS_PAR_DIR/\{}.par \
-o Final.Blake/data/20x30Mbp/\{} \
-g 20 30 --mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--seed 42 \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar --stats stats SE sprime psmc ld crf \
--ceu CEU --yri YRI --vindija Vindija --altai Altai --psmc_pops CEU YRI"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD


#N> seed=42
#G> recomb with hotspots
CMD="$SCRIPT_SIMULATE \
--model qian_demo \
--write_geno1 \
--bin_dir $BIN_DIR \
-i $ACCEPTED_RUNS_PAR_DIR/\{}.par \
-o Final.Blake/data/20x30Mbp.HOTSPOTS/\{} \
-g 20 30 --mut_model 'BinaryMutationModel(False)' \
--do_hotspots --hotspots_seed 42 \
--algo hudson \
--seed 42 \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar --stats stats SE sprime psmc ld crf \
--ceu CEU --yri YRI --vindija Vindija --altai Altai --psmc_pops CEU YRI"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD


#===========================================================================#
#D> 04.01.2023
#T> Simulate 50 runs for 11 published models with variable mutation rates
#G> 10x70Mbp
#G> recomb uniform
#N> variable mutation rates
#R> ocean
#===========================================================================#

fun() {
MODEL=$1
ITER=$2
mkdir -p Final.Blake/data/10x7Mbp.50sim/${MODEL}
CORE="$SCRIPT_SIMULATE \
-o Final.Blake/data/10x7Mbp.50sim/${MODEL}/${MODEL}_${ITER} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i published.est \
--mut_model BinaryMutationModel(False) \
--algo hudson \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed $ITER \
--model $MODEL \
-g 10 7 \
--ceu CEU --yri YRI --vindija Vindija \
--psmc_pops CEU YRI"
##
if [ "$MODEL" == "kamm19" ]; then
SAMPLES="Mbuti:50:0:YRI Sardinian:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "ragsdale19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "jacobs19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NeaA:1:50000:Vindija"
elif [ "$MODEL" == "gower21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea:1:50000:Vindija"
elif [ "$MODEL" == "iasi21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50:Vindija"
elif [ "$MODEL" == "durvasula20" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "fu14" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "yang12" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50000:Vindija"
elif [ "$MODEL" == "skov20" ]; then
SAMPLES="Africa:50:0:YRI Iceland:50:0:CEU Vindija:1:50000:Vindija"
elif [ "$MODEL" == "moorjani16" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "schaefer21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Vindija:1:50000:Vindija"
else
echo "$MODEL not implemented"
fi
##
CMD="$CORE --samples $SAMPLES"
python3.7 $CMD
}
export -f fun

parallel --jobs 25 fun ::: kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21 ::: `seq 1 50`


#===========================================================================#
#D> 04.01.2023
#T> Simulate 50 runs for 11 published models with variable mutation rates
#G> 10x70Mbp
#G> recomb with hotspots
#N> variable mutation rates
#R> ocean
#===========================================================================#

fun() {
MODEL=$1
ITER=$2
mkdir -p Final.Blake/data/10x7Mbp.50sim.HOTSPOTS/${MODEL}
CORE="$SCRIPT_SIMULATE \
-o Final.Blake/data/10x7Mbp.50sim.HOTSPOTS/${MODEL}/${MODEL}_${ITER} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i published.est \
--mut_model BinaryMutationModel(False) \
--algo hudson \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed $ITER \
--model $MODEL \
-g 10 7 \
--do_hotspots --hotspots_seed $ITER \
--ceu CEU --yri YRI --vindija Vindija \
--psmc_pops CEU YRI"
##
if [ "$MODEL" == "kamm19" ]; then
SAMPLES="Mbuti:50:0:YRI Sardinian:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "ragsdale19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "jacobs19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NeaA:1:50000:Vindija"
elif [ "$MODEL" == "gower21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea:1:50000:Vindija"
elif [ "$MODEL" == "iasi21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50:Vindija"
elif [ "$MODEL" == "durvasula20" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "fu14" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "yang12" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50000:Vindija"
elif [ "$MODEL" == "skov20" ]; then
SAMPLES="Africa:50:0:YRI Iceland:50:0:CEU Vindija:1:50000:Vindija"
elif [ "$MODEL" == "moorjani16" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "schaefer21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Vindija:1:50000:Vindija"
else
echo "$MODEL not implemented"
fi
##
CMD="$CORE --samples $SAMPLES"
python3.7 $CMD
}
export -f fun

parallel --jobs 25 fun ::: kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21 ::: `seq 1 50`


#===========================================================================#
#D> 10.01.2023
#T> Simulate 50 runs for 11 published models with variable mutation rates
#G> 20x30Mbp
#G> recomb with hotspots
#N> variable mutation rates
#R> ocean
#===========================================================================#

fun() {
MODEL=$1
ITER=$2
mkdir -p Final.Blake/data/20x30Mbp.50sim.HOTSPOTS/${MODEL}
CORE="$SCRIPT_SIMULATE \
-o Final.Blake/data/20x30Mbp.50sim.HOTSPOTS/${MODEL}/${MODEL}_${ITER} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i published.est \
--mut_model BinaryMutationModel(False) \
--algo hudson \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed $ITER \
--model $MODEL \
-g 20 30 \
--do_hotspots --hotspots_seed $ITER \
--ceu CEU --yri YRI --vindija Vindija \
--psmc_pops CEU YRI"
##
if [ "$MODEL" == "kamm19" ]; then
SAMPLES="Mbuti:50:0:YRI Sardinian:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "ragsdale19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "jacobs19" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NeaA:1:50000:Vindija"
elif [ "$MODEL" == "gower21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea:1:50000:Vindija"
elif [ "$MODEL" == "iasi21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50:Vindija"
elif [ "$MODEL" == "durvasula20" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Neanderthal:1:50000:Vindija"
elif [ "$MODEL" == "fu14" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "yang12" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU NEA:1:50000:Vindija"
elif [ "$MODEL" == "skov20" ]; then
SAMPLES="Africa:50:0:YRI Iceland:50:0:CEU Vindija:1:50000:Vindija"
elif [ "$MODEL" == "moorjani16" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Nea_1:1:50000:Vindija"
elif [ "$MODEL" == "schaefer21" ]; then
SAMPLES="YRI:50:0:YRI CEU:50:0:CEU Vindija:1:50000:Vindija"
else
echo "$MODEL not implemented"
fi
##
CMD="$CORE --samples $SAMPLES"
python3.7 $CMD
}
export -f fun

parallel --jobs 25 fun ::: kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21 ::: `seq 1 50`

tar --exclude "*.gz" --exclude "*.ind" --exclude "*.psmcfa" --exclude "*.geno" --exclude "*.snp" --exclude "*.f0?" --exclude "*.freq" --exclude "*.map" --exclude "*.vcf" -zvcf 20x30Mbp.50sim.HOTSPOTS.tar.gz 20x30Mbp.50sim.HOTSPOTS

#===========================================================================#
#D> 10.01.2023
#T> Simulate 50 runs for 11 published models with variable mutation rates
#G> 20x30Mbp
#G> recomb uniform
#N> variable mutation rates
#R> olympe/calmip
#===========================================================================#

# dispatch on GPU nodes
for MODEL in kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 schaefer21; do
sbatch --job-name "$MODEL" scripts_slurm/Published.20x30.U.sbatch $MODEL
done
# dispatch on smaller nodes
for MODEL in skov20 moorjani16; do
sbatch --job-name "$MODEL" scripts_slurm/Published.20x30.U.SmallNode.sbatch $MODEL
done
# rerun for sims that failed due to memory saturation
for MODEL in kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21; do
sbatch --job-name "$MODEL" scripts_slurm/Published.20x30.U.sbatch $MODEL
done
# check that all models have 50 sims
for MODEL in kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21; do
echo "$MODEL : $(ls $MODEL/*crf.stats | wc -l)"
done


#===========================================================================#
#D> 19.01.2023
#T> Simulate single run for 11 published models with mu=1.2e-8
#G> 20x30Mbp
#G> recomb uniform
#N> mutation rate fixed at 1.2e-8
#R> olympe/calmip
#===========================================================================#

# dispatch on GPU nodes
for MODEL in kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21; do
sbatch --job-name "$MODEL" scripts_slurm/Published.20x30.U.1_2e-8.sbatch $MODEL
done
# check that all models have 50 sims
for MODEL in kamm19 ragsdale19 jacobs19 gower21 iasi21 durvasula20 fu14 yang12 skov20 moorjani16 schaefer21; do
echo "$MODEL : $(ls $MODEL/*crf.stats | wc -l)"
done


#===========================================================================#
#!
#D> 04.01.2023
#T> Simulate the 20 accepted runs
#G> 10x7Mbp
#G> recomb uniform
#R> ocean
#===========================================================================#

CMD="$SCRIPT_SIMULATE \
-o Final.Blake/data/10x7Mbp.50sim/1M.accepted/\{} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i $ACCEPTED_RUNS_PAR_DIR/\{}.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed \{} \
--model qian_demo \
-g 10 7 \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 25 python3.7 $CMD


#===========================================================================#
#!
#D> 04.01.2023
#T> Simulate the 20 accepted runs
#G> 10x7Mbp
#G> recomb with hotspots
#R> ocean
#===========================================================================#

CMD="$SCRIPT_SIMULATE \
-o Final.Blake/data/10x7Mbp.50sim.HOTSPOTS/1M.accepted/\{} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i $ACCEPTED_RUNS_PAR_DIR/\{}.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed \{} \
--model qian_demo \
-g 10 7 \
--do_hotspots --hotspots_seed \{} \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 25 python3.7 $CMD


#===========================================================================#
#D> 05.01.2023
#T> aDNA study : simulate the 20 accepted runs with ancient samples
#G> 20x30Mbp
#G> recomb uniform
#R> ocean
#===========================================================================#

#N> Simulate genetic data with ancient Eurasians sampled in CEU:4
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU4 \
param_files/ocean.U.L.CEU4.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD | tee Final.Blake/further/20x30Mbp.U.CEU4.log

#N> Simulate genetic data with ancient Eurasians sampled in CEU:0,3,6,9
#G> `All` SNPs
#G> diploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U \
param_files/ocean.U.L.di.All.0369.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD | tee Final.Blake/further/20x30Mbp.U.di.All.0369.log

#===============#
#     CEU:4     #
#===============#

# Analyze the aDNA data sampled in CEU:4
#G> `All` SNPs
#G> diploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU4 \
param_files/ocean.U.L.di.All.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 7 python3.7 $CMD | tee Final.Blake/further/20x30Mbp.U.di.All.log

# Analyze the aDNA data sampled in CEU:4
#G> `1M` SNPs
#G> diploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU4 \
param_files/ocean.U.L.di.1M.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 7 python3.7 $CMD | tee Final.Blake/further/20x30Mbp.U.di.1M.log

# Analyze the aDNA data sampled in CEU:4
#G> `Archaic` SNPs
#G> pseudodiploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU4 \
param_files/ocean.U.L.pseudodi.Archaic.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 7 python3.7 $CMD | tee Final.Blake/further/20x30Mbp.U.pseudodi.Archaic.log

# Analyze the aDNA data sampled in CEU:4
#G> `1M` SNPs
#G> pseudodiploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU4 \
param_files/ocean.U.L.pseudodi.1M.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 7 python3.7 $CMD | tee Final.Blake/further/20x30Mbp.U.pseudodi.1M.log


#===============#
#     CEU:3     #
#===============#

# Analyze the aDNA data sampled in CEU:3
#G> `All` SNPs
#G> diploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU3 \
param_files/ocean.U.L.di.All.CEU3.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD

# Analyze the aDNA data sampled in CEU:3
#G> `1M` SNPs
#G> diploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU3 \
param_files/ocean.U.L.di.1M.CEU3.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD

# Analyze the aDNA data sampled in CEU:3
#G> `Archaic` SNPs
#G> pseudodiploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU3 \
param_files/ocean.U.L.pseudodi.Archaic.CEU3.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD

# Analyze the aDNA data sampled in CEU:3
#G> `1M` SNPs
#G> pseudodiploid
CMD="$SCRIPT_FURTHER \
$ACCEPTED_RUNS_PAR_DIR/\{}.par \
Final.Blake/further/20x30Mbp.U.CEU3 \
param_files/ocean.U.L.pseudodi.1M.CEU3.yaml"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 20 python3.7 $CMD


#=================#
#    Schaefer21   #
#=================#

# Simulate and analyze the aDNA data simulated under the published Schaefer et al (2021) model
#N> mutation rate fixed at 1.2e-8
#R> local
python3.7 $SCRIPT_SIMULATE \
-o Final.Blake/further/20x30Mbp.U/stats_aDNA_schaefer21 \
--write_geno1 \
--bin_dir $BIN_DIR \
-i param_files/mutrate_1.2e-8.spar \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model schaefer21 \
-g 20 30 \
--draw param_files/else.dpar \
--samples YRI:50:0:YRI \
CEU:1:0:aDNA.CEU.0.0 \
CEU:1:5000:aDNA.CEU.0.5000 \
CEU:1:10000:aDNA.CEU.0.10000 \
CEU:1:15000:aDNA.CEU.0.15000 \
CEU:1:20000:aDNA.CEU.0.20000 \
CEU:1:25000:aDNA.CEU.0.25000 \
CEU:1:30000:aDNA.CEU.0.30000 \
CEU:1:35000:aDNA.CEU.0.35000 \
CEU:1:40000:aDNA.CEU.0.40000 \
CEU:1:45000:aDNA.CEU.0.45000 \
Vindija:1:50000:Vindija \
Altai:1:130000:Altai \
Nea_Intro:1:50000:NeaIntro
#
mkdir -p Final.Blake/further/20x30Mbp.U/stats_aDNA_schaefer21
CMD="$SCRIPT_SUMSTATS \
--bin_dir $BIN_DIR \
-p param_files/further_stats.spar \
--stats --SE --ld1  \
-i Final.Blake/further/20x30Mbp.U/stats_aDNA_schaefer21 \
-o Final.Blake/further/20x30Mbp.U/stats_aDNA_schaefer21/\{} \
--ceu \{} --yri YRI --vindija Vindija --altai Altai"
parallel --jobs 2 python3.7 $CMD ::: aDNA.CEU.0.0 aDNA.CEU.0.5000 aDNA.CEU.0.10000 aDNA.CEU.0.15000 aDNA.CEU.0.20000 aDNA.CEU.0.25000 aDNA.CEU.0.30000 aDNA.CEU.0.35000 aDNA.CEU.0.40000 aDNA.CEU.0.45000


#===========================================================================#
#D> 24.03.2023
#T> Plot demographic model representation and writes the YAML demo file
#R> local
#===========================================================================#

python3.7 $SCRIPT_SIMULATE \
-o Final.Blake/plots/1014230x \
--bin_dir $BIN_DIR \
-i $ACCEPTED_RUNS_PAR_DIR/1014230.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
-g 1 1 \
-x \
--draw param_files/qian_svg.dpar \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai

CMD="$SCRIPT_SIMULATE \
-o Final.Blake/accepted_runs/\{} \
--bin_dir $BIN_DIR \
-i $ACCEPTED_RUNS_PAR_DIR/\{}.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
-g 1 1 \
-x \
--draw param_files/qian_svg.dpar \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai"
parallel -a $ACCEPTED_RUNS_LIST_FILE --jobs 5 python3.7 $CMD
rm Final.Blake/accepted_runs/*.lprob

#===========================================================================#
# Misc.
#===========================================================================#

# Check if any duplicate position in a *.snp.gz file:
for FILE in *snp.gz; do zcat $FILE | cut -f2,4 -d" " | sed "s/\s/./" | sort | uniq -c | tr -s " " | cut -d" " -f2 | sort | uniq; done


#___
