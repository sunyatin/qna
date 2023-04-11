#________________Notes__________________#
#D> Date
#T> Topic
#G> Genome size or properties
#R> Computational resources
#N> Special note
#_______________________________________#

RUN=1014230
DATA_DIR=Final.Blake/data/20x30Mbp
DATA_OUT=Final.Blake/further/sensitivity

mkdir -p $DATA_OUT



#===========================================================================#
#T> Analyze data with ::: MAC filtering
#R> local
#===========================================================================#

python3.7 $SCRIPT_SUMSTATS \
--assume_no_na \
--bin_dir $BIN_DIR \
-p param_files/stats_MAC_diplo.spar \
--stats --SE --sprime --crf --psmc --ld  \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
-i $DATA_DIR/$RUN \
-o $DATA_OUT/rt.$RUN.MAC

#===========================================================================#
#T> Analyze data with ::: no MAC filtering
#R> local
#===========================================================================#

python3.7 $SCRIPT_SUMSTATS \
--assume_no_na \
--bin_dir $BIN_DIR \
-p param_files/stats_noMAC_diplo.spar \
--stats --SE --sprime --crf --psmc --ld  \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
-i $DATA_DIR/$RUN \
-o $DATA_OUT/rt.$RUN.noMAC

#===========================================================================#
#T> Analyze data with ::: pseudodiploid genotypes
#R> local
#===========================================================================#

python3.7 $SCRIPT_SUMSTATS \
--assume_no_na \
--bin_dir $BIN_DIR \
-p param_files/stats_MAC_pseudodiplo.spar \
--stats --SE --sprime --crf --ld  \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
-i $DATA_DIR/$RUN \
-o $DATA_OUT/rt.$RUN.pseudodiplo
# For LD stats, can only be done with pseudodiplo if *no* MAC filtering:
python3.7 $SCRIPT_SUMSTATS \
--assume_no_na \
--bin_dir $BIN_DIR \
-p param_files/stats_noMAC_pseudodiplo.spar \
--stats --SE --sprime --crf --ld  \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
-i $DATA_DIR/$RUN \
-o $DATA_OUT/rt.$RUN.noMAC.pseudodiplo

#===========================================================================#
#T> Analyze data where ::: Neand sample is Altai
#R> local
#===========================================================================#

python3.7 $SCRIPT_SUMSTATS \
--assume_no_na \
--bin_dir $BIN_DIR \
-p param_files/stats_MAC_diplo.spar \
--stats --SE --sprime --crf --psmc --ld  \
--ceu CEU --yri YRI --vindija Altai --altai Altai \
--psmc_pops CEU YRI \
-i $DATA_DIR/$RUN \
-o $DATA_OUT/rt.$RUN.Altai

#===========================================================================#
#T> Analyze data with ::: ancestral polarization error
#R> local
#===========================================================================#

# Add the random errors
python3.7 scripts/ancestral_errors.py \
-i $DATA_DIR/$RUN \
-o $DATA_OUT/rt.$RUN.AncError3p \
-p 0.03 \
--seed 42

# Analyze
python3.7 $SCRIPT_SUMSTATS \
--assume_no_na \
--bin_dir $BIN_DIR \
-p param_files/stats_MAC_diplo.spar \
--stats --SE --sprime --crf --psmc --ld  \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI \
-i $DATA_OUT/rt.$RUN.AncError3p \
-o $DATA_OUT/rt.$RUN.AncError3p

#===========================================================================#
#T> Compare ancestry-LD population-level vs. individual-level algorithms
#R> local
#N> No MAC-filtering
#===========================================================================#

NCHROM=10

mkdir -p Final.Blake/further/LD_implementation

#N> Population-LD algorithm
POPULATION="$SCRIPT_SIMULATE \
-o Final.Blake/further/LD_implementation/\{} \
--write_geno1 \
--bin_dir $BIN_DIR \
-i Final.Blake/further/sankararaman12_NGFI_Admix.est \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
--seed \{} \
-g $NCHROM 10 \
--samples_within AfW:0:50:0:YRI EuA:0:50:0:CEU Nea:0:1:50000:Vindija && \
python3.7 $SCRIPT_SUMSTATS \
--bin_dir $BIN_DIR \
-p param_files/stats_noMAC_diplo.spar \
--ld \
-i Final.Blake/further/LD_implementation/\{} \
--ceu CEU --yri YRI --vindija Vindija \
--psmc_pops CEU \
-o Final.Blake/further/LD_implementation/\{}.AsSanka12 \
-v --ld_as_sanka12 && \
python3.7 $SCRIPT_SUMSTATS \
--bin_dir $BIN_DIR \
-p param_files/stats_noMAC_diplo.spar \
--ld \
-i Final.Blake/further/LD_implementation/\{} \
--ceu CEU --yri YRI --vindija Vindija \
--psmc_pops CEU \
-o Final.Blake/further/LD_implementation/\{}.RT \
-v && \
bash scripts/computed.sh $NCHROM Final.Blake/further/LD_implementation/\{}"

parallel --jobs 5 python3.7 $POPULATION ::: {1..30}

#N> Individual-LD algorithm
SINGLEIND="$SCRIPT_SIMULATE \
-o Final.Blake/further/LD_implementation/data/\{}.ld1 \
--write_geno1 \
--bin_dir $BIN_DIR \
-i Final.Blake/further/sankararaman12_NGFI_Admix.est \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
--seed \{} \
-g $NCHROM 10 \
--samples_within AfW:0:50:0:YRI EuA:0:1:0:CEU.1 Nea:0:1:50000:Vindija && \
python3.7 $SCRIPT_SUMSTATS \
--bin_dir $BIN_DIR \
-p param_files/stats_noMAC_diplo.spar \
--ld1 \
-i Final.Blake/further/LD_implementation/data/\{}.ld1 \
--ceu CEU.1 --yri YRI --vindija Vindija \
--psmc_pops CEU.1 \
-o Final.Blake/further/LD_implementation/\{}.ld1 \
-v"

parallel --jobs 5 python3.7 $SINGLEIND ::: {1..30}


#===========================================================================#
#T> Estimation F_ST trajectories
#R> local
#===========================================================================#

# Generate log10-spaced sampling times
R --slave -e 'cat(paste0(unique(as.integer(10**seq(0, log10(6926433), length.out=30))), collapse="\n"),"\n")' >Final.Blake/further/TIMES.par
# Manually change time 6926433 to 6500000

# F_ST between AfW:4 and AfE:4
DATA_OUT=Final.Blake/further/Fst_traj_AfW_AfE
mkdir -p $DATA_OUT
CMD="$SCRIPT_SIMULATE \
-o $DATA_OUT/\{} \
--bin_dir ~/bin \
-i $DATA_DIR/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--seed 42 \
--model qian_demo \
-g 10 7 \
--ceu CEU --yri YRI --vindija Vindija \
-s $SCRIPT_SUMSTATS param_files/stats_noMAC_diplo.spar \
--stats stats SE \
--samples_within \
AfW:4:10:\{}:YRI \
AfE:4:10:\{}:CEU \
Nea:4:1:50000:Vindija"
parallel --jobs 7 -a Final.Blake/further/TIMES.par python3.7 $CMD

# F_ST between AfW:0 and AfW:9
DATA_OUT=Final.Blake/further/Fst_traj_AfW0_AfW9
mkdir -p $DATA_OUT
CMD="$SCRIPT_SIMULATE \
-o $DATA_OUT/\{} \
--bin_dir ~/bin \
-i $DATA_DIR/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--seed 42 \
--model qian_demo \
-g 10 7 \
--ceu CEU --yri YRI --vindija Vindija \
-s $SCRIPT_SUMSTATS param_files/stats_noMAC_diplo.spar \
--stats stats SE \
--samples_within \
AfW:0:10:\{}:YRI \
AfW:9:10:\{}:CEU \
Nea:4:1:50000:Vindija"
parallel --jobs 7 -a Final.Blake/further/TIMES.par python3.7 $CMD

# F_ST between AfE:0 and AfE:9
DATA_OUT=Final.Blake/further/Fst_traj_AfE0_AfE9
mkdir -p $DATA_OUT
CMD="$SCRIPT_SIMULATE \
-o $DATA_OUT/\{} \
--bin_dir ~/bin \
-i $DATA_DIR/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--seed 42 \
--model qian_demo \
-g 10 7 \
--ceu CEU --yri YRI --vindija Vindija \
-s $SCRIPT_SUMSTATS param_files/stats_noMAC_diplo.spar \
--stats stats SE \
--samples_within \
AfE:0:10:\{}:YRI \
AfE:9:10:\{}:CEU \
Nea:4:1:50000:Vindija"
parallel --jobs 7 -a Final.Blake/further/TIMES.par python3.7 $CMD

# F_ST between AfE:4 and EuA:4
DATA_OUT=Final.Blake/further/Fst_traj_AfE_EuA
mkdir -p $DATA_OUT
CMD="$SCRIPT_SIMULATE \
-o $DATA_OUT/\{} \
--bin_dir ~/bin \
-i $DATA_DIR/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--seed 42 \
--model qian_demo \
-g 10 7 \
--ceu CEU --yri YRI --vindija Vindija \
-s $SCRIPT_SUMSTATS param_files/stats_noMAC_diplo.spar \
--stats stats SE \
--samples_within \
AfE:4:10:\{}:YRI \
EuA:4:10:\{}:CEU \
Nea:4:1:50000:Vindija"
parallel --jobs 7 -a Final.Blake/further/TIMES.par python3.7 $CMD

#===========================================================================#
#D> 27.01.23
#T> Simulate with empirical recombination maps
#N> Run 1014230
#G> 22x50 Mbp
#R> ocean
#===========================================================================#

#G> Empirical recombination map => HapMap
python3.7 $SCRIPT_SIMULATE \
--write_geno1 \
-o Final.Blake/data/empirical_gmaps/HapMap \
--bin_dir $BIN_DIR \
-i archives/accepted_runs/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
-g 20 50 \
--map data/genetic_maps/HapMap/GRCh37/plink.chrXXX.GRCh37.map \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed 42 \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI

#G> Empirical recombination map => Spence19
python3.7 $SCRIPT_SIMULATE \
--write_geno1 \
-o Final.Blake/data/empirical_gmaps/Spence19 \
--bin_dir $BIN_DIR \
-i archives/accepted_runs/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
-g 20 50 \
--map data/genetic_maps/Spence19/hg19/converted_maps/YRI.chrXXX.map \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar \
--stats stats SE sprime crf psmc ld \
--seed 42 \
--ceu CEU --yri YRI --vindija Vindija --altai Altai \
--psmc_pops CEU YRI

#T> To store the log file containing more info about the resulting mean recombination rate across simulated genomes
#R> local
python3.7 $SCRIPT_SIMULATE \
-x \
-o Final.Blake/data/empirical_gmaps/HapMap. \
--bin_dir ~/bin \
-i archives/accepted_runs/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
-g 20 50 \
--map data/genetic_maps/HapMap/GRCh37/plink.chrXXX.GRCh37.map \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
--seed 42 | tee Final.Blake/data/empirical_gmaps/HapMap..log
#
python3.7 $SCRIPT_SIMULATE \
-x \
-o Final.Blake/data/empirical_gmaps/Spence19. \
--bin_dir ~/bin \
-i archives/accepted_runs/$RUN.par \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model qian_demo \
-g 20 50 \
--map data/genetic_maps/Spence19/hg19/converted_maps/YRI.chrXXX.map \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
--seed 42 | tee Final.Blake/data/empirical_gmaps/Spence19..log


#===========================================================================#
#D> 24.03.2023
#T> Simulate the 20 accepted runs with the Jukes-Cantor 1969 mutation model
#G> 20x30Mbp
#G> JC69 mutation model
#R> ocean
#===========================================================================#

#N> seed=42
CMD="$SCRIPT_SIMULATE \
--model qian_demo \
--write_geno1 \
--bin_dir $BIN_DIR \
-i archives/accepted_runs/\{}.par \
-o Final.Blake/data/20x30Mbp.JC69/\{} \
-g 20 30 --mut_model 'JC69(False)' \
--algo hudson \
--seed 42 \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar --stats stats SE sprime psmc ld crf \
--ceu CEU --yri YRI --vindija Vindija --altai Altai --psmc_pops CEU YRI"
parallel -a Final.Blake/FINAL.accepted --jobs 20 python3.7 $CMD

#G> 10x7Mbp
CMD="$SCRIPT_SIMULATE \
--model qian_demo \
--write_geno1 \
--bin_dir $BIN_DIR \
-i archives/accepted_runs/\{}.par \
-o Final.Blake/data/10x7Mbp.JC69/\{} \
-g 10 7 --mut_model 'JC69(False)' \
--algo hudson \
--seed \{} \
--samples_within AfW:pSample.YRI:50:0:YRI EuA:pSample.CEU:50:0:CEU Nea:pSample.Vindija:1:50000:Vindija Nea:pSample.Altai:1:130000:Altai \
-s $SCRIPT_SUMSTATS param_files/stats_MAC_diplo.spar --stats stats SE sprime psmc ld crf \
--ceu CEU --yri YRI --vindija Vindija --altai Altai --psmc_pops CEU YRI"
parallel -a Final.Blake/FINAL.accepted --jobs 20 python3.7 $CMD

#___
