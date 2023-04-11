#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:59:28 2021
@author: rtournebize

From msp1__zhun.py v1.1.0

0       150421  creation
0.1     160421  removed systematic `quit()` right after history printing
                  added `-x` switch to skip simulation
                  added `--iicr` argument to calculate the IICR using msprime built-in function
                  debugged: it was writing a duplicate of generation_time/seq_length/n_seq when reading *.par file that precontained those
0.2     160421  added `--samples` to specify other sampling strategy + `--stats` to give option to calculate other statistics
0.3     210421  added the possibility to output *.dema1 and *.dema2 files which summarize the demo history for further graphical representations
0.4     260421  added options `--conditions_1` & `--conditions_2` so that we can alter the conditions of pre-rejection externally
                  removed switch `--nolookup` as it is redundant with the previous conditions (which are set as None by default => no check)
                  added the knot/ninsar/thot model calls, for this added a new option `-p` giving the location of an external parameter file
0.4.1   280421  option `-s` is not mandatory anymore; if not provided, no check will be made
0.4.2   020521  minor changes on the way modules are imported (now via importlib) + call to module extras.py  + added `--version`
0.4.3   280521  in sample(): the deme ID was read as str, but should be int to be altered if `deme_id_is_1_based==True`
0.5     020621  option `--map` to provide a genetic map instead of assuming a uniform rec rate with the --rec_rate option
0.5.1   040621  now exports a file "*.samples" containing all haploid individuals that are sampled, if `--verbose`
0.5.2   070621  added option `--draw` to draw the demography using the DemesDraw module
0.5.3   010721  added option `--bin_dir` that is now required for the link with sumstats.py to provide the path of external softwares/modules
0.6.0   020821  option `--samples_float` which can replace `--sample` to specify the sampling deme relatievely, i.e. the sampled deme will be located at int(p*specific_chain_length)
                e.g. --sample_relative 0.4:0:10:0:CEU means: sample at position 0.4 of the chain 0 a total of 10 diploids at age 0 gBP and label those as "CEU"
0.6.1   050921  added a condition so that statistic calculation is not performed if --stats is not provided and default for --stats is now None
0.6.2   120921  if --map, writes -1 for seq_length in the *.par file
0.6.3   221021  added possibility to output the ms command of the demography through option `--ms`

`msp1__zhun` (requires python >=3.7)
1.0.0   101121  first implementation with msprime >=1.0
                removed --samples_float and corresponding sample_relative() function
                removed --iicr
                added --algo option to specify which kind of msp1 algorithm (eg Hudson, Beta, ...)
                added --gconv to introduce gene conversion events
                added --lin_prob_times to specify the parameters of linear sampling for the times of lineage prob calculation
                now writes in *.gz format directly
                removed dependency to {extras.py}
                rewrote write_eigenstrat() to make it quicker
                introduced a seeding option for the simulations
                the rules are now specified in the *.est file and start by the "#R>" code
                now possible to do complex parameter calculation in the *.est file, the formula should start with the "#F>" code
1.0.1   171121  now outputs the ms command in a *.ms file when --draw is present
1.0.2   261121  new option `--samples_within` which allows to specify the deme to sample in according to its position within its own chain
1.0.3   261121  `--samples_within` can now accept relative positioning of the deme within its chain (ie. 0 <= f <= 1)
                also forced reading and transforming the deme sizes (*.n) as integers
1.1.0   071221  reimplemented the import and simulation with recombination maps (RateMaps) using the msprime >1.0 syntax
1.2.0   081221  reimplemented the import of samples:
                    - removed the option `--samples`, now only: `--samples_within` is available
                    - assumes everything (incl input sample sizes) are ***DIPLOIDS***
1.3.0   091221  now outputs the par file by respecting the parameter order (= order of inclusion) since the ordering is super important
1.3.1   141221  small bug when a parameter is defined with "math" -- surnumerary *.split(" ")
1.3.2   181221  added option `--write_geno12` to output a *.geno12.gz file containing all haplotypes
1.3.3   201221  now can provide a variable name for the float proportion of sampling in --sample_within
1.3.4   221221  implemented possibility to provide multiple recombination maps, therefore using hybrid demo models
1.4.0   221221  now through `--mutation_model`, possiblity to implement other mutation models -- consequently had to change the write_eigenstrat() function
1.4.1   231221  in case of multiple recombination maps, determine the map boundaries based on the intersection of all maps (so that the maps are remained aligned, in the previous version there would be map-dependent mismatches)
                renamed default standard algo from "kingman" into "hudson" (as it is) + added model 'dtwf'
1.4.2   140122  added possibility to filter also on `SPRIME.length.mean`, `SPRIME.match.mean` and `DCFS.1`
1.4.3   100322  added switch `--check` which will check that the coalescent process works fine; if not, will throw a *.failed.par
1.4.4   110322  added switch `--force_gpos` (False by default) to force writing the genetic positions in the *.snp file
1.4.5   210322  corrected small bug when using math distrib (#F> switch in par file): missing a split operation on the formula string
                more precision in the help message of --map about the format of the hapmap rate map file
2.0.0   250322  implemented the way to specify recombination rates with crossing-over (as in msHOT, Hellenthal et al 08): options
                --do_hotspots (switch); --hotspots_pars (list of par values); --hostpots_seed (seeding)
                also changed the way seed was implemented, now more rigorous
3.0.0   270422  reorganized main() with new external functions
3.0.1   270422  redesigned the dict() of Maps in get_maps() with proper str keys + changed the returned args--subsequent changes in main() and simulate() + proper ordering of the Maps times
3.0.2   280422  introduced options --ceu, --yri, --vindija, --altai, psmc_pops for continuity with sumstats.py >14.0
3.1.0   290422  implemented seeding for the parameter sampling too /!\
                checked 3.1.0 results vs. 1.4.5
3.1.1   290422  bug with lookup(): since I added sumstats version in the head of *.stats (v3.1.0), that shifted line numbering by 1--pb solved by excluding commented line reading
3.2.0   280422  introduced possibility to run other (eg published) models in a more flexible way
                print version of packages
3.2.1   290422  added option --benchmark
                debug export_all() when plotting demo models of non qian models
4.0.0   030422  implemented a fast EIGENSTRAT writing in case where mutation model imposes biallelic SNPs
4.1.0   281122  new option `--dump_ts` to write the binary file of the TreeSequence objects
                added msprime version check
                changed the location of the SampleSet printing (moved out from the model script)
4.1.1   051222  added option `--no_simplify` to not simplify the TreeSequence object (for debugging purposes)
4.2.0   061222  included possibility to override mut_rate and rec_rate using parameters in the *.par file
4.2.1   121222  changed option `--no_simplify` to `--simplify`, meaning that compared with previous versions, by default from now on,
                the TreeSequences will *NOT* be simplified with ts.simplify(), as simplification can lead to
                node ID changes, which could create confusion + the speedup time brought by simplification is insignificant
4.2.2   211222  more specific error message if call to model.history() failed + check that sampling refers to existing pops in demo model +
                added "Error." as prefix of all sys.exit messages
4.3.0   291222  added `--assume_no_na` option to the call to sumstats.py, as default
4.3.1   040123  introduced an option `--max_seed` to set the upper bound for the overall seeding (should be set as 1e10 to avoid duplicated seeds)
4.3.2   270123  minor bug when importing empirical genetic maps: during truncation, was using index instead of key value for slicing dict of MapRanges
4.3.3   200323  possibility to output different image format for the model
4.3.4   220323  check that only biallelic SNPs in output data


*** Parameters ***

Possible sampling distributions are:
    - unif a b                uniform from a to b
    - log10unif               log10-uniform from a to b
    - seq a b s precision     linear sampling from a to b by step of s
    - log10seq a b s prec     log-10 linear sampling from a to b by step of log10(s)
    - discrete                randomly sample one value out of a user-specified vector of values
    - [one specific number]   use this very number

*** Mutation & recombination rates as parameter
To override --mut_rate and --rec_rate using parameter-style rates, you should write in your *.par or *.est par input file:
    mutation_rate   VALUE_or_DISTRIB_LAW
    recombination_rate   VALUE_or_DISTRIB_LAW
This will replace the values provided as options.

*** Rules ***

To be specified in the *.est file, one rule per line, starting with:     #R>

*** Parameter operations / formula / ratio ***

To be specified in the *.est file, one calculation per line, starting with:     #F>
    Example:
        #F> n.nea.AfW f.0.AfW * X0
    Will be understood as n.nea.AfW = f.0.AfW * X0

*** Package requirements ***

For demes-python | demes:
    python3.7 -m pip install git+https://github.com/popsim-consortium/demes-python.git
For stdpopsim:
    python3.7 -m pip install git+https://github.com/popsim-consortium/stdpopsim.git@190362c

"""

___version___ = '4.3.4'

import sys, argparse, os, subprocess, yaml, time, gzip, copy
import pandas as pd
import numpy as np
import msprime as msp
import demes, stdpopsim
from numpy.random import Generator, PCG64


############################################################
############################################################
############################################################
############################################################
############################################################

def main():
    init_time = time.time()
    if not msp.__version__.startswith("1."):
        sys.exit("Error. msprime version must be v1")
    RunningTime = {"A.init": init_time}
    ARGS, max_iter = get_args()
    # seeding for params
    par_seeds = [None] * int(max_iter)
    if ARGS["seed"] is not None:
        rg0 = Generator(PCG64(ARGS["seed"]))
        par_seeds = rg0.uniform(low = 1, high = ARGS["max_seed"], size = int(max_iter)).astype(np.int64)
    # rules
    PAR = get_pars(ARGS["input"], par_seeds[0])
    Rules = []
    with open(ARGS["input"], "r") as F:
        for line in F:
            if line.startswith("#R>"):
                X = line.strip().replace("#R>", "").split()
                Rules.append(" ".join([x.replace(x, "Par[\""+x+"\"]") if x in set(PAR.keys()) else x for x in X]))
    # parameters
    RunningTime["A1.read"] = time.time()
    is_valid, it = validate(PAR, Rules), 0
    while is_valid==False and it<max_iter:
        PAR = get_pars(ARGS["input"], par_seeds[it])
        is_valid = validate(PAR, Rules)
        it += 1
    if it==max_iter:
        sys.exit("Error. Could not find proper parameter space.")
    if ARGS["output"]+".par" == ARGS["input"]:
        sys.exit("Error. Output par file is identical to input par file.")
    # demography history
    RunningTime["A2.pars"] = time.time()
    model = __import__(ARGS["model"])
    try:
        ARGS, Demography, Samples, PopLabels = model.history(ARGS, PAR)
        if "mutation_rate" in PAR.keys():
            print("\n_/!\__ Warning. Erasing option-given `mut_rate` by a parameter-given one: "+str(PAR["mutation_rate"])+"\n")
            ARGS["mut_rate"] = PAR["mutation_rate"]
        if "recombination_rate" in PAR.keys():
            print("\n_/!\__ Warning. Erasing option-given `recomb_rate` by a parameter-given one: "+str(PAR["recombination_rate"])+"\n")
            ARGS["rec_rate"] = PAR["recombination_rate"]
    except:
        os.rename(ARGS["output"]+".par", ARGS["output"]+".failed.par")
        sys.exit("Error. Faulty demographic model or sampling.")
    # check sampling
    SampledPops = np.unique(np.array([x.population for x in Samples]))
    SamplingOK = np.isin(SampledPops, np.array([x.name for x in Demography.populations]))
    if np.sum(SamplingOK==False) > 0:
        sys.exit("Error. These populations you wanted to sample in do not exist in the demographic model:   "+", ".join(SampledPops[SamplingOK==False]))
    # case of genetic map
    RunningTime["A3.demo_import"] = time.time()
    Maps = get_maps(ARGS)
    # write
    RunningTime["A4.get_map"] = time.time()
    ARGS_PRINT = copy.deepcopy(ARGS)
    ARGS_PRINT["max_seed"] = int(ARGS_PRINT["max_seed"])
    print(yaml.dump(ARGS_PRINT))
    print("Samples.")
    print(Samples)
    print("\nRules.")
    for rule in Rules:
        print("- "+rule.replace("Par[\"", "").replace("\"]", ""))
    with open(ARGS["output"]+".par", "w") as FOUT:
        sl = ARGS["genome"][1]
        if sl > 0:
            sl *= 1e6
        FOUT.write("#qian: "+str(___version___)+"\n")
        FOUT.write("seq_length "+str(int(sl))+"\nn_seq "+str(ARGS["genome"][0])+"\nmut_rate "+str(ARGS["mut_rate"])+"\nrecomb_rate " +str(ARGS["rec_rate"])+"\ngeneration_time "+str(ARGS["generation_time"])+"\n")
        for key in PAR.keys():
            FOUT.write(key+" "+str(PAR[key])+"\n")
    # additional parameters
    PAR2 = None
    if ARGS["param_file"] is not None:
        PAR2 = read_additional_param_file(ARGS["param_file"])
        print("")
        print(yaml.dump(PAR2))
    # export demography
    RunningTime["A5.io_pars"] = time.time()
    if ARGS["draw"] is not None:
        print("\nDrawing demography.")
        export_all(ARGS, PAR, Demography, Samples)
    # simulate (or not)
    RunningTime["B.exportdemo"] = time.time()
    if ARGS["nosimul"]:
        print("Skip simulation.\n")
        quit()
    else:
        simulate(ARGS,
             Demography,
             Samples,
             PopLabels,
             Maps)
    # calculate final statistics
    RunningTime["C.simulatewrite"] = time.time()
    if ARGS["stat_command"][0] is not None and ARGS["stats"] is not None:
        calculate(ARGS)
    # clean
    RunningTime["D.calculate"] = time.time()
    if ARGS["remove_dataset"]:
        os.remove(ARGS["output"]+".geno.gz")
        if ARGS["write_geno1"]:
            os.remove(ARGS["output"]+".geno1.gz")
        os.remove(ARGS["output"]+".snp.gz")
        os.remove(ARGS["output"]+".ind")
    # end
    RunningTime["E.end"] = time.time()
    print(str(np.round((time.time()-init_time)/60, 1))+" min.\n")
    if ARGS["benchmark"]:
        keys_RT = sorted(RunningTime.keys())
        for ikey, key in enumerate(keys_RT):
            if ikey == 0:
                print(key+": "+str(RunningTime[key]))
                continue
            print(key+": "+str(RunningTime[key])+" | "+"%.2f"%( ((RunningTime[keys_RT[ikey]]-RunningTime[keys_RT[ikey-1]])/60.) )+"m")

####################################################################################
####################################################################################

def calculate(ARGS):
    print("\nCalculate statistics.")
    cmd_stat = "python3 "+str(ARGS["stat_command"][0])+" -p "+str(ARGS["stat_command"][1])+" -i "+ARGS["output"]+" --bin_dir "+ARGS["bin_dir"]+" --assume_no_na"
    cmd_stat = cmd_stat+" --"+" --".join(ARGS["stats"])
    for opt in ["ceu", "yri", "vindija", "altai"]:
        if ARGS[opt] != None:
            cmd_stat = cmd_stat+" --"+opt+" "+ARGS[opt]
    if ARGS["psmc_pops"] is not None:
        cmd_stat = cmd_stat+" --psmc_pops "+" ".join(ARGS["psmc_pops"])
    subprocess.run(cmd_stat, shell=True, check=True, stdout=subprocess.DEVNULL)

def simulate(ARGS,
             Demography,
             Samples,
             PopLabels,
             Maps = dict()):
    assert type(Maps) == type(dict()), "Maps should be a dictionary (empty, if no genetic maps provided)"
    # algorithm
    algo = [None, None]
    if ARGS["algo"][0] == "hudson":
        algo = msp.StandardCoalescent()
    elif ARGS["algo"][0] == "dtwf":
        algo = msp.DiscreteTimeWrightFisher()
    elif ARGS["algo"][0] == "beta":
        algo = msp.BetaCoalescent(alpha = float(ARGS["algo"][1]))
    elif ARGS["algo"][0] == "smc":
        algo = msp.SmcPrimeApproxCoalescent()
    else:
        sys.exit("Error. Algorithm",algo,"not implemented.")
    # seeding
    if ARGS["seed"] is None:
        seeds = [None] * int(ARGS["genome"][0])
    else:
        rg0 = Generator(PCG64(ARGS["seed"]))
        seeds = rg0.uniform(low = 1, high = ARGS["max_seed"], size = int(ARGS["genome"][0])).astype(np.int64)
    # sumstats cmd
    cmd_stat = "python3 "+str(ARGS["stat_command"][0])+" -p "+str(ARGS["stat_command"][1])+" -i "+ARGS["output"]+" --bin_dir "+ARGS["bin_dir"]+" --assume_no_na"
    for opt in ["ceu", "yri", "vindija", "altai"]:
        if ARGS[opt] != None:
            cmd_stat = cmd_stat+" --"+opt+" "+ARGS[opt]
    if ARGS["psmc_pops"] is not None:
        assert type(ARGS["psmc_pops"]) == type([]), "--psmc_pops option must be a list of strings"
        cmd_stat = cmd_stat+" --psmc_pops "+" ".join(ARGS["psmc_pops"])
    # simulate
    print("\nSimulation.")
    for chrom in range(int(ARGS["genome"][0])):
        print("\n===================>  "+str(chrom+1))
        ts, n_sim = None, len(Maps)
        if n_sim == 0:
            n_sim = 1
        for sa in range(n_sim):
            if sa == 0:
                initial_state = None
                start_time = None
                if len(Maps) > 1:
                    end_time = float(Maps.keys()[1])
                else:
                    end_time = None
                if len(Maps) > 0:
                    Map = Maps["0"]
                samples = Samples
            else:
                initial_state = ts
                start_time = float(Maps.keys()[sa])
                if sa == n_sim-1:
                    end_time = None
                else:
                    end_time = float(Maps.keys()[sa+1])
                Map = Maps[sa]
                samples = None
            try:
                ts = msp.sim_ancestry(samples=                         samples,
                                      demography=                      Demography,
                                      sequence_length=                 ARGS["genome"][1]*1e6 if (len(Maps)==0) else None,
                                      discrete_genome=                 True,
                                      recombination_rate=              float(ARGS["rec_rate"]) if (len(Maps)==0) else Map[chrom],
                                      gene_conversion_rate=            float(ARGS["gconv"][0]) if ARGS["gconv"][0] is not None else None,
                                      gene_conversion_tract_length=    int(ARGS["gconv"][1]) if ARGS["gconv"][1] is not None else None,
                                      population_size=                 None,
                                      ploidy=                          2,
                                      model=                           algo,
                                      initial_state=                   initial_state,
                                      start_time=                      start_time,
                                      end_time=                        end_time,
                                      record_migrations=               False,
                                      record_full_arg=                 False,
                                      num_labels=                      None,
                                      random_seed=                     seeds[chrom],
                                      num_replicates=                  None,
                                      replicate_index=                 None,
                                      record_provenance=               False)
                msg = "Coalescent check: OK."
            except:
                os.rename(ARGS["output"]+".par", ARGS["output"]+".failed.par")
                sys.exit("Error. Faulty coalescent process.")
            if ARGS["check"]:
                sys.exit(msg)
        if ARGS["dump_ts"]:
            out_ts = ARGS["output"]+"."+str(chrom+1)+".ts"
            if os.path.exists(out_ts) and os.path.isfile(out_ts):
                os.remove(out_ts)
            ts.dump(out_ts)
            if chrom == 0:
                if len(PopLabels) != len(ts.samples()):
                    sys.exit("Error. Incompatibility between PopLabels range and the TreeSequence individual list.")
                seq_along_PopLabels = np.array(list(range( len(PopLabels) )))
                if (seq_along_PopLabels == ts.samples()).all() == False:
                    sys.exit("Error. Incompatibility between PopLabels range and the TreeSequence individual list.")
                with open(ARGS["output"]+".inds.ts", "wt") as TSI:
                    for ts_ind in range(len(PopLabels)):
                        TSI.write(str(ts_ind)+" "+str(PopLabels[ts_ind])+"\n")
        if ARGS["simplify"] == True:
            ts = ts.simplify()
        ts = msp.sim_mutations(ts,
                               rate=                           float(ARGS["mut_rate"]),
                               model=                          eval("msp."+ARGS["mut_model"]),
                               random_seed=                    seeds[chrom])
        if ARGS["force_gpos"]:
            do_write_gpos = True
        else:
            if len(Maps)==0:
                do_write_gpos = False
            else:
                do_write_gpos = True
        if "BinaryMutationModel" in ARGS["mut_model"]:
            assume_biallelic_only = True
        else:
            assume_biallelic_only = False
        write_eigenstrat( treeseq = ts,
                          prefix = ARGS["output"],
                          chromosome_id = chrom+1,
                          pop_labels = PopLabels,
                          diploidize = True,
                          assume_biallelic_only = assume_biallelic_only,
                          recombination_rate = ARGS["rec_rate"] if (len(Maps)==0) else None,
                          recombination_map = None if (len(Maps)==0) else Map[chrom],
                          append = chrom>0,
                          write_geno1 = ARGS["write_geno1"],
                          snp_write_gpos = do_write_gpos,
                          write_geno12 = ARGS["write_geno12"] )
        if (ARGS["stat_command"][0] is not None) and (ARGS["conditions_1"] is not None) and (chrom+1 == 1):
            lookup(cmd_stat, ARGS["output"], ARGS["conditions_1"], ARGS["remove_dataset"], ARGS["write_geno1"])
        if (ARGS["stat_command"][0] is not None) and (ARGS["conditions_2"] is not None) and (chrom+1 == 2):
            lookup(cmd_stat, ARGS["output"], ARGS["conditions_2"], ARGS["remove_dataset"], ARGS["write_geno1"])

def get_args():
    # fixed arguments
    max_iter = 200
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version='%(prog)s: '+str(___version___))
    ## mandatory
    parser.add_argument('-i', '--input', type=str, required=True, help="File of demographic parameters")
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-g', '--genome', type=float, nargs=2, required=True, help="List: [n_seq] [seq_length_Mbp]. If --map is provided, set the sequence length to -1 if you want to simulate the whole chromosome length provided in the genetic map files.")
    parser.add_argument('--model', type=str, required=True)

    ## optional
    parser.add_argument('--mut_rate', type=float, required=False, default=1.2e-8)
    parser.add_argument('--rec_rate', type=float, required=False, default=1.0e-8)
    parser.add_argument('--map', type=str, required=False, default=[None], nargs="*", help="Path to the genetic map files, HapMap format without header, genetic positions in cM. CAUTION! Order of columns: chromosome, Rate(cM/Mbp)[can_be_empty], GeneticPosition(cM), PhysicalPosition(bp). In the file name provided, chromosome number must be replaced by XXX. If the sequence length provided in --genome is not set to -1, the resulting sequence length will be sliced to the specified length. YOU CAN PROVIDE MULTIPLE GENETIC MAPS. In this case, starting from the second map, you should write TIME_IN_ya_WHEN_MAP_STARTS_BACKWARD:PATH_TO_MAPS_chrXXX")
    parser.add_argument('--generation_time', type=float, required=False, default=25)

    parser.add_argument('-s', '--stat_command', type=str, nargs=2, required=False, default=[None, None], help="List: [script_location] [stat_par_file_location]")
    parser.add_argument('--stats', type=str, nargs="*", default=None, help="List, e.g.: stats ld sprime")
    parser.add_argument('--conditions_1', type=str, nargs="*", default=None, help="A list of conditions used to pre-reject the run after simulating the first chromosome. The conditions are formed like, e.g. Fst:min_value:max_value pi.YRI:min_value:max_value")
    parser.add_argument('--conditions_2', type=str, nargs="*", default=None, help="A list of conditions used to pre-reject the run after simulating the second chromosome. The conditions are formed like, e.g. Fst:min_value:max_value pi.YRI:min_value:max_value")

    parser.add_argument('-v', '--verbose', action="store_true", default=False, help="Generate *.dema, *.dema1 and *.dema2 files containing descriptors of the demographic history: *.dema is the msprime .print_history output, *.dema1 is a table of events, *.dema2 is a table of migration changes")
    parser.add_argument('-x', '--nosimul', action="store_true", default=False, help="Skip the simulation process")

    parser.add_argument('-p', '--param_file', type=str, required=False, default=None, help="Location of an external file containing additional parameters for plotting the model")
    parser.add_argument('--write_geno1', action="store_true", required=False, default=False, help="Add this switch to write the haplotype *.geno1 file")
    parser.add_argument('--write_geno12', action="store_true", required=False, default=False, help="Add this switch to write a *.geno12 file containing all the haplotypes")
    parser.add_argument('--remove_dataset', action="store_true", default=False, help="Remove the genetic dataset if the run succeeded")
    parser.add_argument('--samples_within', type=str, nargs="*", default=None, help="List of 5: [chain_name:deme_index_0bases_within_own_chain_to_sample__OR__relative_position_of_deme_within_chain_between_0_and_1:diploid_sample_size:sample_age_gBP:label_of_samples]. NOTE THAT TO BE INTERPRETED AS A RELATIVE POSITION, THE NUMBER MUST HAVE A POINT, i.e. `0.` or `1.` for the bounds.")
    parser.add_argument('--samples', type=str, nargs="*", default=None, help="List of 4: [name_of_deme_to_sample_in:diploid_sample_size:sample_age_gBP:label_of_samples].")
    parser.add_argument('--draw', type=str, default=None, help="Add this switch with a parameter file path to draw the demographic model using Demes and DemesDraw modules, as well as to export the demographic model in *.demo and *.yaml formats.")
    parser.add_argument('--bin_dir', type=str, default="/Users/rtournebize/bin", help="Path to the directory containing external modules/softwares used in sumstats.py")
    parser.add_argument('--ms', action="store_true", required=False, default=False, help="Add this switch to output a *.ms file containing the ms command associated with the demographic scenario")
    parser.add_argument('--algo', type=str, required=False, default=["hudson"], nargs="*", help="Type of algorithm from msprime that is to be used (default: hudson), can also be: smc, dtwf, beta, in that case also specifies the value of 1<alpha<2")
    parser.add_argument('--gconv', type=float, required=False, default=[None, None], nargs=2, help="Parameters of the gene conversion process: gene conversion rate, mean gene conversion tract length in bp")
    parser.add_argument('--lin_prob_times', type=float, required=False, default=[0., 1e6, 1e4], nargs=3, help="Specify the min, max and step value for ages (in ya) to linearly sample for exporting the lineage probability across demes")
    parser.add_argument('--seed', type=int, required=False, default=None)
    parser.add_argument('--mut_model', type=str, required=False, default="BinaryMutationModel(False)", help="Mutation model")
    parser.add_argument('--check', action="store_true", required=False, default=False, help="Add this switch to check if the coalescent genealogical simulation works fine; if true, will write a proper *.par file but will not execute the code further, if not, will write a *.failed.par file.")
    parser.add_argument('--force_gpos', action="store_true", required=False, default=False, help="By default, genetic positions are replaced by . when using a mean recombination rate through --rec_rate option. Add this switch to force writing the genetic positions in Morgans.")
    ###
    parser.add_argument('--do_hotspots', action="store_true", required=False, default=False, help="Add this switch to simulate artificial recombination rates with crossover hotspots as in msHOT - Hellenthal et al 08. If you specify this option, the value of --rec_rate will be ignored.")
    parser.add_argument('--hotspots_seed', type=int, required=False, default=42, help="Seed for the simulation of recombination maps with crossover hotspots.")
    parser.add_argument('--hotspots_pars', type=float, nargs=7, required=False, default=[0.2325e-8, 40, 5, 1, 2, 1, 2.5], help="Parameter values for the simulation of recombination maps with crossover hotspots. There are 7 parameters: [rec_rate_outside_hotspots, average_distance_between_hotspots_kb, min_spacing_between_hotspot_regions_kb, lower_bound_hotspot_region_width_kb, upper_bound_hotspot_region_width_kb, lower_bound_log10lambda, upper_bound_log10lambda]")
    # Sample specification for sumstats.py
    parser.add_argument('--ceu', type=str, required=False, help="Passed to sumstats: CEU population label, as in the *.ind file; alternatively, in the param file: CEU_indices_0based", default = None)
    parser.add_argument('--yri', type=str, required=False, help="Passed to sumstats: YRI population label, as in the *.ind file; alternatively, in the param file: YRI_indices_0based", default = None)
    parser.add_argument('--vindija', type=str, required=False, help="Passed to sumstats: Vindija_Neanderthal population label, as in the *.ind file; alternatively, in the param file: Vindija_indices_0based", default = None)
    parser.add_argument('--altai', type=str, required=False, help="Passed to sumstats: Altai_Neanderthal population label, as in the *.ind file; alternatively, in the param file: Altai_indices_0based", default = None)
    parser.add_argument('--psmc_pops', type=str, required=False, nargs="*", help="Passed to sumstats: Labels of the populations to analyze with PSMC, will select the first ind of each pop; alternatively, in the param file: psmc_samples_0based", default = None)
    parser.add_argument('--benchmark', action="store_true", help="Add this switch if you want to print a summary of running time by task.")
    parser.add_argument('--dump_ts', action="store_true", help="Add this switch to write the TreeSequence, will create output.X.ts objects where X corresponds to each chromosome")
    parser.add_argument('--simplify', action="store_true", help="Add this switch to simplify the TreeSequence in order to (possibly) speed up the script execution")
    parser.add_argument('--max_seed', type=np.int64, nargs=1, required=False, default=np.int64(1e6), help="Maximum allowed seed, should be high enough, on the order of 1e10") # 4.3.1
    # store arguments in a dictionary, all values will be string by default
    ARGS = vars(parser.parse_args())
    print("\n")
    print("                 __  __")
    print("                 __  __")
    print("                 __  __")
    print("                 ______")
    print("                 __  __")
    print("                 __  __\n")
    print("v"+str(___version___)+"\n")
    # versions
    Versions = "msprime: "+str(msp.__version__)+" | "+"stdpopsim: "+str(stdpopsim.__version__)+" | "+"demes: "+str(demes.__version__)
    if ARGS["draw"] is not None:
        import demesdraw
        Versions += " | demesdraw: "+str(demesdraw.__version__)
    print(Versions+"\n")
    # create path
    root = os.path.dirname(ARGS["output"])
    if os.path.exists(root)==False and root!="":
        os.mkdir(root)
    # checks
    if ARGS["map"][0] is not None and ARGS["do_hotspots"] is True:
        sys.exit("Error. You cannot specify --map and --do_hotspots altogether.")
    if ARGS["samples_within"] is not None and ARGS["samples"] is not None:
        sys.exit("Error. You cannot specify altogether --samples and --samples_within.")
    if ARGS["samples_within"] is None and ARGS["samples"] is None:
        sys.exit("Error. Samples are missing.")
    # transform
    ARGS["genome"][0] = int(ARGS["genome"][0])
    ARGS["genome"][1] = int(ARGS["genome"][1])
    # if just checking
    if ARGS["check"] == True:
        ARGS["genome"][0] = 1
        ARGS["genome"][1] = float(1e-6)
        ARGS["map"][0] = None
        ARGS["do_hotspots"] = False
        ARGS["rec_rate"] = 1e-8
    return ARGS, max_iter

def get_maps(ARGS):
    #__ if reading from files
    Maps, MapRanges = dict(), dict()
    if ARGS["map"][0] is not None:
        assert ":" not in ARGS["map"][0], "The first map cannot have a time value, because by default it starts at present-time."
        if len(ARGS["map"]) > 0:
            assert np.all([":" in x for x in ARGS["map"][1:]]), "Ancient maps must have a time value associated to them, i.e. TIME:PATH."
        for gmap in ARGS["map"]:
            if ":" in gmap:
                tmap, gmap = gmap.split(":")
            else:
                tmap = 0.
            tmap = int(np.true_divide(float(tmap), ARGS["generation_time"]))
            print("\n> Map starting "+str(tmap)+" ga (backward in time) with path: "+gmap)
            Maps[str(tmap)], MapRanges[str(tmap)] = get_map(ARGS, gmap)
        # intersect maps
        Ranges = intersect_maps(ARGS, MapRanges)
        # slice maps
        Maps, Mean_Rhos = slice_maps(Maps, Ranges)
        mean_rho = np.mean(np.array(Mean_Rhos))
        print("\n\nMean rate: "+str(mean_rho))
        ARGS["rec_rate"] = float(mean_rho)
        print("\nThe default recombination rate will be replaced by the empirical rate averaged over the maps: "+str(ARGS["rec_rate"])+"\n")
        print("\n\n")
    #__ if generating artificial map with crossover hotspots
    elif ARGS["do_hotspots"]:
        TempMaps, mean_rho = [], 0.
        rg0 = Generator(PCG64(np.int64(ARGS["hotspots_seed"])))
        seeds = rg0.uniform(low = 1, high = ARGS["max_seed"], size = int(ARGS["genome"][0])).astype(np.int64)
        for q in range(ARGS["genome"][0]):
            TempMaps += [msHOT_Map_SingleChrom(   seeds[q],
                                                  SeqLen_kb =           int(float(ARGS["genome"][1])*1e3),
                                                  RecRate_OutHotspot =  ARGS["hotspots_pars"][0],
                                                  InterSpacing_kb =     ARGS["hotspots_pars"][1],
                                                  MinSpacing_kb =       ARGS["hotspots_pars"][2],
                                                  Width_kb =            [ARGS["hotspots_pars"][3], ARGS["hotspots_pars"][4]],
                                                  Log10Lambdas =        [ARGS["hotspots_pars"][5], ARGS["hotspots_pars"][6]]
                                                  )]
            mean_rho += TempMaps[-1].mean_rate
            Maps["0"] = TempMaps
        mean_rho = np.true_divide(mean_rho, float(ARGS["genome"][0]))
        print("\n\nMean rate: "+str(mean_rho))
        ARGS["rec_rate"] = float(mean_rho)
        print("\nThe default recombination rate will be replaced by the empirical rate averaged over the maps: "+str(ARGS["rec_rate"])+"\n")
        print("\n\n")
    ###
    if len(Maps) > 0:
        MapTimes = np.sort( np.array(list(Maps.keys())).astype(int) ).astype(str)
        Maps = {key: Maps[key] for key in MapTimes}
    # Maps = dict() (empty) if no gMaps provided as request
    return Maps

def U(a, b):
    return np.random.uniform(a, b)

def LU(a, b):
    return 10**np.random.uniform(np.log10(a), np.log10(b))

def SEQ(lb, ub, step, log10):
    if log10:
        if step <= 1:
            sys.exit("Error. step must be >1")
        x = np.arange(np.log10(lb), np.log10(ub), step = np.log10(step))
    else:
        x = np.arange(lb, ub, step = step)
    if len(x) == 0:
        sys.exit("Error. empty seq vector")
    if log10:
        rx = 10**x
    else:
        rx = step * np.round(np.true_divide(x, step))
        rx[0] = lb
    return np.random.choice(rx, 1)[0]

def get_map(ARGS, gmap_root):
    assert "XXX" in gmap_root, "XXX must replace the chromosome value in the genetic map file name, --map option."
    Map, MapRanges = [], []
    print("q start end length")
    for q in range(int(ARGS["genome"][0])):
        chrom = q + 1
        gfile = gmap_root.replace("XXX", str(chrom))
        gmap = msp.RateMap.read_hapmap(gfile,
                           sequence_length = None,
                           has_header = False,
                           position_col = 3,
                           rate_col = None,
                           map_col = 2)
        bad_int = gmap.missing_intervals()
        start, end = 0., gmap.sequence_length
        if bad_int.shape[0] == 0: pass
        elif bad_int.shape[0] == 1:
            if bad_int[0,0] == 0.: start = bad_int[0,1] + 1
            elif bad_int[0,1] == end: end = bad_int[0,0] - 1
            else: sys.exit("Error. Wrong genetic map.")
        elif bad_int.shape[0] == 2:
            if bad_int[0,0] == 0.: start = bad_int[0,1] + 1
            else: sys.exit("Error. Wrong genetic map.")
            if bad_int[1,1] == end: end = bad_int[1,0] - 1
            else: sys.exit("Error. Wrong genetic map.")
        else:
            sys.exit("Error. Too many missing intervals for: "+gfile+"\nThere must be less or exactly 2 corresponding to the two ends of the chromosome. Check the HapMap file.")
        MapRanges.append([start, end])
        print(str(chrom)+" "+"%.2f"%(start*1e-6)+" "+"%.2f"%(end*1e-6)+" "+"%.2f"%((end-start)*1e-6))
        Map.append(gmap)
    MapRanges = np.array(MapRanges)
    return Map, MapRanges

def intersect_maps(ARGS, MapRanges):
    Range = []
    print("\nIntersected ranges:")
    print("start end length | start end length [keeping only the user-defined length]")
    for q in range(MapRanges["0"].shape[0]):
        rg = [MapRanges[tmap][q,:] for tmap in MapRanges.keys()]
        rg = np.array(rg)
        rg = np.array([np.max(rg[:,0]), np.min(rg[:,1])])
        pr = str(q+1)+" "+"%.2f"%(rg[0]*1e-6)+" "+"%.2f"%(rg[1]*1e-6)+" "+"%.2f"%((rg[1]-rg[0])*1e-6)
        if rg[0] > rg[1]:
            sys.exit("Error. Not intersectable ranges for chromosome "+str(q))
        if float(ARGS["genome"][1]) > 0.:
            assert (rg[0] + float(ARGS["genome"][1])*1e6) < rg[1], "The required length of the sequence leads to a longer chromosome compared to the genetic map."
            rg[1] = rg[0] + float(ARGS["genome"][1])*1e6
        pr = pr+" | "+"%.2f"%(rg[0]*1e-6)+" "+"%.2f"%(rg[1]*1e-6)+" "+"%.2f"%((rg[1]-rg[0])*1e-6)
        print(pr)
        Range.append(rg)
    Range = np.array(Range)
    return Range

def slice_maps(Maps, MapRanges):
    Mean_Rhos = []
    for key in Maps.keys():
        total_lengths = 0
        for q in range(MapRanges.shape[0]):
            Maps[key][q] = Maps[key][q].slice(left = MapRanges[q,0], right = MapRanges[q,1], trim = True)
            total_lengths = total_lengths + np.array([Maps[key][q].sequence_length, Maps[key][q].get_cumulative_mass(Maps[key][q].sequence_length)])
        mean_rho = np.true_divide(total_lengths[1], total_lengths[0])
        Mean_Rhos.append(mean_rho)
    return Maps, Mean_Rhos

def get_pars(input, seed):
    #df = pd.read_csv(input, sep="\t+", header=None, names=["name", "distrib", "lb", "ub"], dtype=str, engine="python", skipinitialspace=True, na_filter=False, skip_blank_lines=True, comment="#").to_numpy()
    df = []
    with open(input, "r") as F:
        for line in F:
            if line.startswith("#F>"):
                line = line.replace("#F>", "").strip().replace("\t", " ")
                line = " ".join(line.split())
                line = line.split(" ")
                line = [line[0], "math", " ".join(line[1:])]
                df.append(line)
            elif line.startswith("#"):
                 continue
            elif line.strip()=="":
               continue
            else:
                line = line.strip().replace("\t", " ")
                line = " ".join(line.split())
                line = line.split(" ")
                df.append(line)
    ###
    Par = dict()
    for i in range(len(df)):
        if seed is not None:
            np.random.seed(seed = seed + i)
        x = df[i]
        if x[1] in Par.keys():
            Par[x[0]] = Par[x[1]]
        elif x[1][0].isdigit():
            Par[x[0]] = float(x[1])
        else:
            if x[1]=="unif":
                Par[x[0]] = U(float(x[2]), float(x[3]))
            elif x[1]=="log10unif":
                Par[x[0]] = LU(float(x[2]), float(x[3]))
            elif x[1]=="seq":
                Par[x[0]] = np.round(SEQ(float(x[2]), float(x[3]), float(x[4]), False), int(x[5]))
            elif x[1]=="log10seq":
                Par[x[0]] = np.round(SEQ(float(x[2]), float(x[3]), float(x[4]), True), int(x[5]))
            elif x[1]=="discrete":
                Par[x[0]] = np.random.choice(np.array(x[2:]).astype(float), 1)[0]
            elif x[1]=="math":
                x = [x[0],x[1]] + (x[2].split())
                ctd = " ".join([a.replace(a, "Par[\""+a+"\"]") if a in set(Par.keys()) else a for a in x[2:]])
                Par[x[0]] = eval(ctd)
            else:
                print(x)
                sys.exit("Error. Distribution not implemented.")
        np.random.seed(seed = None)
    Par.pop('seq_length', None)
    Par.pop('n_seq', None)
    Par.pop('mut_rate', None)
    Par.pop('recomb_rate', None)
    Par.pop('generation_time', None)
    return Par

def read_additional_param_file(location):
    Par = dict()
    with open(location, "r") as F:
        for line in F:
            if line.startswith("#") or line.strip()=="":
               continue
            line = line.strip().split()
            value = line[1:]
            if len(value) == 1:
                value = value[0]
            Par[line[0]] = value
    return Par

def validate(Par, Rules):
    isok = True
    for rule in Rules:
        isok = isok and eval(rule)
    return isok

def check_statistics(stat_file, conditions, min_SNP_D = 2000):
    STAT = []
    with open(stat_file, "r") as FS:
        for line in FS:
            if line.startswith("#"):
                continue
            line = line.strip().split(" ")
            STAT.append(line)
    sprime_stats = stat_file.replace(".stats", ".sprime2.stats")
    if os.path.exists(sprime_stats):
        SPRIME = []
        with open(sprime_stats, "r") as FS:
            for line in FS:
                if line.startswith("#"):
                    continue
                line = line.strip().split(" ")
                SPRIME.append(line)
    ###
    isok = True
    if conditions[0] is not None:
        C = np.array([x.split(":") for x in conditions])
        for i in range(len(C)):
            c = C[i]
            if len(c)!=3:
                sys.exit("Error. There should be 3 arguments in the condition "+c[0])
            if c[0]=="Fst":
                n, x = int(STAT[4][1]), float(STAT[4][2])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[Fst] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            elif c[0]=="D":
                n, x = int(STAT[3][1]), float(STAT[3][2])
                if n > min_SNP_D: # APPLY CONDITION ONLY IF THERE IS A SUFFICIENT NUMBER OF SNPs
                    so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                    isok = isok and so
                print("[D] ("+STAT[3][1]+") "+"%.3f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            elif c[0]=="AFS.CEU.1":
                n, x = int(STAT[0][1]), float(STAT[0][2])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[AFS.CEU.1] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            elif c[0]=="AFS.YRI.1":
                n, x = int(STAT[1][1]), float(STAT[1][2])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[AFS.YRI.1] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            elif c[0]=="pi.CEU":
                x = float(STAT[6][2])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[pi.CEU] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            elif c[0]=="pi.YRI":
                x = float(STAT[7][2])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[pi.YRI] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            ###
            elif c[0]=="DCFS.1":
                x = float(STAT[2][2])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[DCFS.1] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            ###
            elif c[0]=="SPRIME.match.mean":
                x = float(SPRIME[1][8])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[SPRIME.match.mean] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            elif c[0]=="SPRIME.length.mean":
                x = float(SPRIME[1][5])
                so = np.logical_and(x>=float(c[1]), x<=float(c[2]))
                isok = isok and so
                print("[SPRIME.length.mean] "+"%.2f"%(x)+str(np.where(so, "+", "-")), end = " | ")
            else:
                sys.exit("Error. Statistical condition is not implemented.")
    return isok

def lookup(cmd_stat, output, conditions, remove_dataset, write_geno1):
    cmd_stat += " --stats"
    if np.sum(np.array(["SPRIME" in cd for cd in conditions])) > 0:
        cmd_stat += " --sprime"
    subprocess.run(cmd_stat, shell=True, check=True, stdout=subprocess.DEVNULL)
    isok = check_statistics(output+".stats", conditions)
    if isok==False:
        os.remove(output+".geno.gz")
        if write_geno1:
            os.remove(output+".geno1.gz")
        os.remove(output+".snp.gz")
        os.remove(output+".ind")
        if remove_dataset:
            os.remove(output+".par")
            os.remove(output+".stats")
            if os.path.exists(output+".sprime2.stats"):
                os.remove(output+".sprime2.stats")
        else:
            os.rename(output+".par", output+".failed.par")
            os.rename(output+".stats", output+".failed.stats")
            if os.path.exists(output+".sprime2.stats"):
                os.rename(output+".sprime2.stats", output+".failed.sprime2.stats")
        print("(-) Failed criteria. Exiting.")
        sys.exit()
    else:
        os.remove(output+".stats")
        if os.path.exists(output+".sprime2.stats"):
            os.remove(output+".sprime2.stats")
        print("(+) Passed criteria. Pursuing.")

def check_instances(objs, instances):
    if not isinstance(objs, list):
        objs = list(objs)
    if not isinstance(instances, list):
        instances = list(instances)
    if len(objs) != len(instances):
        sys.exit("Error. objs must have same length as instances")
    for ix in range(len(objs)):
        if not isinstance(objs[ix], instances[ix]):
            sys.exit("Error. "+objs[ix]+" not of instance "+instances[ix])

def get_color(col, n, sub_n):
        cv = np.array(["#ffffff00"]*int(n))
        cv[np.linspace(0, n-1, num=sub_n, endpoint=True).astype(int)] = col
        cv = list(cv)
        return cv

# draw demography using DemesDraw module + export yaml + print_history into a file
def export_all(ARGS,
               PAR,
               Demography,
               Samples):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import demesdraw
    ######################################################################
    ############ Import plotting parameters
    Pars = pd.read_csv(ARGS["draw"], header=None, sep="\s+", skipinitialspace=True, na_filter=False, skip_blank_lines=True, comment="#").to_numpy()
    ParPlot = { x[0]: float(x[1]) for x in Pars if (x[0]!="is_qian") and (x[0]!="img.format") }
    ParPlot["is_qian"] = str(Pars[Pars[:,0]=="is_qian",1][0])
    if "img.format" in Pars[:,0]:
        ParPlot["img.format"] = Pars[Pars[:,0]=="img.format",1][0]
    else:
        ParPlot["img.format"] = "png"
    ######################################################################
    ############ Write history
    dp = Demography.debug()
    DEMA = open(ARGS["output"]+".demo", "w")
    dp.print_history(output = DEMA)
    DEMA.close()
    ######################################################################
    ############ Write location probabilities of sampled lineages
    if ParPlot["is_qian"] == "True":
        gt = float(ARGS["generation_time"])
        T = np.round(np.arange(ARGS["lin_prob_times"][0], ARGS["lin_prob_times"][1], step = ARGS["lin_prob_times"][2])) / gt
        X = dp.lineage_probabilities(T, sample_time = 0.)
        SampledDemes = [s.population for s in Samples]
        PD = {pop.name: pop.id for pop in Demography.populations}
        SampledDemes = [PD[s] for s in SampledDemes]
        with open(ARGS["output"]+".lprob", "w") as LP:
            LP.write("time.ga sampled.deme.index focal.deme proba.lineage.located.in.focal.deme\n")
            for t, x in zip(T, X):
                for s in SampledDemes:
                    xx = x[s,:].astype(str)
                    for ip, xxx in enumerate(xx):
                        LP.write(str(t)+" "+str(s)+" "+str(ip)+" "+xxx+"\n")
    ######################################################################
    ############ Plot the model
    #cols = ["#FF847C", "#99B898", "#E84A5F", "#2A363B"]
    plt.rcParams.update({'font.size': ParPlot["demo.font.size"]})
    if ParPlot["is_qian"] == "False":
        ### yaml>
        if type(Demography) == msp.demography.Demography:
            graph = Demography.to_demes()
        elif type(Demography) == demes.demes.Graph:
            graph = Demography
            Demography = msp.Demography.from_demes(graph)
        else:
            sys.exit("Error. Not implemented.")
        demes.dump(graph, ARGS["output"]+".yaml") # yaml-formatted demographic history
        # plot
        fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(ParPlot["demo.x.size"], ParPlot["demo.y.size"]), constrained_layout=False)
        plt.xticks(rotation=90)
        ax = demesdraw.tubes(graph,
                             ax = axs,
                             log_time = True)
        ax.figure.savefig(ARGS["output"]+".logT."+ParPlot["img.format"])
    elif ParPlot["is_qian"] == "True":
        ### yaml>
        graph = Demography.to_demes()
        demes.dump(graph, ARGS["output"]+".yaml") # yaml-formatted demographic history
        #
        fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(ParPlot["demo.x.size"], ParPlot["demo.y.size"]), constrained_layout=False)
        num_lines_per_migration = int(ParPlot["demo.n.migr.arrows"])
        tube_width = np.max( np.array([PAR[k] for k in PAR.keys() if k.startswith("n.0.")]).astype(float) ) * ParPlot["demo.tube.rescale"]
        nPops = len(Demography.populations)
        DemeSizes = [PAR[k] for k in PAR.keys() if k.startswith("C.")]
        nChains = len(DemeSizes)
        cols = sns.husl_palette(nChains, l = .4).as_hex()
        PopColors = []
        for i in range(nChains):
            PopColors += get_color(cols[i], int(DemeSizes[i]), nPops)
        PopColors = {key.name: PopColors[i] for i, key in enumerate(Demography.populations)}
        PopPositions = {key.name: tube_width*i for i, key in enumerate(Demography.populations)}
        ### ms> Not that the order is that given by the graph object!
        # Can only be done, currently (Nov 24th 2021) with a github dev version of demes
        """
        print("__/!\_ Beware! The order of populations and samples in the ms command depends on their order in the demes::graph object, which itself is determined by starting_time, so this can be completely reshuffled compared to the expected AfW-AfE-EuA-Nea order! Please check the *.msind file for an explicit listing of the sample order that is processed and output by ms.")
        Dm = np.array([dm["name"] for dm in graph.asdict_simplified()["demes"]])
        SampledMs = [0] * len(Demography.populations)
        for s in Samples:
            ids = np.where(Dm=="pop_"+str(s.population))[0]
            assert len(ids)==1, "No pop_*"
            SampledMs[ids[0]] = int(s.num_samples * int(ARGS["ploidy"]))
        with open(ARGS["output"]+".msind", "w") as MS:
            MS.write("#[Ind_ID] [Sex] [Pop_Index_In_Original_MSPRIME_Command_0based] [Pop_Index_In_MS_Command_1based]\n")
            iid = 1
            for pid, s in enumerate(SampledMs):
                for ni in range(0, s, ARGS["ploidy"]):
                    MS.write(str(iid)+" U "+str(Dm[pid]).replace("pop_", "")+" "+str(pid+1)+"\n")
                    iid += 1
        with open(ARGS["output"]+".ms", "w") as MS:
            N0 = 10000.
            mutrate = float(ARGS["mut_rate"])
            recrate = float(ARGS["rec_rate"])
            seqL = float(ARGS["genome"][1]) * 1e6
            MS.write(str(int(sum(SampledMs)))+" "+str(int(ARGS["genome"][0]))+" ")
            MS.write("-t "+str(4.*N0*mutrate*seqL)+" ")
            MS.write("-r "+str(4.*N0*recrate*(seqL-1.))+" "+str(int(seqL))+" ")
            MS.write(demes.to_ms(graph, # ms command
                        N0 = N0,
                        samples = SampledMs ))
            MS.write(" -p 9\n")
        """
        plt.xticks(rotation=90)
        ax = demesdraw.tubes(graph,
                             ax = axs,
                             colours = PopColors,
                             positions = PopPositions,
                             num_lines_per_migration = num_lines_per_migration,
                             log_time = True )
        ax.figure.savefig(ARGS["output"]+".logT."+ParPlot["img.format"])
    else:
        sys.exit("Error. Missing is_qian parameter in the draw par file.")

# By default, writes genetic positions in Morgans
def write_eigenstrat(treeseq,
                     prefix,
                     chromosome_id,
                     pop_labels,  # pop_labels should give the population labels of each haploid sample
                     diploidize,
                     assume_biallelic_only,
                     recombination_rate = 1e-8,
                     recombination_map = None,
                     append = False,
                     encode_ancestral_allele_as_0 = True,
                     verbose = False,
                     snp_write_gpos = False,
                     snp_write_alleles = False,
                     write_geno1 = False,
                     write_geno12 = False,
                     RunningTime = dict() ):
    ##### Checks
    assert not (diploidize==False and write_geno1==True), "Cannot diploidize while exporting a .geno1 file: contradictory."
    assert not (diploidize==False and write_geno12==True), "Cannot diploidize while exporting a .geno12 file: contradictory."
    if append:
        mode = "at"
    else:
        mode = "wt"
    leftChromosome = np.arange(0, len(pop_labels), int(2))
    assert len(pop_labels) == len(treeseq.samples()), "Length of pop_label should be equal to the number of columns in the genotype matrix."
    if verbose:
        if encode_ancestral_allele_as_0 == True:
            print("Ancestral allele is encoded as 0 i.e. 0=Anc/Anc, 2=Der/Der")
        else:
            print("Ancestral allele is encoded as 1 i.e. 0=Der/Der, 2=Anc/Anc")

    ##### File connections
    SNP = gzip.open(prefix+".snp.gz", mode)
    GENO = gzip.open(prefix+".geno.gz", mode)
    if diploidize:
        if write_geno1:
            GENO1 = gzip.open(prefix+".geno1.gz", mode)
        if write_geno12:
            GENO12 = gzip.open(prefix+".geno12.gz", mode)

    #### If genotypes are assumed to be biallelic only, then:
    if assume_biallelic_only:
        print("<fast write: assuming all SNPs strictly biallelic>")
        ################
        # Positions
        ################
        pos = treeseq.tables.sites.position.astype(np.int32) + 1 # tskit positioning is 0-based indexed
        nSNP, nSites = len(pos), len(pos)
        if snp_write_gpos:
            if recombination_map is None:
                gpos = pos * recombination_rate
            else:
                gpos = recombination_map.get_cumulative_mass(pos) # rates in RateMap are in Morgans already
            gpos = np.array(['%.8f' % x for x in gpos]) # genetic position in Morgans
        else:
            gpos = '.'
        ## alleles
        if snp_write_alleles:
            aa = 'A T'
        else:
            aa = ''
        snp = pd.DataFrame({ 'rs': '.',
                           'chrom': str(chromosome_id),
                           'gpos': gpos,
                           'pos': pos.astype(str),
                           'aa': aa })
        snp.to_csv( SNP,
                    sep = ' ',
                    header = False,
                    index = False)
        ################
        # Genotypes
        ################
        geno = treeseq.genotype_matrix()
        if len(np.unique(geno)) != 2:
            sys.exit("Error. Alleles must be strictly biallelic.")
        if encode_ancestral_allele_as_0 == False:
            geno = 1 - geno
        ## geno1
        if write_geno1:
            np.savetxt(GENO1, geno[:,leftChromosome], fmt='%s', delimiter='')
        ## geno12
        if write_geno12:
            np.savetxt(GENO12, geno, fmt='%s', delimiter='')
        ## geno
        if diploidize:
            np.savetxt(GENO, geno[:,leftChromosome] + geno[:,leftChromosome+1], fmt='%s', delimiter='')
        else:
            np.savetxt(GENO, geno, fmt='%s', delimiter='')

    #### Write only the biallelic variants
    else:
        nSNP, nSites = 0, 0
        for i, site in enumerate(treeseq.variants()):
            nSites += 1
            pos = int(site.site.position) + 1
            anc = site.site.ancestral_state
            gt = site.genotypes
            uq = set(gt)
            nuq = len(uq)
            anc_code = np.where(np.array(site.alleles)==anc)[0][0]
            if nuq==1:
                if gt[0]==anc_code:
                    pass
                else:
                    gt = gt*0+1
            elif nuq==2:
                if anc_code in uq:
                    gt = np.array([0 if g==anc_code else 1 for g in gt])
                else:
                    continue
            else:
                continue
            ###
            if encode_ancestral_allele_as_0 == False:
                gt = 1 - gt
            if diploidize:
                gt_diplo = gt[leftChromosome] + gt[leftChromosome+1]
                gt_diplo = gt_diplo.astype(str)
            gt = gt.astype(str)
            if diploidize and write_geno1:
                GENO1.write("".join(gt[leftChromosome])+"\n")
            if diploidize and write_geno12:
                GENO12.write("".join(gt)+"\n")
            if diploidize:
                GENO.write("".join(gt_diplo)+"\n")
            else:
                GENO.write("".join(gt)+"\n")
            # *.snp
            if snp_write_gpos:
                if recombination_map is None:
                    gpos = pos * recombination_rate
                else:
                    gpos = recombination_map.get_cumulative_mass(pos) # rates in RateMap are in Morgans already
                gpos = '%.8f'%gpos
            else:
                gpos = "."
            if snp_write_alleles:
                AA = anc+" . "
            else:
                AA = ""
            snp = ". "+str(chromosome_id)+" "+gpos+" "+str(pos)+" "+AA
            SNP.write(snp+"\n")
            nSNP += 1

    ##### Stats
    print("Nb. SNPs [raw, biallelic, %raw]:  "+str(nSites)+"  |  "+str(nSNP)+"  |  "+str( np.round(100*np.true_divide(nSNP, nSites),1) )+"%")

    ##### Close file connections
    SNP.close()
    GENO.close()
    if diploidize and write_geno1:
        GENO1.close()
    if diploidize and write_geno12:
        GENO12.close()

    ##### Write *.ind file
    with open(prefix+'.ind', 'w') as IND:
        if diploidize:
            for i in leftChromosome:
                IND.write(str(int(i/2+1))+' U '+str(pop_labels[i])+'\n')
        else:
            for i in range(len(pop_labels)):
                IND.write(str(i+1)+' U '+str(pop_labels[i])+'\n')

def msHOT_Map_SingleChrom(seed, SeqLen_kb, RecRate_OutHotspot, InterSpacing_kb, MinSpacing_kb, Width_kb, Log10Lambdas):
    rg = Generator(PCG64(seed))
    nHotspots = rg.poisson(lam = np.true_divide(float(SeqLen_kb), float(InterSpacing_kb)))
    #Lefts = np.sort(rg.choice(range(int(SeqLen_kb*1e3)), size=nHotspots, replace=False)).astype(np.int64)
    Lefts = np.sort(rg.uniform(0., float(SeqLen_kb*1e3), size=nHotspots)).astype(np.int64)
    Widths = (1e3 * rg.uniform(Width_kb[0], Width_kb[1], size = len(Lefts))).astype(np.int64)
    Lbds = 10 ** rg.uniform(Log10Lambdas[0], Log10Lambdas[1], size = len(Lefts))
    Position, Rate = [0], [RecRate_OutHotspot]
    for i in range(len(Lefts)):
        if (Position[-1] + MinSpacing_kb) >= Lefts[i]:
            continue
        Position += [Lefts[i]]
        Rate += [Lbds[i]*RecRate_OutHotspot]
        Position += [Lefts[i]+Widths[i]+1]
        Rate += [RecRate_OutHotspot]
    if Position[-1] >= int(SeqLen_kb*1e3) - 1:
        Position[-1] = int(SeqLen_kb*1e3) - 1
        Position += [int(SeqLen_kb*1e3)]
    else:
        Position += [int(SeqLen_kb*1e3)]
    rate_map = msp.RateMap(position = Position, rate = Rate)
    return rate_map

if __name__ == "__main__":
    main()
