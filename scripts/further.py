#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 23:44:56 2022
@author: rtournebize

Command:
    python3(.7) further.py parfile_to_simulate.par output_directory configuration_file.yaml

Will generate four directories:
    `data/` containing simulated dataset
    `stats_aDNA/` containing stats results of the aDNA analysis
    `stats_IBD/` containing stats results of the IBD analysis
    `stats_Modern/` containing stats results of the classical analysis on modern samples

300322  0.0
050422  0.1  now outputs the exception errors
150622  1.0  changes to accomodate the new simulate.py:::4.0.0 and sumstats.py:::14.2 scripts
211222  1.1  possibility to give a suffix for stats_* directories
050123  1.2  added option `--assume-no-na` to sumstats.py call to speed up exec
090123  1.3  conditional outdir (through midir) specification

"""

__version__ = "1.3"

import numpy as np
import subprocess, copy, os, yaml, sys, logging
import pandas as pd


#---------------- Functions ----------------#

def print_log(string):
    print(string)
    logging.warning(string)

def get_Args():
    parfile, outdir, config = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(config, "r") as F:
        Args = yaml.safe_load(F)
    # Derivation
    Args["parfile"] = parfile
    Args["outdir"] = outdir
    Args["run"] = int(Args["parfile"].split("/")[-1].replace(".par", ""))
    Args["seed"] = Args["run"]
    return Args

def mkdir_p(path):
    if os.path.exists(path) == False:
        os.makedirs(path)
    else:
        print(path+" already exists")

def get_Samples(Args):
    Times = np.array(Args["sampling_times"]).astype(float).astype(int)
    Locations = np.array(Args["sampling_locations"]).astype(float).astype(int)
    #
    M = pd.read_csv(Args["parfile"], sep = " ", header = None).to_numpy()
    T_OOA = float(M[M[:,0]=="T.ooa",1])
    D_OOA = float(M[M[:,0]=="j.ooa.EuA_into_AfE",1])
    Samples = ["AfW:pSample.YRI:50:0:YRI", "EuA:pSample.CEU:50:0:CEU", "Nea:pSample.Vindija:1:50000:Vindija", "Nea:pSample.Altai:1:130000:Altai"]
    # aDNA
    for location in Locations:
    	for time in Times:
            if (T_OOA - D_OOA * location) > time:
                Samples += [":".join(["EuA",str(location),"1",str(time),"aDNA.CEU."+str(location)+"."+str(time)])]
    # Present-day
    nDiploids = 15
    for location in Locations:
    	Samples += [":".join(["EuA",str(location),str(nDiploids),"0","Modern.CEU."+str(location)+".0"])]
    for location in Locations:
    	Samples += [":".join(["AfW",str(location),str(nDiploids),"0","Modern.YRI."+str(location)+".0"])]
    return Samples

#---------------- Prepare ----------------#

Args = get_Args()
Samples = get_Samples(Args)
Pops = [x.split(":")[4] for x in Samples]
mkdir_p(Args["outdir"])
logging.basicConfig(filename = Args["outdir"]+"/"+str(Args["run"])+".log",
                    filemode = "a",
                    format = '%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt = '%H:%M:%S',
                    level = logging.DEBUG)
logging.info("Start")


#---------------- Simulate ----------------#

if bool(Args["do_simul"]):
    logging.info("Simulate")
    print("Number of sampled populations:    "+str(len(Samples)))
    # command
    cmd = Args["python"]+" "+Args["script_qian"]
    cmd += " --model qian_demo"
    cmd += " --write_geno1"
    cmd += " --bin_dir "+Args["bin_dir"]
    cmd += " -i "+Args["parfile"]
    cmd += " -o "+Args["outdir"]+"/data/"+str(Args["run"])
    cmd += " -g "+str(Args["genome"][0])+" "+str(Args["genome"][1])
    cmd += " --mut_model "+Args["mut_model"]
    cmd += " --algo hudson"
    cmd += " --seed "+str(Args["seed"])
    cmd += " --samples_within "+" ".join(Samples)
    if bool(Args["do_hotspots"]):
        cmd += " --do_hotspots"
    # launch
    print(cmd)
    try:
        subprocess.run(cmd, shell = True, check = True)
    except Exception as e:
        print_log("Error while __simulating__:\n")
        print_log(e)
    print("-"*50)

#---------------- Sumstats ----------------#

if bool(Args["do_stats_aDNA"]) or bool(Args["do_stats_IBD"]) or bool(Args["do_stats_Modern"]):
    Inds = pd.read_csv(Args["outdir"]+"/data/"+str(Args["run"])+".ind", header = None, sep = " ").to_numpy()[:,2]
    RefStatsPar = dict()
    with open(Args["ref_stats_file"], "r") as F:
        for line in F:
            if line.strip() != "":
                line = line.split()
                RefStatsPar[line[0]] = " ".join(line[1:])

# aDNA
if bool(Args["do_stats_aDNA"]):
    logging.info("Stats ::: aDNA")
    midir = "stats_aDNA"
    if "outdir_stats_suffix" in Args.keys() and len(Args["outdir_stats_suffix"]) > 0:
        midir = midir+"."+Args["outdir_stats_suffix"]
    mkdir_p(Args["outdir"]+"/"+midir)
    cmd_ini = "python3 "+Args["script_sumstats"]+" --assume_no_na"
    cmd_ini += " --stats --SE --f4"
    if bool(Args["do_only_aDNA_D"]) == False:
        cmd_ini += " --ld1"
    cmd_ini += " --bin_dir "+Args["bin_dir"]
    cmd_ini += " -i "+(Args["outdir"]+"/data/"+str(Args["run"]))
    for pop in Pops:
        if pop.startswith("aDNA.CEU"):
            S = copy.deepcopy(RefStatsPar)
            stats_file = Args["outdir"]+"/"+midir+"/"+str(Args["run"])+"__"+pop+".par"
            S["AFS_num_inds"] = str(int(min([10, np.sum(Inds==pop)])))
            with open(stats_file, "w") as F:
                for key in S.keys():
                    F.write(key+"\t"+S[key]+"\n")
            # command
            cmd = cmd_ini
            cmd += " -p "+stats_file
            cmd += " -o "+(stats_file.replace(".par", ""))
            cmd += " --ceu "+pop
            cmd += " --yri "+"YRI"
            cmd += " --vindija "+"Vindija"
            cmd += " --altai "+"Altai"
            cmd += " --psmc_pops "+pop
            # launch
            try:
                subprocess.run(cmd, shell = True, check = True)
            except Exception as e:
                print_log("Error on __aDNA__ stats calculation for: "+stats_file+"\nCommand:\n"+cmd+"\nError:\n")
                print_log(e)
                print_log("\nPursuing.\n")
    print("-"*50)

# IBD
if bool(Args["do_stats_IBD"]):
    logging.info("Stats ::: IBD")
    midir = "stats_IBD"
    if "outdir_stats_suffix" in Args.keys() and len(Args["outdir_stats_suffix"]) > 0:
        midir = midir+"."+Args["outdir_stats_suffix"]
    mkdir_p(Args["outdir"]+"/"+midir)
    cmd_ini = "python3 "+Args["script_sumstats"]+" --assume_no_na"
    cmd_ini += " --stats --SE --sprime --crf --ld --f4"
    cmd_ini += " --bin_dir "+Args["bin_dir"]
    cmd_ini += " -i "+(Args["outdir"]+"/data/"+str(Args["run"]))
    for pop1 in Pops:
        for pop2 in Pops:
            if pop1.startswith("Modern.CEU") and pop2.startswith("Modern.YRI"):
                S = copy.deepcopy(RefStatsPar)
                stats_file = Args["outdir"]+"/"+midir+"/"+str(Args["run"])+"__"+pop1+"__"+pop2+".par"
                S["AFS_num_inds"] = str(int(min([10, np.sum(Inds==pop1), np.sum(Inds==pop2)])))
                with open(stats_file, "w") as F:
                    for key in S.keys():
                        F.write(key+"\t"+S[key]+"\n")
                # command
                cmd = cmd_ini
                cmd += " -p "+stats_file
                cmd += " -o "+(stats_file.replace(".par", ""))
                cmd += " --ceu "+pop1
                cmd += " --yri "+pop2
                cmd += " --vindija "+"Vindija"
                cmd += " --altai "+"Altai"
                cmd += " --psmc_pops "+pop1
                # launch
                try:
                    subprocess.run(cmd, shell = True, check = True)
                except Exception as e:
                    print_log("Error on __IBD__ stats calculation for: "+stats_file+"\nCommand:\n"+cmd+"\nError:\n")
                    print_log(e)
                    print_log("\nPursuing.\n")
    print("-"*50)


# Modern
if bool(Args["do_stats_Modern"]):
    logging.info("Stats ::: Modern")
    midir = "stats_Modern"
    if "outdir_stats_suffix" in Args.keys() and len(Args["outdir_stats_suffix"]) > 0:
        midir = midir+"."+Args["outdir_stats_suffix"]
    mkdir_p(Args["outdir"]+"/"+midir)
    #
    S = copy.deepcopy(RefStatsPar)
    stats_file = Args["outdir"]+"/"+midir+"/"+str(Args["run"])+".par"
    S["AFS_num_inds"] = "10"
    with open(stats_file, "w") as F:
        for key in S.keys():
            F.write(key+"\t"+S[key]+"\n")
    # command
    cmd = "python3 "+Args["script_sumstats"]+" --assume_no_na"
    cmd += " --stats --SE --sprime --crf --psmc --ld --f4"
    cmd += " --bin_dir "+Args["bin_dir"]
    cmd += " -p "+stats_file
    cmd += " -i "+(Args["outdir"]+"/data/"+str(Args["run"]))
    cmd += " -o "+(stats_file.replace(".par", ""))
    cmd += " --ceu CEU"
    cmd += " --yri YRI"
    cmd += " --vindija Vindija"
    cmd += " --altai Altai"
    cmd += " --psmc_pops CEU"
    # launch
    try:
        subprocess.run(cmd, shell = True, check = True)
    except Exception as e:
        print_log("Error on __Modern__ stats calculation for: "+stats_file+"\nCommand:\n"+cmd+"\nError:\n")
        print_log(e)
        print_log("\nPursuing.\n")
    print("-"*50)

logging.info("End")


#___