#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:24:33 2020
@author: rtournebize

2     011220   introduced PSMC + LD statistics; all statistics were checked against stats_v2.py results: OK
3     171220   introduced Sprime
4     211220   implemented the haploid version of CRF (as it should be)
5     140121   CRF: now outputs 1 for Neand geno if >=1 derived allele (example files say the contrary but the sup mat says so)
                 S' & CRF & LD : option for MAC-sensitive SNP-downsampling (to mimick 1000G properties), cf Sankararaman 14
                 S' : removed min_MatchRate and min_nSNP_per_segment filtering (they are not used in the paper in fact)
                 S' : outputs a 3rd line: the genome-wide match rate estimate (= single value for whole genome)
                 CRF : considers only the SNPs that are polymorphic in CEU (as done in Sankararaman 14)
                 CRF : introgression rate variation now computed across haploid individuals (instead of across chromosomes)
                 CRF : introduced option to filter out archaic segments < CRF_min_segment_length_cM
                 LD : now considers only the SNPs that are polymorphic in CEU (as done in Sankararaman 12)
6     190121   introduced the pi statistic (manually checked that pi=simulated_theta in a simple constant-No panmic model => OK)
7     220321   added uncertainty computation for the D statistic (added option --SE + option --block to specify block length)
                 added jackknife calculation for ancestry-LD (equivalent Sankararaman's `computedD`) + outputs *.cov.txt.gz only if --verbose
                 linked to remi_functions_v3 (instead of _v2) + outputs *.sprime.txt.gz only if --verbose
                 new ancestry-LD configuration was checked with simulations, cf. checks/LD => ok
                 added a sys.exit() in case of the presence of a 9 in the dataset -- we currently do not allow missing data
8     010421   added f3 calculation (outputs into *.stats)
                 changed Sprime records and outputs into a novel *.sprime2.stats output which contains different new statistics (eg. Detection Rate) and gives estimates
                 by varying the Minimum Match Rate for fragments to be considered introgressed (ie. to remove False Positives)
9     020521   changed the *.psmc.stats output to a single run: will systematically output the last iteration + few minor changes + added `--version`
10    190521   found that CRF could only work if the -o path was relative
                 CRF: reformatted the whole script (incl. debugged when no MAC downsampling (was reading single chromosome)) + no subset on *.crf.filtered.snps (already performed) + read gpos directly from *.SNP.snp
                 LD: added expfit with jacquelin inferred starting values + export all the jackknife runs
11    200521   stats: added a switch so that DCFS = NA if the number of target samples is greater or equal to `ndiploids_for_dcfs` (usually :=5)
                 LD: added the possibility to compute the single-sample LD statistic (Moorjani et al 2016) => introduced the --ld1 switch to call this mode
                 pseudodiploidization: functionality was moved out of the stats section (to be more general) => before, we could not do --ld without --stats
11.1  250521   bug in the condition leading to call pseudodiploidize(): `(do_stat)` instead of `(PAR["D_do_pseudodiploid"]=="YES" and do_stats)`
11.2  110621   LD: now outputs the number of SNP pairs and the NRMSD in the output (tried a convergence assessment approach, but not discriminative)
11.3  140621   can read only a subset of the SNPs using the `--chrom` option
                  LD: outputs on last line the sequence of age estimates for increasing number of chromosomes
11.4  180621   LD: added `--ld_Nea_ancestral` to ascertain SNPs where Nea has the ancestral allele (instead of derived, cf. Sankararaman 2012)
12.0  190621   major rewriting of the LD part: now --ld and --ld1 are in separate sections
               canceled the 11.4 change (ie removed the --ld_Nea_ancestral option)
               LD: now outputs 3 files corresponding to the three ascertainment schemes proposed in Sankararaman 2012
               LD: new parameter slot: "LD_sample_haploid_Neanderthal" in the stats parfile, which allows sampling a single chromosome from the geno1 file
                 |__. This is useful to exactly reproduce the Sankararaman's simulation procedure.
               stats: new parameter in par file "AFS_num_inds 10" (by default: 10) to restrict the number of inds on which to compute the AFSs
               general: added option `--bin_dir`, indicating the path to the dir containing the softwares (eg. sprime.jar) specified in the stat parfile
               LD & LD1: the block size is now defined as the seq length (instead of the number of SNPs polymorphic in CEU (v<12.0))
               checked by comparing results with sumstats 1.4: ok
               "psmc_samples_0based" parameter: changed the syntax: now write SampleIdx0Based:Label => this will write a file `output.Label.psmc2.stats`
12.1  240621   LD: bug in the Moorjani_0 ascertainment (missing variable name) + in --ld1, existed incremental_assessment but impossible for single sample LD
12.2  080921   Fst: following the Fst line in the *.stats output, now add one line: Fst calculated for SNPs only polymorphic in YRI (as in Bhatia 13)
                    added the option --Fst_ascertainment in relation to this (by default: False, meaning that the alternative ascertainment is not applied prior to Fst calculation)
13.0  120921   now able to handle chromosome of heterogeneous sizes -- in which case, the seq_length but be set to -1 in the *.par file
               the script will then use the last position of the chromosome as the chrom size -- note that improvement can be made setting chrom size as the last position in the genetic map (more accurate but in the end should change nearly nothing)
13.1  281121   now outputs also the median length for *.sprime2.stats
13.2  081221   was missing two nan to print per line in *.sprime2.stats, if no tracts detected
13.3  220322   now feeds read_fixedwidth with file handler instead of file name, because was throwing an error in latest pandas version (python3.10)
               generate_CRF_input: print the mean density of effective SNPs fed to CRF
14.0  270422   possibility now to specify the pops to analyze by argparse options
14.1  270422   included the --f4 option to calculate f4-ratios, this was checked in comparison with admixr results using simulations (also D-stats)
               checked 14.1 results vs. 13.3
14.2  020522   f4-ratio: moving window calculation with np.nansum instead of np.sum (in 14.1)
14.3  200622   now if `-v`, outputs the length of all S' segments
14.3.1  031122  microbug on get_sprime_stats() whereby script would through an error when there was no fragment at a particular min threshold (because forgot to return l, m in the exception)
15.0  111222   integrated the {extras.py} script v6.5 to avoid ext dependency
16.0  211222   added a new step to **subset sites** (either randomly to retain only `n` sites or using ascertainments);
               occurs _after_ the EIGENSTRAT data import step and _after_ the `na_rm` step (cf below)
17.0  291222   added option `--na_rm` which removes any SNP containing any `9` (missing genotype) and moved this step just before the SNP-subsetting step (`--subset`)
               added option `--assume_no_na` to speed up the calculation
               added a condition at the start of do_f4 to use the allele count matrices generated from do_stats, to speed up code
17.1  291222   removed some functions imported from previous extras.py script, which are useless here
18.0  291222   added option `--ancestral` to specify one single individual from the input *.ind file which specifies the ancestral genotypes
18.1  030122   debugged a small bug which occurred sometimes when generating PSMC output, when the very last bin was hetero (in prev version, would throw a stop error)
19.0  160123   introduced option `--ld_as_sanka12` which allows to calculate the ancestry-LD nearly exactly as in Sankararaman 2012; checked using simulations with computed vs. this script and agreement to 99%
               under `--ld_as_sanka12`, what changes:
                   - uses the biased covariance stat (divided by N) [instead of the unbiased (divided by N+1)]
                   - binned distances have a +1-incremented index [instead of non-incremented index]
                   - default value (for missing cov values) is 0 [instead of nan]
                   - do not set all cov values to nan if more than 50% of the bins have missing or infinite cov values [instead of setting them all to nan]
19.0.1  030523  added option to argparse to print out default values when using --help

______________________

_____OBSERVATIONS_____
______________________

240521. MAC-sensitive downsampling should not be performed for single-sample LD stat. Pseudodiploidization has very minor impact on LD decay rate (for both single and multi-sample stat). MAC-sensitive downsampling for multi-sample LD stat has little impact on decay rate estimates. Filtering out monomorphs has little impact on the decay rate estimates.

______________________

________PAR FILE______
______________________
Options:

- snp_subsetting ::: optional ::: if not added (or if set as `NO`), will never subset SNPs
                                    otherwise, if equal to a number, will subset SNPs to this number
                                    otherwise, if a string, will subset SNPs using an ascertainment scheme:
                                        * Archaic_array : at least one Neanderthal allele differs from
                                                          the majority allele in a panel of 24 YRI samples
                                                          (Fu et al 2015, Moorjani et al 2016)
                                 *NOTE* The subsetting is down *after* the whole data importation & optional sequence subsetting (`--chom`)
                                        & *before* MAC subsetting and pseudodiplo conversion

"""

___version___ = '19.0.1'

import numpy as np
import sys, re, argparse, time, allel, os, subprocess, gzip, copy
from os import path
import pandas as pd
from decimal import Decimal
from scipy.spatial.distance import pdist
from scipy.optimize import curve_fit

np.seterr(divide = 'ignore', invalid = 'ignore')

############################################################
############################################################
############################################################
############################################################
############################################################

def global_params():
    ndiploids_for_dcfs = 5
    # { minor_allele_count: proba_of_SNP_acceptance } from Sankararaman et al 2014 (SI2.1)
    acceptance_rates = {0: 0,
                    1: 0.25,
                    2: 0.5,
                    3: 0.75,
                    4: 0.8,
                    5: 0.9,
                    6: 0.95,
                    7: 0.96,
                    8: 0.97,
                    9: 0.98,
                    10: 0.99 }
    return ndiploids_for_dcfs, acceptance_rates

###############################################################################

def parse_options():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s: '+str(___version___))
    parser.add_argument('-i', '--input_prefix', type=str, required=True, help="Prefix of the input genetic data files")
    parser.add_argument('-o', '--output_prefix', type=str, required=False, default=None, help="Prefix of the output files")
    parser.add_argument('-p', '--parameter_file', type=str, required=True, help="Path to the statistical parameter file")
    parser.add_argument('--stats', action="store_true", default=False, help="Compute classical summary statistics?")
    parser.add_argument('--psmc', action="store_true", default=False, help="Compute PSMC curves?")
    parser.add_argument('--ld', action="store_true", default=False, help="Compute the ancestry LD on a multi-sample mode.")
    parser.add_argument('--ld1', action="store_true", default=False, help="Compute the ancestry LD on a single-sample mode.")
    parser.add_argument('--sprime', action="store_true", default=False, help="Compute S' analysis?")
    parser.add_argument('--crf', action="store_true", default=False, help="Compute CRF analysis?")
    parser.add_argument('--SE', action="store_true", default=False, help="Compute block jackknife standard errors?")
    parser.add_argument('--f4', action="store_true", default=False, help="Compute f4-ratios?")
    parser.add_argument('--block', type=str, default="5000kb", required=False, help="Block length for jackknife, either *kb for length in kb or *snp for length in a fixed number of SNPs, where * is a numeric value. A block length of 5 Mbp is used for D calculation in Fu et al. 2015, SI 16.")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, help="Add this option to output additional files.")
    parser.add_argument('--chrom', type=int, nargs="*", default=[None], required=False, help="List of chromosomes to analyze.")
    parser.add_argument('--bin_dir', type=str, required=True, help="Path to the directory containing the binaries/softwares.")
    parser.add_argument('--Fst_ascertainment', action="store_true", default=False, help="Add this switch to calculate Hudson Fst only with SNPs polymorphic in YRI.")
    # Samples
    parser.add_argument('--ceu', type=str, required=False, help="CEU population label, as in the *.ind file; alternatively, in the param file: CEU_indices_0based", default = None)
    parser.add_argument('--yri', type=str, required=False, help="YRI population label, as in the *.ind file; alternatively, in the param file: YRI_indices_0based", default = None)
    parser.add_argument('--vindija', type=str, required=False, help="Vindija_Neanderthal population label, as in the *.ind file; alternatively, in the param file: Vindija_indices_0based", default = None)
    parser.add_argument('--altai', type=str, required=False, help="Altai_Neanderthal population label, as in the *.ind file; alternatively, in the param file: Altai_indices_0based", default = None)
    parser.add_argument('--psmc_pops', type=str, required=False, nargs="*", help="Labels of the populations to analyze with PSMC, will select the first ind of each pop; alternatively, in the param file: psmc_samples_0based", default = None)
    parser.add_argument('--na_rm', action="store_true", default=False, help="Remove missing genotypes?") # v17.0
    parser.add_argument('--assume_no_na', action="store_true", default=False, help="Should the script assume there is no missing genotypes in the data?") # v17.1
    parser.add_argument('--ancestral', type=str, required=False, default = None, nargs=1, help="Population label, as in the *.ind file, of the single individual reference specifying the ancestral genotypes: note that we will polarize all alleles only if the ancestral genotypes are unambiguously 0 or 2") # v18.0
    parser.add_argument('--ld_as_sanka12', action='store_true', default=False, help="Add this switch to calculate the population-wise ancestry-LD with the same algorithm and parameters than Sankararaman et al. 12 original study, i.e. (1) using biased covariance stat (divide by N) instead of unbiased (divide by N+1); (2) binning distances with +1-incremented index instead of non-incremented index; (3) default value as 0 instead of nan; (4) do not set all cov values to nan if more than 50 percents of the bins have missing or infinite cov values") # v19.0
    return parser.parse_args()

def downsample(geno, snp, genetic_pos, acceptance_rates, generator, geno1 = None, columns = None):
    if columns is not None:
        nind = len(columns)
        MAC = np.sum(geno[:,columns], axis = 1).astype(int)
        MAC[MAC>nind] = 2*nind-MAC[MAC>nind]
    else:
        nind = geno.shape[1]
        MAC = np.sum(geno, axis = 1).astype(int)
        MAC[MAC>nind] = 2*nind-MAC[MAC>nind]
    MAC[MAC>10] = 10
    probs = np.array([acceptance_rates[x] for x in MAC])
    probs = generator.binomial(n = 1, p = probs)
    subset = np.array([False]*len(probs))
    subset[probs==1] = True
    print("Proportion of MAC-informed SNPs accepted:     "+str(np.round(np.true_divide(np.sum(subset), len(subset))*100,1))+"%")
    if geno1 is None:
        return geno[subset,:], snp[subset,:], genetic_pos[subset]
    else:
        return geno[subset,:], snp[subset,:], genetic_pos[subset], geno1[subset,:]

def count(geno, sam_labels, label, ncol = None):
    # matrix of 2 columns
    # first column:  count of the derived alleles
    # second column: count of the ancestral alleles
    cols = np.where(sam_labels==label)[0]
    if ncol is not None:
        cols = cols[0:int(ncol)]
    x = np.sum(geno[:,cols], axis = 1)
    x = np.column_stack((x, 2*len(cols)-x))
    return x

def union(list_of_doublets):
    b = []
    for begin,end in sorted(list_of_doublets):
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    return b

def P(x, n=6):
    f = "{0:."+str(int(n))+"f}"
    return f.format(x)

def populate(file, PAR):
    with open(file, "r") as FIN:
        for line in FIN:
            if line.startswith("#") or line=="\n":
                continue
            line = line.strip().replace("\t", " ")
            line = re.split("\s+", line)
            if ("_0based" in line[0]) and (line[0].startswith("psmc_")==False):
                if len(line)==2:
                    PAR[line[0]] = np.array([int(line[1])])
                else:
                    PAR[line[0]] = np.array([int(x) for x in line[1:]])
            elif ("_0based" in line[0]) and (line[0].startswith("psmc_")==True):
                if ":" not in line[1]:
                    sys.exit("New since v12: `psmc_samples_0based` must be of the form SampleIndex0Based:Label.")
                if len(line)==2:
                    PAR[line[0]] = np.array([line[1]])
                else:
                    PAR[line[0]] = np.array(line[1:])
            else:
                PAR[line[0]] = line[1]
    return PAR

def eigenstrat2psmcfa(geno, snp, index, seq_lengths, output_prefix, bin_size = 100):
    bin_size = float(bin_size)
    ss = snp[geno[:,index]==1,:]
    seqs = uniq(ss[:,0])

    count, FOUT = 0, open(output_prefix+".psmcfa", "w")
    for seq in seqs:
        n_bins = int(seq_lengths[seq] / bin_size) + (seq_lengths[seq] % bin_size != 0)
        count += 1
        sss = ss[ss[:,0]==seq,1]
        pos = np.true_divide(sss, bin_size)
        pos = np.floor(pos).astype(np.int64)
        pos[pos >= np.int64(n_bins)] = np.int64(n_bins - 1) # 18.1
        A = np.repeat("T", n_bins)
        A[pos] = "K"
        del sss, pos
        FOUT.write(">{}\n".format(count))
        for i in range(len(A)):
            if i>0 and i%60==0:
                FOUT.write("\n")
            FOUT.write(A[i])
        FOUT.write('\n')
    FOUT.close()

def gpos_check(PAR, genetic_pos, snp, output_unit = "cM", verbose = True):
    if PAR["genetic_pos_unit"].lower().startswith("m"):
        if output_unit=="cM":
            unit = Decimal("100.0")
        else:
            unit = Decimal("1.0")
        if verbose and genetic_pos is not None:
            print("Unit of genetic pos:   Morgans")
    elif PAR["genetic_pos_unit"].lower().startswith("c"):
        if output_unit=="cM":
            unit = Decimal("1.0")
        else:
            unit = Decimal("0.01")
        if verbose and genetic_pos is not None:
            print("Unit of genetic pos:   centiMorgans")
    else:
        sys.exit("Unit of genetic positions not recognized.")

    # if no genetic_position provided in the snp file
    if genetic_pos is None:
        if verbose:
            print("No genetic positions provided: inferring using the uniform recombination rate r="+str(PAR["recomb_rate"])+".")
        rec = Decimal(str(PAR["recomb_rate"]))
        if output_unit=="cM":
            genetic_pos = np.array([Decimal("100")*(x*rec) for x in snp[:,1]])
        else:
            genetic_pos = np.array([x*rec for x in snp[:,1]])

    # if provided
    else:
        genetic_pos = unit*genetic_pos

    return genetic_pos

def generate_Sprime_input(geno, snp, genetic_pos, ind, samples2, output):
    # replacement dictionary
    GT = {2: "1/1", 1: "0/1", 0: "0/0", 9: "./."}
    G_NEA = dict()

    ceu, yri = np.where(samples2==0)[0], np.where(samples2==1)[0]
    ind2 = "\t".join(ind[np.concatenate((ceu, yri)),0])
    gg = geno[:,np.concatenate((ceu, yri))]
    gNEA = geno[:,samples2==2]

    # keep only the polymorphs in (CEU, YRI)
    sub = np.sum(gg, axis = 1)
    sub = np.logical_and(sub > 0, sub < 2*gg.shape[1])
    gg, pp, gp, gNEA = gg[sub,:], snp[sub,:], genetic_pos[sub], gNEA[sub,0]

    # remove duplicate positions
    x = np.array([x1+":"+x2 for x1, x2 in zip(pp[:,0].astype(str), pp[:,1].astype(str))])
    unique, sub = np.unique(x, return_index=True); del unique
    sub = np.sort(sub)
    gg, pp, gp, gNEA = gg[sub,:], pp[sub,:], gp[sub], gNEA[sub]

    # export
    with open(output+".map", "w") as MAP:
        with open(output+".vcf", "w") as VCF:
            VCF.write("##fileformat=VCFv4.0\n")
            VCF.write("##source=rt.py\n")
            VCF.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            VCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+ind2+"\n")
            for i in range(gg.shape[0]):
                VCF.write(str(pp[i,0])+"\t"+str(pp[i,1])+"\t"+str(pp[i,0])+":"+str(pp[i,1])+"\tA\tT\t100\t.\t.\tGT\t")
                VCF.write("\t".join([GT[x] for x in gg[i,:]]))
                VCF.write("\n")
                MAP.write(str(pp[i,0])+" "+str(pp[i,0])+"."+str(pp[i,1])+" "+"%.8f"%(gp[i])+" "+str(int(pp[i,1]))+"\n")
                # populate the dictionary
                G_NEA[ str(pp[i,0])+":"+str(pp[i,1]) ] = gNEA[i]

    with open(output+".YRI.ind", "w") as fIND:
        for i in range(len(yri)):
            fIND.write(ind[yri[i],0]+"\n")

    return G_NEA

def rm_if_exists(file):
    if path.exists(file):
        os.remove(file)

def summary(x):
    a = ["n", "Min", "Q2.5p", "Q25p", "Median", "Mean", "Std", "Q75p", "Q97.5p", "Max"]
    if len(x)==0:
        b = [0]+["None"]*(len(a)-1)
    else:
        b = [len(x), np.min(x), np.quantile(x, 0.025), np.quantile(x, 0.25), np.median(x), np.mean(x),
             np.sqrt(np.var(x)), np.quantile(x, 0.75), np.quantile(x, 0.975), np.max(x)]
    return a, b

def NRMSD(y, yfit):
    rmsd = np.sqrt(np.mean((np.array(y) - np.array(yfit))**2))
    nrmsd = rmsd * (np.max(yfit) - np.min(yfit))**-1
    return nrmsd

def interweave(a, b):
    c = np.empty((a.size + b.size,), dtype=a.dtype)
    c[0::2] = a
    c[1::2] = b
    return c

def get_haplotypes(geno, geno1):
    geno2 = geno - geno1
    return interweave(geno1, geno2)

def generate_CRF_input(PAR, geno, geno1, snp, genetic_pos, ind, samples2, seq, output, bin_dir, seq_lengths):
    # currently prints all YRI and NEA samples in the temporary *.geno
    sub = np.where(snp[:,0]==seq)[0]
    gg, pp, gp, gg1 = geno[sub,:], snp[sub,:], genetic_pos[sub], geno1[sub,:]
    seq_len = seq_lengths[seq]
    assert seq_len >= pp[-1,1], "Error in generate_CRF_input: seq_lengths[seq] < last SNP position in bp"
    # remove duplicate positions
    unique, sub = np.unique(pp[:,1], return_index=True); del unique
    sub = np.sort(sub)
    gg, pp, gp, gg1 = gg[sub,:], pp[sub,:], gp[sub], gg1[sub,:]
    with open(output+".CEU.ind", "w") as fIND:
        ii = ind[samples2==0,:]
        for i in range(ii.shape[0]):
            fIND.write(ii[i,0]+" U "+ii[i,1]+"\n" + ii[i,0]+".2 U "+ii[i,1]+"\n")
    fCEU, fYRI, fNEA = open(output+".CEU.geno", "w"), open(output+".YRI.geno", "w"), open(output+".NEA.geno", "w")
    fqYRI, fqNEA, fSNP = open(output+".YRI.freq", "w"), open(output+".NEA.freq", "w"), open(output+".SNP.snp", "w")
    nCEU = np.sum(samples2==0)
    nSNP = 0
    for i in range(gg.shape[0]):
        # output only the SNPs that are polymorphic in CEU
        dac = np.sum(gg[i,samples2==0])
        if dac==0 or dac==2*nCEU:
            continue
        # SNP file
        fSNP.write(str(pp[i,0])+":"+str(pp[i,1])+" "+str(pp[i,0])+" "+str(gp[i]*Decimal('0.01'))+" "+str(pp[i,1])+" A T\n")
        # CEU
        x = get_haplotypes(gg[i,samples2==0], gg1[i,samples2==0])
        fCEU.write("".join(x.astype(str))+"\n")
        # YRI
        x = get_haplotypes(gg[i,samples2==1], gg1[i,samples2==1])
        fYRI.write("".join(x.astype(str))+"\n")
        fqYRI.write(str(pp[i,0])+":"+str(pp[i,1])+" "+str(np.true_divide(np.sum(x), len(x)))+"\n")
        # NEA
        x = get_haplotypes(gg[i,samples2==2], gg1[i,samples2==2])
        fNEA.write("".join(x.astype(str))+"\n")
        #fNEA.write(str(int(np.sum(x)>0))+"\n")
        fqNEA.write(str(pp[i,0])+":"+str(pp[i,1])+" "+str(np.true_divide(np.sum(x), len(x)))+"\n")
        nSNP += 1
    fCEU.close(); fYRI.close(); fNEA.close(); fqNEA.close(); fqYRI.close(); fSNP.close()
    print("Chrom "+str(seq)+": "+str(nSNP)+" SNPs => "+str(np.round(np.true_divide(nSNP, seq_len)*1e3, 1))+" SNPs/kb")
    # parameter file
    with open(bin_dir+"/"+PAR["CRF_parFile"], "r") as fCRF:
        with open(output+".CRF.par", "w") as oCRF:
            oCRF.write("PREFIX:            "+output+"\n")
            oCRF.write("chr:               "+str(seq)+"\n")
            oCRF.write("parameters:        "+str(bin_dir+"/"+PAR["CRF_parEst"])+"\n")
            for line in fCRF:
                oCRF.write(line)

def zero_runs(a):
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def get_pos(string):
    return int(string.split(":")[1])

def get_gpos(snp, gpos, seq, phys_pos):
    z = np.where(np.logical_and(snp[:,0]==seq, snp[:,1]==phys_pos))[0]
    if len(z)==0:
        sys.exit("Error multiple gpos for one ppos")
    elif len(z)==1:
        z = gpos[z]
    else:
        z = gpos[z[0]]
    return z

def pseudodiploidize(M, generator):
    geno = copy.deepcopy(M)
    rep = geno==1
    geno[rep] = generator.choice([0,2], size=np.sum(rep), replace=True).astype(np.int8)
    return geno

def write_psmc_stats(file, output, run, mut_rate, bin_size):
    read = False
    X, Y = [], []
    with open(file, "r") as F:
        for line in F:
            if read and line.startswith("//"):
                break
            if line.startswith("RD\t"+str(int(run))):
                read = True
            if read:
                if line.startswith("RS\t"):
                    linarr = line.split("\t")
                    X.append(float(linarr[2]))
                    Y.append(float(linarr[3]))
                if line.startswith("PA\t"):
                    theta = float(line.split("\t")[1].split(" ")[1])
    X, Y = np.array(X), np.array(Y)
    No = np.true_divide(theta, 4. * float(mut_rate) * float(bin_size))
    X = 2. * X * No
    Y *= No
    with open(output, "w") as FOUT:
        FOUT.write("# mut_rate: "+str(mut_rate)+" run: "+str(run)+" theta: "+str(theta)+"\n")
        FOUT.write("# gBP Ne\n")
        for i in range(len(X)):
            FOUT.write(str(int(round(X[i])))+"\t"+str(int(round(Y[i])))+"\n")

def get_sprime_stats(x, min_mr, genome_size, min_nsnp = 3):
    x = x[x[:,0]>=float(min_nsnp),:]
    x = x[x[:,2]>=float(min_mr),:]
    if x.shape[0]==0:
        return [min_mr, 0]+[np.nan]*11, np.nan, np.nan
    l, m = x[:,1], x[:,2]
    dr = np.true_divide(np.sum(l), genome_size)
    # non-overlapping segments
    s = 0
    for q in np.unique(x[:,5]):
        z = x[x[:,5]==q,:]
        z = union([(z[i,3], z[i,4]) for i in range(z.shape[0])])
        s += np.sum(np.array([zz[1]-zz[0]+1 for zz in z]))
    dr_noov = np.true_divide(s, genome_size)
    # end
    gm = np.true_divide(np.sum(x[:,0]*m), np.sum(x[:,0]))
    return [min_mr, x.shape[0], dr, dr_noov, gm,
     np.mean(l), np.quantile(l, 0.025), np.quantile(l, 0.975),
     np.mean(m), np.quantile(m, 0.025), np.quantile(m, 0.975),
     np.median(l), np.median(m)], l, m

def main():
    ndiploids_for_dcfs, acceptance_rates = global_params()

    start_time = time.time()

    options = parse_options()
    prefix = options.input_prefix
    output = options.output_prefix
    parfile = options.parameter_file
    verbeux = options.verbose
    Chrom = options.chrom
    # statistic switches
    do_stats = options.stats
    do_psmc = options.psmc
    do_ld = options.ld
    do_ld1 = options.ld1
    if do_ld1 and do_ld:
        sys.exit("You can provide only one LD mode (either --ld or --ld1)")
    do_sprime = options.sprime
    do_crf = options.crf
    do_SE = options.SE
    do_f4 = options.f4
    bin_dir = options.bin_dir
    # block stats
    block_length = options.block
    # v12.2
    do_Fst_ascertainment = options.Fst_ascertainment
    # v17.0
    na_rm = options.na_rm
    assume_no_na = options.assume_no_na
    # v18.0
    ancestral = options.ancestral
    if assume_no_na and na_rm:
        sys.exit("You cannot specify altogether --na_rm and --assume_no_na")

    print("v"+str(___version___))
    print("allel: "+str(allel.__version__)+"\n")

    if output is None:
        output = prefix

    if not path.exists(prefix+".geno") and not path.exists(prefix+".geno.gz"):
        sys.exit(prefix+".geno(.gz) does not exist")
    if not path.exists(prefix+".snp") and not path.exists(prefix+".snp.gz"):
        sys.exit(prefix+".snp(.gz) does not exist")
    if not path.exists(prefix+".ind") and not path.exists(prefix+".ind.gz"):
        sys.exit(prefix+".ind(.gz) does not exist")
    if path.exists(prefix+".geno") and path.exists(prefix+".geno.gz"):
        sys.exit(prefix+".geno and .geno.gz both exist")
    if path.exists(prefix+".snp") and path.exists(prefix+".snp.gz"):
        sys.exit(prefix+".snp and .snp.gz both exist")
    if not path.exists(parfile):
        sys.exit(parfile+" does not exist")

    # read parameters
    PAR = {"CEU_indices_0based": None, "YRI_indices_0based": None, "Vindija_index_0based": None, "Altai_index_0based": None, "seed": None, "psmc_samples_0based": None}
    PAR = populate(parfile, PAR)
    PAR = populate(prefix+".par", PAR)
    PAR["ld_as_sanka12"] = options.ld_as_sanka12

    # checks (v14.0)
    IndFile = pd.read_csv(prefix+".ind", sep = "\s+", usecols = [2], engine = "c", dtype = str, na_filter = None, header = None).to_numpy()[:,0]
    ## all
    optSam = {"CEU": options.ceu, "YRI": options.yri, "Vindija": options.vindija, "Altai": options.altai}
    for key in optSam.keys():
        if key=="Vindija" or key=="Altai":
            sp = "index"
        else:
            sp = "indices"
        if PAR[key+"_"+sp+"_0based"] is None and optSam[key]==None and key!="Altai":
            sys.exit("Missing specification of "+key+" samples")
        elif PAR[key+"_"+sp+"_0based"] is not None and optSam[key]!=None:
            sys.exit("You cannot specify altogether "+key+"_"+sp+"_0based in the par file and in the command line.")
        elif PAR[key+"_"+sp+"_0based"] is None and optSam[key]!=None:
            PAR[key+"_"+sp+"_0based"] = np.where(IndFile==optSam[key])[0]
            if len(PAR[key+"_"+sp+"_0based"]) == 0:
                sys.exit("We found 0 samples for "+key)
    ## f4
    if do_f4: assert PAR["Altai_index_0based"] is not None and PAR["Vindija_index_0based"] is not None, "If --f4, you must specify two Neanderthal samples."
    if do_f4: assert do_stats, "If --f4, you must also add the switch --stats."
    ## psmc
    if options.psmc_pops is None and PAR["psmc_samples_0based"] is None and do_psmc:
        sys.exit("Missing specification of the PSMC samples")
    elif options.psmc_pops is not None and PAR["psmc_samples_0based"] is not None:
        sys.exit("You cannot specify altogether psmc_samples_0based in the par file and in the command line.")
    elif options.psmc_pops is not None and PAR["psmc_samples_0based"] is None:
        PAR["psmc_samples_0based"] = []
        for p in options.psmc_pops:
            print(p)
            PAR["psmc_samples_0based"].append(str(np.where(IndFile==p)[0][0])+":"+p)

    if PAR["seed"]=="None":
        seed = None
    else:
        seed = int(PAR["seed"])
    rng = np.random.default_rng(seed = seed)

    if not "LD_sample_haploid_Neanderthal" in PAR.keys():
        PAR["LD_sample_haploid_Neanderthal"] = "NO"

    nCEU, nYRI, nVindija = len(PAR["CEU_indices_0based"]), len(PAR["YRI_indices_0based"]), len(PAR["Vindija_index_0based"])
    if PAR["Altai_index_0based"] is None:
        nAltai = 0
    else:
        nAltai = len(PAR["Altai_index_0based"])
    samples1 = np.concatenate((PAR["CEU_indices_0based"], PAR["YRI_indices_0based"], PAR["Vindija_index_0based"]))
    if PAR["Altai_index_0based"] is not None:
        samples1 = np.concatenate((samples1, PAR["Altai_index_0based"]))
    samples2 = np.array([0]*nCEU + [1]*nYRI + [2]*nVindija + [3]*nAltai)
    # CEU: 0
    # YRI: 1
    # Vindija: 2
    # Altai: 3

    if ancestral is not None:
        ancestral_idx = np.where(IndFile==ancestral)[0]
        if len(ancestral_idx) != 1:
            sys.exit("Error. There must be only a single ancestral individual.")
        ancestral_idx = ancestral_idx[0]
        samples1 = np.concatenate((samples1, np.array([ancestral_idx])))
    else:
        ancestral_idx = None

    print("\n____SUMMARY STATISTICS___\n")
    print("Parameter file:        "+parfile)
    print("Input prefix:          "+prefix)
    print("Output prefix:         "+output)
    if Chrom[0] is not None:
        print("Sequences to analyze:  "+" ".join([str(x) for x in Chrom]))
    print("YRI idx:               "+"(n="+str(len(PAR["YRI_indices_0based"]))+") "+" ".join([str(x) for x in PAR["YRI_indices_0based"]]))
    print("CEU idx:               "+"(n="+str(len(PAR["CEU_indices_0based"]))+") "+" ".join([str(x) for x in PAR["CEU_indices_0based"]]))
    print("Vindija idx:           "+"(n="+str(len(PAR["Vindija_index_0based"]))+") "+" ".join([str(x) for x in PAR["Vindija_index_0based"]]))
    if nAltai > 0:
        print("Altai idx:             "+"(n="+str(len(PAR["Altai_index_0based"]))+") "+" ".join([str(x) for x in PAR["Altai_index_0based"]]))
    print("Ancestral idx:         "+str(ancestral_idx))
    print("Pseudodiploid for D:   "+str(PAR["D_do_pseudodiploid"]))
    print("Seed:                  "+str(seed))
    print("Calculating SEs:       "+str(do_SE))
    print("Block length:          "+block_length)
    print("Software directory:    "+bin_dir)
    if path.exists(bin_dir+"/"+PAR["Sprime"])==False or path.exists(bin_dir+"/"+PAR["CRF_script"])==False or path.exists(bin_dir+"/"+PAR["CRF_parFile"])==False:
        sys.exit(bin_dir+"/"+PAR["Sprime"]+" or "+bin_dir+"/* do(es) not exist.")
    print("Ascertained Fst:       "+str(do_Fst_ascertainment))
    print("Assume no NA:          "+str(assume_no_na))
    print("Remove SNPs with NA:   "+str(na_rm))
    print("")

    ind, snp, genetic_pos, geno = read_eigenstrat(prefix,
                             read_geno = True,
                             read_snp = True,
                             read_ind = True,
                             ancestral_allele_encoded_as_0 = True,
                             columns = samples1)

    # subset, if requested
    if Chrom[0] is not None:
        chrom_ok = np.in1d(snp[:,0], np.array(Chrom), assume_unique = False)
        snp = snp[chrom_ok,:]
        if genetic_pos is not None:
            genetic_pos = genetic_pos[chrom_ok]
        geno = geno[chrom_ok,:]

    seq_names = np.unique(snp[:,0], return_index = True)
    seq_end = snp[np.concatenate((seq_names[1][1:]-1, np.array([-1]))),:][:,1]
    seq_names = seq_names[0]
    n_seq = len(seq_names)

    # estimate sequence lengths, if heterogeneous
    if float(PAR["seq_length"]) == -1.:
        seq_lengths = { seq_names[sn]: float(seq_end[sn]) for sn in range(n_seq) }
    else:
        seq_lengths = { seq_names[sn]: float(PAR["seq_length"]) for sn in range(n_seq) }
    print("Sequence lengths:")
    for sn in seq_names:
        print(str(sn)+": "+str(np.round(seq_lengths[sn]*1e-6, 2))+" Mbp")

    # estimate genome size
    genome_size = float(sum(list(seq_lengths.values())))
    print("=========================")
    print(str(np.round(genome_size*1e-6, 2)),"Mbp ("+str(int(genome_size))+" bp)\n")

    if (do_ld==True) or (do_sprime==True) or (do_crf==True) or (do_ld1==True):
        genetic_pos = gpos_check(PAR, genetic_pos, snp, output_unit = "cM", verbose = True)

    if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
        print("Reading *.geno1.")
        if not path.exists(prefix+".geno1") and not path.exists(prefix+".geno1.gz"):
            sys.exit(prefix+".geno1(.gz) does not exist")
        if path.exists(prefix+".geno1") and path.exists(prefix+".geno1.gz"):
            sys.exit(prefix+".geno1 and .geno1.gz both exist")
        if path.exists(prefix+".geno1.gz"):
            FGENO = gzip.open(prefix+".geno1.gz", "rt")
        else:
            FGENO = open(prefix+".geno1", "r")
        geno1 = read_fixedwidth(FGENO, dtype = np.int8, columns = samples1)
        FGENO.close()
        if Chrom[0] is not None:
            geno1 = geno1[chrom_ok,:]
            del chrom_ok
        if geno.shape != geno1.shape:
            sys.exit("geno and geno1 must have the same dimensions.")

    print("Reading done.\n")

    # polarize alleles
    if ancestral is not None:
        print("Polarizing alleles according to ancestral reference:")
        anc_gt = geno[:,-1]
        anc_gt_0 = anc_gt==0
        anc_gt_2 = anc_gt==2
        #
        geno[anc_gt_2,0:-1] = 2-geno[anc_gt_2,0:-1]
        geno[geno<0] = 9
        ok = np.logical_or(anc_gt_0, anc_gt_2)
        nori = geno.shape[0]
        geno = geno[ok,0:-1]
        if genetic_pos is not None:
            genetic_pos = genetic_pos[ok]
        snp = snp[ok,:]
        if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
            geno1 = geno1[ok,0:-1]
            geno1[geno1<0] = 9
        del ok, anc_gt, anc_gt_0, anc_gt_2
        print("   Removed:   "+str(nori-geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(nori-geno.shape[0],nori),1))+"%)")
        print("   Polarized: "+str(geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(geno.shape[0],nori),1))+"%)\n")
        #
        samples1 = samples1[0:-1]

    # handle missing data in the *.geno file
    if assume_no_na:
        print("_/!\__ Assumes that there are NO missing genotypes in the input *geno files.")
    else:
        # remove all SNPs containing at least one missing genotype (9)
        bads = np.any(geno==9, axis = 1)
        if np.any(bads==True):
            if na_rm:
                nori = geno.shape[0]
                ok = (bads==False)
                geno = geno[ok,:]
                if genetic_pos is not None:
                    genetic_pos = genetic_pos[ok]
                snp = snp[ok,:]
                if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
                    geno1 = geno1[ok,:]
                del ok, bads
                print("Removing SNPs having missing genotypes:")
                print("   Removed:   "+str(nori-geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(nori-geno.shape[0],nori),1))+"%)")
                print("   Retained:  "+str(geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(geno.shape[0],nori),1))+"%)\n")
            else:
                sys.exit("Error. Presence of missing data. We currently do not handle missing genotypes.\n")
        else:
            pass

    # if requested, subset sites
    if "snp_subsetting" in PAR.keys():
        if PAR["snp_subsetting"].upper() == "NO":
            print("Not subsetting sites.\n")
            pass
        else:
            if PAR["snp_subsetting"].isdigit():
                if int(PAR["snp_subsetting"]) >= snp.shape[0]:
                    subSites = range(snp.shape[0])
                else:
                    subSites = np.sort(rng.choice(range(snp.shape[0]), int(PAR["snp_subsetting"])))
            elif PAR["snp_subsetting"] == "Archaic_array":
                if not 3 in samples2:
                    sys.exit("You must sample an Altai Neand __and__ a Vindija Neand for `Archaic` snp_subsetting")
                subSites = ascertain(geno,
                            method = PAR["snp_subsetting"],
                            target_idx = np.where(samples2==0)[0],
                            outgroup_idx = np.where(samples2==1)[0],
                            archaic_idx = np.concatenate((np.where(samples2==2)[0], np.where(samples2==3)[0])))
            else:
                sys.exit("snp_subsetting not implemented.")
            #
            nori = geno.shape[0]
            geno = geno[subSites,:]
            if genetic_pos is not None:
                genetic_pos = genetic_pos[subSites]
            snp = snp[subSites,:]
            if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
                geno1 = geno1[subSites,:]
            del subSites
            print("Subsetting sites based on: `snp_subsetting="+PAR["snp_subsetting"]+"`")
            print("   Removed:   "+str(nori-geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(nori-geno.shape[0],nori),1))+"%)")
            print("   Retained:  "+str(geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(geno.shape[0],nori),1))+"%)\n")
    else:
        print("Not subsetting sites.\n")

    # MAC-sensitive SNP-downsampling
    # To mimick 1000G, the MAC should be computed across sapiens populations only (excluding Neanderthal)
    if do_sprime or do_crf or do_ld or do_ld1:
        if PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES":
            columns = np.concatenate((np.where(samples2==0)[0], np.where(samples2==1)[0]))
            if do_crf or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
                xs_geno, xs_snp, xs_genetic_pos, xs_geno1 = downsample(geno, snp, genetic_pos, acceptance_rates, rng, geno1 = geno1, columns = columns)
            else:
                xs_geno, xs_snp, xs_genetic_pos = downsample(geno, snp, genetic_pos, acceptance_rates, rng, geno1 = None, columns = columns)

    # Pseudodiploidize
    if (do_stats) or (PAR["LD_do_pseudodiploid"]=="YES" and (do_ld or do_ld1)):
        geno_pd = pseudodiploidize(geno, rng)

    ################################################################################################
    if do_stats:

        if "AFS_num_inds" not in PAR.keys():
            PAR["AFS_num_inds"] = 10
        else:
            PAR["AFS_num_inds"] = int(PAR["AFS_num_inds"])
        print("Nb diploids for AFS:   "+str(PAR["AFS_num_inds"]))

        # count matrices have 2 columns: the first gives the count of the derived allele, the second the count of the ancestral allele
        ac_ceu = count(geno, samples2, 0)
        ac_yri = count(geno, samples2, 1)
        ac_nea = count(geno, samples2, 2)
        ac_ref = ac_nea*0
        ac_ref[:,1] = 2

        FOUT = open(output+".stats", "w")
        FOUT.write("#sumstats: "+str(___version___)+"\n")

        ################################################################################################
        # 2D-AFS & 1D-AFS

        if nCEU > PAR["AFS_num_inds"]:
            ac_ceu_2 = count(geno, samples2, 0, ncol = PAR["AFS_num_inds"])
            n1 = PAR["AFS_num_inds"]
        elif nCEU < PAR["AFS_num_inds"]:
            sys.exit("Error. AFS_num_inds must be <= nCEU")
        else:
            ac_ceu_2 = ac_ceu
            n1 = nCEU

        if nYRI > PAR["AFS_num_inds"]:
            ac_yri_2 = count(geno, samples2, 1, ncol = PAR["AFS_num_inds"])
            n2 = PAR["AFS_num_inds"]
        elif nCEU < PAR["AFS_num_inds"]:
            sys.exit("Error. AFS_num_inds must be <= nYRI")
        else:
            ac_yri_2 = ac_yri
            n2 = nYRI

        afs = allel.joint_sfs(ac_ceu_2[:,0], ac_yri_2[:,0], n1=2*n1, n2=2*n2)
        afs_ceu = np.sum(afs, axis = 1)
        afs_ceu = afs_ceu[1:-1]
        FOUT.write("AFS_CEU "+str(np.sum(afs_ceu))+" ")
        FOUT.write(" ".join([P(x) for x in np.true_divide(afs_ceu, np.sum(afs_ceu))])+"\n")
        afs_yri = np.sum(afs, axis = 0)
        afs_yri = afs_yri[1:-1]
        FOUT.write("AFS_YRI "+str(np.sum(afs_yri))+" ")
        FOUT.write(" ".join([P(x) for x in np.true_divide(afs_yri, np.sum(afs_yri))])+"\n")

        ################################################################################################
        # DCFS

        ceu_idx = np.where(samples2==0)[0]
        if len(ceu_idx) >= ndiploids_for_dcfs:
            """
            YRI and NEA are pseudodiploidized following the Yang et al (2012) procedure:
                "In order to build the dcfs, at each site, we sampled one read
                at random for the Neanderthal and called it ancestral if it
                matched the reconstructed ancestor or derived if it did not.
                Similarly for the Yoruba, at each position, we sampled one
                chromosome at random and compared it to the recon-
                structed ancestor."
            """
            g_yris = geno_pd[:,np.where(samples2==1)[0]]
            g_nea = geno_pd[:,np.where(samples2==2)[0][0]]
            seeds = rng.integers(low = 0, high = 100000, size = nYRI)
            for j in range(nYRI):
                g_yri = g_yris[:,j]
                sub = (g_yri==0) & (g_nea==2)
                rng2 = np.random.default_rng(seed = seeds[j])
                acd = geno[:,rng2.choice(ceu_idx, size = ndiploids_for_dcfs, replace = False)]
                acd = np.sum(acd[sub,:], axis = 1)
                dcfs = allel.sfs(acd, n=2*ndiploids_for_dcfs)[1:-1]
                if j==0:
                    DCFS = dcfs
                else:
                    DCFS += dcfs
            dcfs = np.true_divide(DCFS, nYRI)
            FOUT.write("DCFS "+str(int(np.sum(dcfs)))+" ")
            FOUT.write(" ".join([P(x) for x in np.true_divide(dcfs, np.sum(dcfs))])+"\n")
        else:
            FOUT.write("DCFS NA "+" ".join(["NA"]*(2*ndiploids_for_dcfs-1))+"\n")

        ################################################################################################
        # D-statistic

        # block length
        if do_SE:
            if block_length.endswith("kb"):
                bl = block_length.replace("kb", "")
                if not bl.isdigit():
                    sys.exit("D> Block length must be of the form FLOATkb or INTsnp")
                bl = float(bl)
                bl = int(bl * 1000. * np.true_divide(geno.shape[0], genome_size))
                print("Block size (#SNPs):    "+str(bl))
            elif block_length.endswith("snp"):
                bl = block_length.replace("snp", "")
                if not bl.isdigit():
                    sys.exit("Block length must be of the form FLOATkb or INTsnp")
                bl = int(bl)
            else:
                sys.exit("Block length must be of the form FLOATkb or INTsnp")

        # ABBA-BABA = D statistic
        if PAR["D_do_pseudodiploid"].upper() == "YES":
            label = "D_pseudodiploid"
            ac_ceu_pd = count(geno_pd, samples2, 0)
            ac_yri_pd = count(geno_pd, samples2, 1)
            ac_nea_pd = count(geno_pd, samples2, 2)
            if not do_SE:
                num, den = allel.patterson_d(ac_ceu_pd, ac_yri_pd, ac_nea_pd, ac_ref)
                D = np.true_divide(np.sum(num), np.sum(den))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(D)+"\n")
            else:
                num, den = allel.patterson_d(ac_ceu_pd, ac_yri_pd, ac_nea_pd, ac_ref)
                d1, d2, d3, d4, d5 = allel.average_patterson_d(ac_ceu_pd, ac_yri_pd, ac_nea_pd, ac_ref, bl)
                print("Number of blocks:      "+str(len(d4)))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(d1)+" "+P(d2)+" "+P(d3)+"\n")
        else:
            label = "D_diploid"
            if not do_SE:
                num, den = allel.patterson_d(ac_ceu, ac_yri, ac_nea, ac_ref)
                D = np.true_divide(np.sum(num), np.sum(den))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(D)+"\n")
            else:
                num, den = allel.patterson_d(ac_ceu, ac_yri, ac_nea, ac_ref)
                d1, d2, d3, d4, d5 = allel.average_patterson_d(ac_ceu, ac_yri, ac_nea, ac_ref, bl)
                print("Number of blocks:      "+str(len(d4)))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(d1)+" "+P(d2)+" "+P(d3)+"\n")

        ################################################################################################
        # Fst
        # _NOTA BENE_ The Fst measure is insensitive to the monomorphic SNPs
        # As a consequence, the Fst here calculated is identical as if calculated on the SNPs polymorphic in both YRI and CEU
        num, den = allel.hudson_fst(ac_ceu, ac_yri)
        fst = np.true_divide(np.sum(num), np.sum(den))
        FOUT.write("Fst "+str(int(len(den)))+" "+P(fst)+"\n")

        # v12.2:
        if do_Fst_ascertainment:
            # for SNP only polymorphic in YRI:
            ok = np.sum(ac_yri==0, axis = 1) == 0
            num, den = allel.hudson_fst(ac_ceu[ok,:], ac_yri[ok,:])
            fst = np.true_divide(np.sum(num), np.sum(den))
            FOUT.write("Fst_Polymorphic_YRI "+str(int(len(den)))+" "+P(fst)+"\n")

        ################################################################################################
        # Output 2D-AFS
        afs2d = np.array(afs)
        afs2d[0,0], afs2d[-1,-1] = 0, 0
        afs2d = afs2d.flatten()
        FOUT.write("2D_AFS "+str(np.sum(afs2d))+" ")
        FOUT.write(" ".join([P(x) for x in np.true_divide(afs2d, np.sum(afs2d))])+"\n")

        ################################################################################################
        # Nucleotide diversity (pi) = proportion of average pairwise difference between chromosomes
        pi_ceu = np.true_divide(np.sum(allel.mean_pairwise_difference(ac_ceu)), genome_size)
        pi_yri = np.true_divide(np.sum(allel.mean_pairwise_difference(ac_yri)), genome_size)
        FOUT.write("pi_CEU_per_kb "+str(int(genome_size))+" "+str(P(pi_ceu*1e3))+"\n")
        FOUT.write("pi_YRI_per_kb "+str(int(genome_size))+" "+str(P(pi_yri*1e3))+"\n")

        FOUT.close()

    ################################################################################################
    if do_f4:
        print("\n____F4-ratio____")
        if do_stats:
            ac_vindija = ac_nea
            ac_altai = count(geno, samples2, 3)
        else:
            ac_ceu = count(geno, samples2, 0)
            ac_yri = count(geno, samples2, 1)
            ac_vindija = count(geno, samples2, 2)
            ac_altai = count(geno, samples2, 3)
            ac_ref = ac_vindija*0
            ac_ref[:,1] = 2
        #
        if PAR["D_do_pseudodiploid"].upper() == "YES":
            label = "F4ratio_pseudodiploid"
            ac_ceu_pd = count(geno_pd, samples2, 0)
            ac_yri_pd = count(geno_pd, samples2, 1)
            ac_vindija_pd = count(geno_pd, samples2, 2)
            ac_altai_pd = count(geno_pd, samples2, 3)
            f4a, _ = allel.patterson_d(ac_altai_pd, ac_ref, ac_ceu_pd, ac_yri_pd)
            f4b, _ = allel.patterson_d(ac_altai_pd, ac_ref, ac_vindija_pd, ac_yri_pd)
            # This is exactly equivalent:
            #f4a, _ = allel.patterson_d(ac_ceu_pd, ac_yri_pd, ac_altai_pd, ac_ref)
            #f4b, _ = allel.patterson_d(ac_vindija_pd, ac_yri_pd, ac_altai_pd, ac_ref)
        else:
            label = "F4ratio_diploid"
            f4a, _ = allel.patterson_d(ac_altai, ac_ref, ac_ceu, ac_yri)
            f4b, _ = allel.patterson_d(ac_altai, ac_ref, ac_vindija, ac_yri)
            # This is exactly equivalent:
            #f4a, _ = allel.patterson_d(ac_ceu, ac_yri, ac_altai, ac_ref)
            #f4b, _ = allel.patterson_d(ac_vindija, ac_yri, ac_altai, ac_ref)
        num = allel.moving_statistic(f4a, statistic=np.nansum, size=bl)
        den = allel.moving_statistic(f4b, statistic=np.nansum, size=bl)
        _, se, vj = allel.stats.misc.jackknife((num, den), statistic = lambda n, d: np.true_divide(np.nansum(n), np.nansum(d)) )
        m = np.true_divide(np.nansum(f4a), np.nansum(f4b))
        z = np.true_divide(m, se)
        with open(output+".stats", "a") as FOUT:
            FOUT.write(label+" . "+'%.5f'%(m)+" "+'%.5f'%(se)+" "+'%.5f'%(z)+"\n")

    ################################################################################################
    if do_psmc:
        PAR["psmc_r"] = str(np.true_divide(float(PAR["mut_rate"]), float(PAR["recomb_rate"])))
        print("\n____PSMC____")
        print("PSMC sample idx:    "+" ".join([str(x) for x in PAR["psmc_samples_0based"]]))
        print("PSMC bin size:      "+str(PAR["psmc_binsize"]))
        print("PSMC pattern:       "+str(PAR["psmc_pattern"]))
        print("PSMC t:             "+str(PAR["psmc_t"]))
        print("PSMC N:             "+str(PAR["psmc_N"]))
        print("PSMC r:             "+str(PAR["psmc_r"]))
        print("Generation time:    "+str(PAR["generation_time"])+"")

        for sample in PAR["psmc_samples_0based"]:
            print("Analyzing sample "+str(sample))
            sample, slabel = sample.split(":")
            sample = int(sample)
            index = np.where(samples1==sample)[0]
            if len(index)==0:
                sys.exit("The PSMC sample(s) should be present in previous CEU|YRI|NEA sample arrays.")
            else:
                index = index[0]
            outfile = output+"."+str(slabel)
            eigenstrat2psmcfa(geno, snp, index = index, seq_lengths = seq_lengths, output_prefix = outfile, bin_size = PAR["psmc_binsize"])
            cmd = bin_dir+"/"+PAR["soft_psmc"]+" -p \""+PAR["psmc_pattern"]+"\" -t"+PAR["psmc_t"]+" -N"+PAR["psmc_N"]+" -r"+PAR["psmc_r"]+" -o "+outfile+".psmc.stats "+outfile+".psmcfa"
            subprocess.run(cmd, shell=True, check=True)
            write_psmc_stats(outfile+".psmc.stats", outfile+".psmc2.stats", PAR["psmc_N"], PAR["mut_rate"], PAR["psmc_binsize"])
            os.remove(outfile+".psmcfa")
            os.remove(outfile+".psmc.stats")

    ################################################################################################
    if do_ld:
        incremental_assessment = True
        ###
        print("\n____LD____")
        if np.sum(samples2==0) < 2:
            sys.exit("--ld mode can only be used when there are at least two target diploid samples.")
        print("*** MULTI-SAMPLE MODE ***")
        print("Pseudodiploid for LD:             "+str(PAR["LD_do_pseudodiploid"]))
        print("MAC-sensitive SNP-downsampling:   "+str(PAR["do_MAC_sensitive_SNP_downsampling"]))
        print("Start (cM):                       "+str(PAR["min_cM"]))
        print("End (cM):                         "+str(PAR["max_cM"]))
        print("Step (cM):                        "+str(PAR["binsize_cM"]))
        print("Selecting haploid Neanderthal:    "+str(PAR["LD_sample_haploid_Neanderthal"]))
        print("Algo like Sankararaman12:         "+str(PAR["ld_as_sanka12"]))

        if PAR["LD_do_pseudodiploid"].upper() == "YES":
            if PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES":
               sys.exit("This configuration is not implemented yet!")
            else:
                ld_geno, ld_snp, ld_gpos = geno_pd, snp, genetic_pos
        else:
            if PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES":
                if PAR["LD_sample_haploid_Neanderthal"].upper() == "YES":
                    ld_geno, ld_snp, ld_gpos = copy.deepcopy(xs_geno), xs_snp, xs_genetic_pos
                    ld_geno[:,np.where(samples2==2)[0]] = xs_geno1[:,np.where(samples2==2)[0]]
                else:
                    ld_geno, ld_snp, ld_gpos = xs_geno, xs_snp, xs_genetic_pos
            else:
                if PAR["LD_sample_haploid_Neanderthal"].upper() == "YES":
                    ld_geno, ld_snp, ld_gpos = copy.deepcopy(geno), snp, genetic_pos
                    ld_geno[:,np.where(samples2==2)[0]] = geno1[:,np.where(samples2==2)[0]]
                else:
                    ld_geno, ld_snp, ld_gpos = geno, snp, genetic_pos

        # block sizes
        BS = [uniq(ld_snp[:,0])]
        BS.append(np.array([seq_lengths[bsn] for bsn in BS[0]]))
        n_bs = n_seq

        # ancestry ascertainment
        modes_to_analyze = ["0", "1", "2"]
        modes_to_analyze = ["0"]
        for asc_mode in modes_to_analyze:
            asc = ascertain(ld_geno,
                            method = "Sankararaman_"+asc_mode,
                            target_idx = np.where(samples2==0)[0],
                            outgroup_idx = np.where(samples2==1)[0],
                            archaic_idx = np.where(samples2==2)[0])
            g, q, p = ld_geno[asc,:], ld_snp[asc,:], ld_gpos[asc]
            g, q = g[:,np.where(samples2==0)[0]], q[:,0]
            print("\nNumber of SNPs ascertained:       "+str(g.shape[0]))

            # calculate per-chromosome
            LD_nSNP, IA = 0., []
            for block in range(n_bs):
                ok = (q==BS[0][block])
                m = np.array(average_snp_covariance(g[ok,:], q[ok], p[ok], PAR["ld_as_sanka12"], float(PAR["min_cM"]), float(PAR["max_cM"]), float(PAR["binsize_cM"])))
                LD_nSNP += np.sum(m[:,2])
                m = np.column_stack((np.array([BS[0][block]]*m.shape[0]), m))
                if block==0:
                    X = m[:,1]
                    M = m
                else:
                    if m.shape[0]!=len(X):
                        sys.exit("Error")
                    M = np.row_stack((M, m))
                # if incremental assessment
                if incremental_assessment:
                    aa, bb = LD_Y(X, M, PAR["ld_as_sanka12"])
                    IA.append(aa[1])

            # calculate jackknife
            global_estimate, global_xyfit = LD_Y(X, M, PAR["ld_as_sanka12"])
            nrmsd = NRMSD(global_xyfit[:,1], global_xyfit[:,2])
            YB = np.array([np.nan]*n_bs).astype(float)
            XYFIT = np.column_stack((np.array(["all"]*global_xyfit.shape[0]), global_xyfit))
            for block in range(n_bs):
                ok = (M[:,0].astype(int)!=BS[0][block])
                coefs, xyfit = LD_Y(X, M[ok,:], PAR["ld_as_sanka12"])
                YB[block] = coefs[1]
                xyfit = np.column_stack((np.array([str(block+1)]*xyfit.shape[0]), xyfit))
                XYFIT = np.row_stack((XYFIT, xyfit))
            if not np.isnan(global_estimate[1]):
                jk_mean, jk_se = jackknife(q_sizes = BS[1], values = YB, estimate = global_estimate[1])
            else:
                jk_mean, jk_se = np.nan, np.nan
            LD_nSNP = int(LD_nSNP)

            # export the global curve + the jackknife curves
            if verbeux:
                with gzip.open(output+"."+asc_mode+".cov.txt.gz", "wt") as FOUT:
                    FOUT.write("jackknife_run bin_cM average_cov exp_fit\n")
                    for i in range(XYFIT.shape[0]):
                        FOUT.write(XYFIT[i,0]+" "+'%.5f'%(float(XYFIT[i,1]))+" "+'%.5f'%(float(XYFIT[i,2]))+" "+'%.5f'%(float(XYFIT[i,3]))+"\n")

            # write the final statistics
            with open(output+"."+asc_mode+".cov.stats", "w") as FOUT:
                if np.isnan(global_estimate[1]) or np.isnan(jk_mean):
                    FOUT.write("A NA\nt NA\nc NA\nTime_gBSample NA\nTime_gBSample.jk.MEAN NA\nTime_gBSample.jk.SE NA\nTime_gBSample.jk.IC95.low NA\nTime_gBSample.jk.IC95.up NA\nnSNP "+str(LD_nSNP)+"\nNRMSD NA\n")
                else:
                    FOUT.write("A "+'%.6f'%global_estimate[0]+"\nt "+'%.6f'%(-1.*global_estimate[1])+"\nc "+'%.6f'%global_estimate[2]+"\nTime_gBSample "+'%.3f'%(100.*global_estimate[1])+"\nTime_gBSample.jk.MEAN "+'%.3f'%(100.*jk_mean)+"\nTime_gBSample.jk.SE "+'%.3f'%(100.*jk_se)+"\nTime_gBSample.jk.IC95.low "+'%.3f'%(100.*jk_mean-1.96*100.*jk_se)+"\nTime_gBSample.jk.IC95.up "+'%.3f'%(100.*jk_mean+1.96*100.*jk_se)+"\nn_SNP_pairs "+str(LD_nSNP)+"\nNRMSD "+'%.5f'%nrmsd+"\n")
                    if incremental_assessment:
                        FOUT.write("Incremental.Assessment.Time_gBSample "+" ".join(['%.3f'%(100.*aa) for aa in IA])+"\n")


    ################################################################################################
    if do_ld1:
        print("\n____LD____")
        if np.sum(samples2==0) != 1:
            sys.exit("--ld1 mode can only be used when there is a single diploid sample.")
        print("*** SINGLE-SAMPLE MODE ***")
        print("In this mode, the MAC-downsampling cannot be performed and was therefore set to NO.")
        print("Pseudodiploid for LD:             "+str(PAR["LD_do_pseudodiploid"]))
        print("MAC-sensitive SNP-downsampling:   NO")
        print("Start (cM):                       "+str(PAR["min_cM"]))
        print("End (cM):                         "+str(PAR["max_cM"]))
        print("Step (cM):                        "+str(PAR["binsize_cM"]))

        if PAR["LD_do_pseudodiploid"].upper() == "YES":
            ld_geno, ld_snp, ld_gpos = geno_pd, snp, genetic_pos
        else:
            ld_geno, ld_snp, ld_gpos = geno, snp, genetic_pos

        # ancestry ascertainment
        asc = ascertain(ld_geno,
                        method = "Moorjani_0",
                        target_idx = np.where(samples2==0)[0],
                        outgroup_idx = np.where(samples2==1)[0],
                        archaic_idx = np.where(samples2==2)[0])
        g, q, p = ld_geno[asc,:], ld_snp[asc,:], ld_gpos[asc]
        g, q = g[:,np.where(samples2==0)[0]], q[:,0]

        # calculate covariance
        IA = []
        global_estimate, jk_mean, jk_se, XYFIT = single_sample_covariance(PAR, output, g, q, p, bin_dir, seq_lengths, float(PAR["min_cM"]), float(PAR["max_cM"]), float(PAR["binsize_cM"]))
        nrmsd = NRMSD(XYFIT[XYFIT[:,0]=="all",2].astype(float), XYFIT[XYFIT[:,0]=="all",3].astype(float))
        LD_nSNP = "NA"

        # export the global curve + the jackknife curves
        if verbeux:
            with gzip.open(output+".cov.txt.gz", "wt") as FOUT:
                FOUT.write("jackknife_run bin_cM average_cov exp_fit\n")
                for i in range(XYFIT.shape[0]):
                    FOUT.write(XYFIT[i,0]+" "+'%.5f'%(float(XYFIT[i,1]))+" "+'%.5f'%(float(XYFIT[i,2]))+" "+'%.5f'%(float(XYFIT[i,3]))+"\n")

        # write the final statistics
        with open(output+".cov.stats", "w") as FOUT:
            if np.isnan(global_estimate[1]) or np.isnan(jk_mean):
                FOUT.write("A NA\nt NA\nc NA\nTime_gBSample NA\nTime_gBSample.jk.MEAN NA\nTime_gBSample.jk.SE NA\nTime_gBSample.jk.IC95.low NA\nTime_gBSample.jk.IC95.up NA\nnSNP 0\nNRMSD NA\n")
            else:
                FOUT.write("A "+'%.6f'%global_estimate[0]+"\nt "+'%.6f'%(-1.*global_estimate[1])+"\nc "+'%.6f'%global_estimate[2]+"\nTime_gBSample "+'%.3f'%(100.*global_estimate[1])+"\nTime_gBSample.jk.MEAN "+'%.3f'%(100.*jk_mean)+"\nTime_gBSample.jk.SE "+'%.3f'%(100.*jk_se)+"\nTime_gBSample.jk.IC95.low "+'%.3f'%(100.*jk_mean-1.96*100.*jk_se)+"\nTime_gBSample.jk.IC95.up "+'%.3f'%(100.*jk_mean+1.96*100.*jk_se)+"\nn_SNP_pairs "+str(LD_nSNP)+"\nNRMSD "+'%.5f'%nrmsd+"\n")


    ################################################################################################
    if do_sprime:
        print("\n____S'____")
        print("Min score:                        "+str(PAR["Sprime_minscore"]))
        print("MAC-sensitive SNP-downsampling:   "+str(PAR["do_MAC_sensitive_SNP_downsampling"]))

        if PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES":
            G_NEA = generate_Sprime_input(xs_geno, xs_snp, xs_genetic_pos, ind, samples2, output)
        else:
            G_NEA = generate_Sprime_input(geno, snp, genetic_pos, ind, samples2, output)

        # run S'
        cmd = "java -jar "+bin_dir+"/"+PAR["Sprime"]+" gt="+output+".vcf outgroup="+output+".YRI.ind"+" map="+output+".map out="+output+" mu="+str(PAR["mut_rate"])+" minscore="+str(PAR["Sprime_minscore"])+" >/dev/null"
        subprocess.call(cmd, shell=True)

        # read the scores
        sp = pd.read_csv(output+".score", sep = "\t", usecols = [0,1,5,6], dtype = np.int32, na_filter = None, header = 0).to_numpy()
        RS = pd.read_csv(output+".score", sep = "\t", usecols = [2], dtype = str, na_filter = None, header = 0).to_numpy()[:,0]
        segs = np.unique(sp[:,2])

        # compute statistics
        for i, seg in enumerate(segs):
            x = sp[sp[:,2]==seg,:]
            if len(np.unique(x[:,0]))!=1:
                sys.exit("Error. Sprime. A putative introgressed segments spans over more than one chromosome.")
            rs = RS[sp[:,2]==seg]
            # match rate
            g_nea = np.array([G_NEA[r] for r in rs])
            mr = np.abs(x[:,3]*2 - g_nea)
            mr = np.true_divide(np.sum(mr<=1), len(mr))
            # length
            s, e = np.min(x[:,1]), np.max(x[:,1])
            ls = e-s+1
            # concat
            ad = [x.shape[0], ls, mr, s, e, x[0,0]]
            if i==0:
                if len(segs)==1:
                    X = np.array([ad])
                else:
                    X = np.array(ad)
            else:
                X = np.row_stack((X, np.array(ad)))

        # outputs
        with open(output+".sprime2.stats", "w") as fS:
            fS.write("MinFragmentMatchRate nFragments DetectionRate NonOverlapDetectionRate MatchRateOverAllFragments LengthBpMean LengthBpQ2.5 LengthBpQ97.5 MatchRateMean MatchRateQ2.5 MatchRateQ97.5 LengthBpMedian MatchRateMedian\n")
            for minMR in [0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
                if len(segs)==0:
                    fS.write('%.2f'%(minMR)+" 0 "+" ".join(["nan"]*11)+"\n")
                else:
                    x, all_l, all_m = get_sprime_stats(X, minMR, genome_size, 3)
                    if verbeux and minMR==0.:
                        with gzip.open(output+".sprime2.txt.gz", "wt") as fT:
                            fT.write("\n".join([str(int(ala))+" "+'%.2f'%(ama) for (ala,ama) in zip(all_l,all_m)])+"\n")
                    fS.write('%.2f'%(x[0])+" "+str(int(x[1]))+" ")
                    fS.write(" ".join(['%.3f'%(z) for z in x[2:5]])+" ")
                    fS.write(" ".join([str(np.round(z)) for z in x[5:8]])+" ")
                    fS.write(" ".join(['%.3f'%(z) for z in x[8:]])+"\n")

        # cleaning
        os.remove(output+".vcf")
        os.remove(output+".map")
        os.remove(output+".YRI.ind")
        os.remove(output+".score")
        os.remove(output+".log")


    ################################################################################################
    if do_crf:
        print("\n____CRF____")
        print("Min probability:                  "+str(PAR["CRF_threshold"]))
        print("CRF par file:                     "+str(PAR["CRF_parFile"]))
        print("MAC-sensitive SNP-downsampling:   "+str(PAR["do_MAC_sensitive_SNP_downsampling"]))
        print("Minimum segment length (cM):      "+str(PAR["CRF_min_segment_length_cM"]))
        # note that by default in Sankararaman 14, CRF_min_segment_length_cM = 0.02 cM

        if os.path.isabs(output):
            sys.exit("When using CRF, you must provide an output path that is relative.")

        def crf_cleaning(output):
            rm_if_exists(output+".CEU.ind")
            rm_if_exists(output+".CRF.par")
            rm_if_exists(output+".CEU.geno")
            rm_if_exists(output+".YRI.geno")
            rm_if_exists(output+".YRI.freq")
            rm_if_exists(output+".NEA.geno")
            rm_if_exists(output+".NEA.freq")
            rm_if_exists(output+".SNP.snp")
            rm_if_exists(output+".crf")
            rm_if_exists(output+".crf.filtered.snps")
            rm_if_exists(output+".f01")
            rm_if_exists(output+".f02")

        if PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES":
            seqs = uniq(xs_snp[:,0])
        else:
            seqs = uniq(snp[:,0])

        SEGMENTS, N_ARCHAIC, N_ALL = [], np.array([0]*nCEU*2), np.array([0]*nCEU*2)
        # `SEGMENTS` is the vector of min.cM-filtered archaic segment lengths detected across all chromosomes and all individuals
        # `N_ARCHAIC` is the vector of #introgressed SNPs per individual
        # `N_ALL` is the vector of #SNPs per individual
        print("\n")
        for seq in seqs:
            crf_cleaning(output)
            if PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES":
                generate_CRF_input(PAR, xs_geno, xs_geno1, xs_snp, xs_genetic_pos, ind, samples2, seq, output, bin_dir, seq_lengths)
            else:
                generate_CRF_input(PAR, geno, geno1, snp, genetic_pos, ind, samples2, seq, output, bin_dir, seq_lengths)
            cmd = bin_dir+"/"+PAR["CRF_script"]+" -p "+output+".CRF.par >/dev/null"
            subprocess.call(cmd, shell=True)
            o_pos = pd.read_csv(output+".SNP.snp", sep = "\s+", usecols = [2,3], engine = "c", dtype = str, header = None).to_numpy()
            o_gpos = o_pos[:,0].astype(float)
            o_pos = o_pos[:,1].astype(int)
            prob = pd.read_csv(output+".crf", sep = "\t", dtype = float, header = None).to_numpy()[:,range(nCEU*2)]
            for j in range(prob.shape[1]):
                    # below threshold encoded as 1 -- above, as 0
                    pnaa = (prob[:,j] < float(PAR["CRF_threshold"])).astype(int)
                    # introgression rate
                    N_ARCHAIC[j] += (prob.shape[0] - np.sum(pnaa))
                    N_ALL[j] += prob.shape[0]
                    # putative archaic fragment length
                    pas = zero_runs(pnaa)
                    if len(pas) >= 1:
                        for i in range(pas.shape[0]):
                            length_bp = o_pos[pas[i,1]-1] - o_pos[pas[i,0]]
                            length_cM = 100. * (o_gpos[pas[i,1]-1] - o_gpos[pas[i,0]])
                            if length_cM < float(PAR["CRF_min_segment_length_cM"]):
                                continue
                            else:
                                SEGMENTS.append(length_bp)

        alpha = np.true_divide(N_ARCHAIC, N_ALL)
        crf_cleaning(output)

        # exporting
        with open(output+".crf.stats", "w") as fS:
            a, b = summary(np.array(SEGMENTS))
            fS.write("Variable "+" ".join([str(x) for x in a])+"\n")
            if b[0]==0:
                fS.write("Putative_Archaic_Segment_Length_bp_Across_Inds "+" ".join(["NA" for x in b])+"\n")
            else:
                fS.write("Putative_Archaic_Segment_Length_bp_Across_Inds "+" ".join([str(int(x)) for x in b])+"\n")
            a, b = summary(alpha)
            if b[0]==0:
                fS.write("Ancestry_Proportion_Across_Inds "+" ".join(["NA" for x in b])+"\n")
            else:
                fS.write("Ancestry_Proportion_Across_Inds "+str(int(b[0]))+" "+" ".join(['%.5f'%(x) for x in b[1:]])+"\n")
        if verbeux:
            with gzip.open(output+".crf.txt.gz", "wt") as fS:
                fS.write("\n".join(np.array(SEGMENTS).astype(np.int64).astype(str))+"\n")

    ################################################################################################
    print("\n____End____")
    print("Analyzed:   "+prefix)
    print("Output:     "+output)
    print(str(np.round((time.time()-start_time)/60, 1))+" min\n")


def LD_Y(X, m, ld_as_sanka12):
    if ld_as_sanka12:
        Y = np.zeros(len(X)).astype(float)
        for i in range(len(X)):
            ok = (m[:,1]==X[i])
            Y[i] = np.true_divide(np.sum(m[ok,2]), np.sum(m[ok,3]))
        coefs, xyfit = expfit(X, Y, isbad = None)
        xyfit[np.isnan(xyfit[:,1]),1] = 0.
    else:
        Y = np.array([np.nan]*len(X)).astype(float)
        for i in range(len(X)):
            ok = (m[:,1]==X[i])
            Y[i] = np.true_divide(np.sum(m[ok,2]), np.sum(m[ok,3]))
        bad = np.isnan(Y) | np.isinf(Y)
        if np.sum(bad)>=(0.5*len(X)):
            return [np.nan]*3, np.column_stack((X, np.array([np.nan]*len(X)), np.array([np.nan]*len(X))))
        coefs, xyfit = expfit(X, Y, isbad = bad)
    return coefs, xyfit

# Ascertain "ancestry-informative" SNPs
def ascertain(g, method, target_idx, outgroup_idx, archaic_idx, max_DAF = 0.1):
    # returns a vector where True = "SNP checks the ascertainment scheme", False otherwise
    print("Ascertainment mode:               "+method)

    if method == "Sankararaman_0":
        # Ascertainment scheme used in Sankararaman et al. 2012
        """
        We pick SNPs that are derived in N (at least one of the reads that maps to the SNP carries
        the derived allele), are polymorphic in E and have a derived allele frequency in E < 0.1. This
        ascertainment enriches for SNPs that arose in the N lineage and introgressed into E (in addition to
        SNPs that are polymorphic in the N H ancestor and are segregating in the present-day population).
        We chose a cutoff of 0.10 based on an analysis that computes the excess of the number of sites
        where Neandertal carries the derived allele compared to the number of sites where Denisova carries
        the derived allele stratified by the derived allele frequency in European populations ( (n−d)
        where s is the total number of polymorphic SNPs in Europeans).
        """
        a1 = np.sum(g[:,archaic_idx], axis = 1) >= 1
        a2 = np.true_divide(np.sum(g[:,target_idx], axis = 1), 2.*len(target_idx))
        a2 = (a2>0.) & (a2<float(max_DAF))
        asc = a1 & a2

    elif method == "Sankararaman_1":
        """
        Ascertainment 1: SNPs for which Neandertal carries a derived allele, E is polymorphic and Y does not carry a derived allele.
        """
        a1 = np.sum(g[:,archaic_idx], axis = 1) >= 1
        a2 = np.sum(g[:,target_idx], axis = 1) > 0
        a3 = np.sum(g[:,outgroup_idx], axis = 1) == 0
        asc = a1 & a2 & a3

    elif method == "Sankararaman_2":
        """
        Ascertainment 2: SNPs for which Neandertal carries a derived allele, E is polymorphic and
        Y does not carry a derived allele and SNPs for which Neandertal carries a derived allele,
        E does not carry a derived allele and Y is polymorphic.
        """
        a1 = np.sum(g[:,archaic_idx], axis = 1) >= 1
        a2 = np.sum(g[:,target_idx], axis = 1) > 0
        a3 = np.sum(g[:,outgroup_idx], axis = 1) == 0
        a4 = np.sum(g[:,target_idx], axis = 1) == 0
        a5 = np.sum(g[:,outgroup_idx], axis = 1) > 0
        asc = (a1 & a2 & a3) | (a1 & a4 & a5)

    elif method == "Moorjani_0":
        # Ascertainment scheme used in Moorjani et al. 2016 ("ascertainment 0"):
        """
        That is, sites where Neanderthals carry at least one
        derived allele (relative to chimpanzees) and all individuals in a
        panel of sub-Saharan Africans [which have little or no evidence of
        Neanderthal ancestry (18) carry the ancestral allele (17)]. We chose
        this ascertainment (referred to as “ascertainment 0”) because it
        minimizes the signal of background correlation, while amplifying
        the signal of Neanderthal ancestry (9).
        """
        a1 = np.sum(g[:,archaic_idx], axis = 1) >= 1
        a2 = np.sum(g[:,outgroup_idx], axis = 1) == 0
        asc = a1 & a2

    elif method == "Archaic_array":
        """
        At least one Neanderthal allele differs from the majority allele
        in a panel of 24 YRI samples (Fu et al 2015, Moorjani et al 2016)
        """
        n_YRI = 24
        outgroup_idx = outgroup_idx[range(n_YRI)]
        daf_outgroup = np.true_divide(np.sum(g[:,outgroup_idx], axis = 1), 2.*n_YRI)
        daf_nea = np.true_divide(np.sum(g[:,archaic_idx], axis = 1), 2.*len(archaic_idx))
        asc = np.repeat(False, g.shape[0], axis = 0)
        asc1 = (daf_outgroup>0.5) & (daf_nea<1.)
        asc2 = (daf_outgroup<0.5) & (daf_nea>0.)
        asc = asc | asc1 | asc2

    else:
        sys.exit("Not implemented.")

    return asc

#_____________________________________________________________________________#
#_____________________________________________________________________________#
#_____________________________________________________________________________#
#                               EXTRAS v6.5                                   #
#_____________________________________________________________________________#
#_____________________________________________________________________________#
#_____________________________________________________________________________#

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# Replace genotypes encoded as 2 into 0, and encoded as 0 into 2
def invert_ancestrality(x):
    y = copy.deepcopy(x)
    dix = {0:2, 2:0, 1:1, 9:9}
    return np.vectorize(dix.get)(y)

# Quickly read fixed-width file (typically a .geno file)
def read_fixedwidth(file, dtype, columns):
    X = pd.read_fwf(file, colspecs=[(x,x+1) for x in columns], dtype = dtype, header = None, na_filter = False)
    X = X.to_numpy()
    return X

# Reshape EIGENSTRAT matrices
def reshape_EIGENSTRAT(x, isgeno = False):
    if len(x.shape) == 1:
        x = np.array([x])
    if isgeno == True and x.shape[0] == 1:
        x = np.transpose(x)
    return x

# Read an EIGENSTRAT file
def read_eigenstrat(prefix,
                    read_geno = True, read_snp = True, read_ind = True,
                    columns = None,
                    ancestral_allele_encoded_as_0 = True,
                    verbose = False):

    if type(columns) == type(None):
        read_specific_cols = False
    else:
        read_specific_cols = True
    if type(columns) == type(list()):
        columns = np.array(columns)

    if path.exists(prefix+".snp") and path.exists(prefix+".snp.gz"):
        sys.exit("There cannot be both a .snp and a .snp.gz file.")
    if path.exists(prefix+".geno") and path.exists(prefix+".geno.gz"):
        sys.exit("There cannot be both a .geno and a .geno.gz file.")
    if path.exists(prefix+".snp.gz"):
        snp_outfix = ".snp.gz"
    else:
        snp_outfix = ".snp"
    if path.exists(prefix+".geno.gz"):
        geno_outfix = ".geno.gz"
    else:
        geno_outfix = ".geno"

    if verbose:
        print(prefix)
        print("Read geno: "+str(read_geno))
        print("Read snp: "+str(read_snp))
        if read_specific_cols:
            print("Read ind: "+str(True)+" (to indicate selected sample names)")
        else:
            print("Read ind: "+str(read_ind))
        if read_specific_cols:
            print("Reading columns: "+", ".join([str(x) for x in columns]))

    if read_geno:
        if ancestral_allele_encoded_as_0 == True:
            print("[anc:=0] [na:=9]")
        else:
            print("[anc:=1 => anc:=0] [na:=9]")

    if read_specific_cols:
        read_ind = True

    no_ind = True
    if read_ind:
        ind = pd.read_csv(prefix+".ind", sep = "\s+", usecols = [0,2], engine = "c", dtype = str, na_filter = None, header = None).to_numpy()
        no_ind = False
    else:
        ind = None

    if read_specific_cols:
        inds = ind[columns,:]
        df = np.column_stack((columns, inds[:,0], inds[:,1]))
        df = pd.DataFrame(df)
        df.columns = ["0.Index.IndFile", "SampleName.IndFile", "Population"]
        print("_________________________________________________\n")
        print(df)
        print("_________________________________________________")
        del df

    if read_snp:
        snp = pd.read_csv(prefix+snp_outfix, sep = "\s+", usecols = [1,2,3], names = ["chrom", "gpos", "ppos"], engine = "c", dtype = {"chrom": str, "gpos": str, "ppos": str}, na_filter = None, header = None, float_precision = 'round_trip', nrows = 1)
        if snp["gpos"][0] == ".":
            snp = pd.read_csv(prefix+snp_outfix, sep = "\s+", usecols = [1,3], names = ["chrom", "ppos"], engine = "c", dtype = {"chrom": np.int32, "ppos": np.int32}, na_filter = None, header = None, float_precision = 'round_trip')
            genetic_pos = None
            snp = snp.to_numpy()
        else:
            snp = pd.read_csv(prefix+snp_outfix, sep = "\s+", usecols = [1,2,3], names = ["chrom", "gpos", "ppos"], engine = "c", dtype = {"chrom": np.int32, "gpos": str, "ppos": np.int32}, na_filter = None, header = None, float_precision = 'round_trip')
            genetic_pos = np.array([Decimal(x) for x in snp["gpos"]])
            snp = snp[["chrom", "ppos"]].to_numpy()
    else:
        genetic_pos, snp = None, None

    if read_geno:
        if no_ind:
            if geno_outfix.endswith(".gz"):
                with gzip.open(prefix+geno_outfix, "rt") as f:
                    first_line = f.readline().strip()
            else:
                with open(prefix+geno_outfix) as f:
                    first_line = f.readline().strip()
            nind = len(first_line)
        else:
            nind = ind.shape[0]
        if geno_outfix.endswith(".gz"):
            FGENO = gzip.open(prefix+geno_outfix, "rt")
        else:
            FGENO = open(prefix+geno_outfix, "r")
        if read_specific_cols:
            geno = read_fixedwidth(FGENO, dtype = np.int8, columns = columns)
        else:
            geno = read_fixedwidth(FGENO, dtype = np.int8, columns = list(range(nind)))
        FGENO.close()
        geno = reshape_EIGENSTRAT(geno, isgeno = True)
    else:
        geno = None

    if no_ind == False and read_specific_cols:
        ind = ind[columns,]

    if read_geno and ancestral_allele_encoded_as_0 == False:
        geno = invert_ancestrality(geno)

    # checks
    if read_geno==True and read_snp==True:
        if geno.shape[0] != snp.shape[0]:
            sys.exit("Not equal number of SNPs")
    if read_geno==True and read_ind==True:
        if geno.shape[1] != ind.shape[0]:
            sys.exit("Not equal number of individuals")

    return ind, snp, genetic_pos, geno

# unique() but keep the order of elements (by default unique sorts them)
def uniq(x):
    indexes = np.unique(x, return_index=True)[1]
    return np.array([x[index] for index in sorted(indexes)])

# Compute the average weighted covariance decay as in Moorjani et al 2012; Sankararaman et al. 2012
def average_snp_covariance(geno, chrom, gpos, ld_as_sanka12, min_cM = 0.02, max_cM = 1.0, binsize_cM = 1e-3):

    # geno is a 2d-array of genotypes, with SNPs across rows, and individuals across columns
    # chrom is an array giving the chromosome of the SNPs
    # gpos is an array giving the genetic positions of the SNPs in centiMorgans (cM)
    # _NOTA BENE_ Handling of missing data (9) has not been implemented yet!

    bins = np.arange(min_cM, max_cM, step = binsize_cM)
    cov_sum, cov_n = np.repeat(float(0), len(bins)), np.repeat(float(0), len(bins))
    seqs = uniq(chrom)
    for seq in seqs:
        sub = np.where(chrom==seq)[0]
        if len(sub) <= 2:
            continue
        g, p = geno[sub,:], gpos[sub]
        p = p.reshape(p.shape[0], -1)
        try:
            gdist = pdist(p)
        except:
            #print("/!\ Warning! Error when using Decimal numbers. Converting to float.")
            pp = copy.deepcopy(p)
            pp = pp.astype(float)
            gdist = pdist(pp)
        # covariance
        if ld_as_sanka12:
            cov_bias = True
        else:
            cov_bias = False
        covar = np.cov(g, bias = cov_bias)
        covar = covar[np.triu_indices(covar.shape[0], k = 1)]
        # distance binning
        if ld_as_sanka12:
            ok = np.where((gdist>=(min_cM-binsize_cM)) & (gdist<(max_cM-binsize_cM)))[0]
            covar, gdist = covar[ok], gdist[ok]
            #gdist = np.floor(np.true_divide(gdist-min_cM, binsize_cM)).astype(np.int32) + np.int32(1)
            gdist = np.floor(np.true_divide(gdist-min_cM+binsize_cM, binsize_cM)).astype(np.int32)
        else:
            ok = np.where((gdist>=min_cM) & (gdist<max_cM))[0]
            covar, gdist = covar[ok], gdist[ok]
            gdist = np.floor(np.true_divide(gdist-min_cM, binsize_cM)).astype(np.int32)
        np.add.at(cov_sum, gdist, covar)
        np.add.at(cov_n, gdist, 1)

    cov_mean = np.true_divide(cov_sum, cov_n)

    return np.column_stack((bins, cov_sum, cov_n, cov_mean))

# Delete-m jackknife
def jackknife(q_sizes, values, estimate):
    # Busing et al. 1999
    # `q_sizes` is a vector containing the sizes of each jackknife block
    # `value` is a vector containing the values of the statistic computed for each jackknife block
    # `estimate` is the value of the statistic calculated over the whole data

    if len(q_sizes) != len(values):
        sys.exit("q_sizes should have the same length as value")

    tot_size = np.sum(q_sizes)
    n_q = float(len(values))
    a = np.sum(values * (1. - np.true_divide(q_sizes, tot_size)))
    mean = n_q * estimate - a
    h = np.true_divide(tot_size, q_sizes)
    b = h*estimate - (h-1.)*values - n_q*estimate + a
    var = np.true_divide(1, n_q) * np.sum( (h-1.)**-1 * b**2 )
    sd = np.sqrt(var)

    return mean, sd

# Jacquelin exponential fit
def expfit(x, y, isbad):
    ori_x = copy.deepcopy(x)
    ori_y = copy.deepcopy(y)
    if isbad is not None:
        x = x[isbad==False]
        y = y[isbad==False]
    def exp1(x, a, b, c):
        return(a + b * np.exp(c*x))
    try:
        n = len(x)
        S = [0]
        for k in np.arange(1, n):
            S += [ S[-1] + 1/2 * (y[k]+y[k-1])*(x[k]-x[k-1]) ]
        S, x, y = np.array(S), np.array(x), np.array(y)
        M1 = [[sum( (x-x[0])**2 ), sum( (x-x[0])*S )],[sum( (x-x[0])*S ), sum( S**2 )]]
        M2 = [sum( (y-y[0])*(x-x[0]) ), sum( (y-y[0])*S )]
        M1, M2 = np.array(M1), np.array(M2)
        K = np.dot(np.linalg.inv(M1), M2)
        c = K[1]
        N1 = [[n, sum( np.exp(c*x) )], [sum( np.exp(c*x) ), sum( np.exp(2*c*x) )]]
        N2 = [sum( y ), sum( y*np.exp(c*x) )]
        AB = np.dot(np.linalg.inv(N1), N2)
        a, b = AB[0], AB[1]
        if np.isnan(a):
            return np.array([np.nan]*3), np.column_stack((ori_x, ori_y, np.array([np.nan]*len(ori_x))))
        try:
            popt, pcov = curve_fit(exp1, x, y, maxfev = 5000, p0 = [a, b, c])
            yfit = exp1(x, popt[0], popt[1], popt[2])
            ori_x[isbad==False] = x
            ori_y[isbad==False] = y
            xy = np.column_stack((ori_x, ori_y, yfit))
            return np.array([popt[1], -1*popt[2], popt[0]]), xy
        except:
            return np.array([np.nan]*3), np.column_stack((ori_x, ori_y, np.array([np.nan]*len(ori_x))))
    except:
        return np.array([np.nan]*3), np.column_stack((ori_x, ori_y, np.array([np.nan]*len(ori_x))))

# Function to calculate the ancestry LD curve using Moorjani (2016) dating method.
def single_sample_covariance(PAR, output, geno, chrom, gpos, bin_dir, seq_lengths, min_cM = 0.02, max_cM = 1.0, binsize_cM = 1e-3):
    seqs = uniq(chrom)
    if geno.shape[1] != 1:
        sys.exit("When using --ld1, you must provide one single sample.")
    with open(output+".LD.geno", "w") as F:
        for i in range(geno.shape[0]):
            # write with the proper EIGENSTRAT format where 2 = ancestral/ancestral, 0 = derived/derived
            F.write(str(2-geno[i,0])+"\n")
    with open(output+".LD.snp", "w") as F:
        for i in range(chrom.shape[0]):
            F.write(". "+str(chrom[i])+" "+str(gpos[i]*Decimal('0.01'))+" . A T\n")
    with open(output+".LD.ind", "w") as F:
        F.write("1 U 1")
    with open(output+".LD.par", "w") as F:
        F.write("genotypename: "+output+".LD.geno\n")
        F.write("snpname: "+output+".LD.snp\n")
        F.write("indivname:  "+output+".LD.ind\n")
        F.write("maxdis: "+str(max_cM*0.01)+"\n")
        F.write("binsize: "+str(binsize_cM*0.01)+"\n")
        F.write("idname: 1\n")
        F.write("output: "+output+".LD.out\n")
        F.write("jackknife: YES\n")
        F.write("numchrom: "+str(int(np.max(seqs)))+"\n")
    cmd = bin_dir+"/"+PAR["LD_Moorjani_bin"]+" -p "+output+".LD.par"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    # global estimate
    xy = pd.read_csv(output+".LD.out", sep = "\t", dtype = float, na_filter = None, header = None, comment = "#").to_numpy()
    xy = xy[np.logical_and(xy[:,0]>=min_cM, xy[:,0]<max_cM),:]
    global_estimate, global_xyfit = expfit(xy[:,0], xy[:,1], isbad = None)
    XYFIT = np.column_stack((np.array(["all"]*xy.shape[0]), global_xyfit))
    # jackknife-specific estimates
    COEFS = []
    for seq in seqs:
        #N_SNPs.append(np.sum(chrom==seq))
        xy = pd.read_csv(output+".LD.out:"+str(seq), sep = "\t", dtype = float, na_filter = None, header = None, comment = "#").to_numpy()
        xy = xy[np.logical_and(xy[:,0]>=min_cM, xy[:,0]<max_cM),:]
        coefs, xyfit = expfit(xy[:,0], xy[:,1], isbad = None)
        COEFS.append(coefs)
        xyfit = np.column_stack((np.array([seq]*xy.shape[0]), xyfit))
        XYFIT = np.row_stack((XYFIT, xyfit))
    COEFS = np.array(COEFS)
    # calculate jackknife statistics
    if not np.isnan(global_estimate[1]):
        jk_mean, jk_se = jackknife(q_sizes = np.array([seq_lengths[sns] for sns in seqs]), values = COEFS[:,1], estimate = global_estimate[1])
    else:
        jk_mean, jk_se = np.nan, np.nan
    for suffix in ["geno", "snp", "ind", "par", "out"]+["out:"+str(seq) for seq in seqs]:
        os.remove(output+".LD."+suffix)
    return global_estimate, jk_mean, jk_se, XYFIT

if __name__=="__main__":
    main()
