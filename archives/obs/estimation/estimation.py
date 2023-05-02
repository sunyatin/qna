#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 20:00:13 2021
@author: rtournebize
"""

import sys
sys.path.append('~/gitdir/qna/scripts')
from sumstats import get_sprime_stats, summary
import pandas as pd
import numpy as np


#########################################################################################################
#########################################################################################################
###################################              Sprime

prefix = "~/gitdir/qna/archives/obs/estimation/Sprime/CEU.chr"

first = True
genome_size = 0
for chrom in range(1,23):
    print(chrom)
    
    # read the scores
    sp = pd.read_csv(prefix+str(chrom)+".ND_match", sep = "\t", usecols = [0,1,5,6,7], dtype = np.int32, na_filter = None, header = 0).to_numpy()
    Match = pd.read_csv(prefix+str(chrom)+".ND_match", sep = "\t", usecols = [8], dtype = str, na_filter = None, header = 0).to_numpy()
    segs = np.unique(sp[:,2])
    
    genome_size += np.max(sp[:,1]) + 1 # note that this is a very rough estimate of the actual size of the genome analyzed
    
    # compute statistics
    for i, seg in enumerate(segs):
        idx = sp[:,2]==seg
        x = sp[idx,:]
        match = Match[idx,:]
        if len(np.unique(x[:,0]))!=1:
            sys.exit("Error. Sprime. A putative introgressed segments spans over more than one chromosome.")
        #rs = RS[sp[:,2]==seg]
        # match rate
        #g_nea = np.array([G_NEA[r] for r in rs])
        #mr = np.abs(x[:,3]*2 - g_nea)
        #mr = np.true_divide(np.sum(mr<=1), len(mr))
        n1 = np.sum(Match=="match")
        n0 = np.sum(Match=="mismatch")
        if n0+n1 == 0:
            mr = 0.
        else:
            mr = np.true_divide(n1, n0+n1)
        # length
        s, e = np.min(x[:,1]), np.max(x[:,1])
        ls = e-s+1
        # concat
        ad = [x.shape[0], ls, mr, s, e, x[0,0]]
        if first:
            if len(segs)==1:
                X = np.array([ad])
            else:
                X = np.array(ad)
            first = False
        else:
            X = np.row_stack((X, np.array(ad)))
    
# outputs
with open("~/gitdir/qna/archives/obs/estimation/obs.sprime2.stats", "w") as fS:
    fS.write("MinFragmentMatchRate nFragments DetectionRate NonOverlapDetectionRate MatchRateOverAllFragments LengthBpMean LengthBpQ2.5 LengthBpQ97.5 MatchRateMean MatchRateQ2.5 MatchRateQ97.5 LengthBpMedian MatchRateMedian\n")
    for minMR in [0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        if len(segs)==0:
            fS.write('%.2f'%(minMR)+" 0 "+" ".join(["nan"]*9)+"\n")
        else:
            x = get_sprime_stats(X, minMR, genome_size, 3)
            fS.write('%.2f'%(x[0])+" "+str(int(x[1]))+" ")
            fS.write(" ".join(['%.3f'%(z) for z in x[2:5]])+" ")
            fS.write(" ".join([str(np.round(z)) for z in x[5:8]])+" ")
            fS.write(" ".join(['%.3f'%(z) for z in x[8:]])+"\n")


#########################################################################################################
#########################################################################################################
###################################              CRF

prefix = "~/gitdir/qna/archives/obs/estimation/CRF/CEU.hapmap/summaries/haplotypes/chr-"
output = "~/gitdir/qna/archives/obs/estimation/obs.crf.stats"


CRF_min_segment_length_cM = 0.02
CRF_threshold = 0.90

SEGMENTS, ALPHA = [], []
for chrom in range(1,23):
    print(chrom)
    
    df = pd.read_csv(prefix+str(chrom)+".thresh-90.length-0.00.haplotypes.gz", sep = "\s+", usecols = [1,4,5,6], dtype = np.float, header = None, na_filter = None, comment = "#").to_numpy()
    #     ind    length_bp    length_cM    proba
    
    df = df[df[:,2]>=CRF_min_segment_length_cM,:]
    df = df[df[:,3]>=CRF_threshold,:]
    
    Inds = np.unique(df[:,0])
    
    for ind in Inds:
        X = df[df[:,0]==ind,:]
        for x in X:
            SEGMENTS.append(x[1])
    
ALPHA = np.array([]) # Cannot be really determined from the files provided, cf. article

# exporting
with open(output, "w") as fS:
    a, b = summary(np.array(SEGMENTS))
    fS.write("Variable "+" ".join([str(x) for x in a])+"\n")
    if b[0]==0:
        fS.write("Putative_Archaic_Segment_Length_bp_Across_Inds "+" ".join(["NA" for x in b])+"\n")
    else:
        fS.write("Putative_Archaic_Segment_Length_bp_Across_Inds "+" ".join([str(int(x)) for x in b])+"\n")
    a, b = summary(ALPHA)
    if b[0]==0:
        fS.write("Ancestry_Proportion_Across_Inds "+" ".join(["NA" for x in b])+"\n")
    else:
        fS.write("Ancestry_Proportion_Across_Inds "+str(int(b[0]))+" "+" ".join(['%.5f'%(x) for x in b[1:]])+"\n")


        
                
                
