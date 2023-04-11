#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> ms command of Model C.1 (which has the highest goodness-of-fit)
Gene flow into the modern human ancestor from a lineage that splits from the modern human ancestor 24K generations ago
Durvasula et al. (2020)

"In model C.1, we simulate gene flow from a branch population separates from the modern human and archaic
common ancestor prior to their split from each other (696Ky B.P.) and then introgresses into the modern-
human ancestor (63.8Ky B.P; Figure S2). We find that this model can produce the U-shaped spectrum we
observe in the data, particularly at introgression fractions between 0.05-0.1 (Figure S27).
We found the admixture time was 2,700 generations B.P. (95% HPD: 2,200-3,700), the admixture fraction
was 0.05 (95% HPD: 0.03-0.07), the split time was 26,000 generation B.P. (95% HPD: 24,700-27,000), and
the effective population size was 26,000 (95% HPD: 22,000-29,000). Simulations from the best fit values
resulted in a KS p−value of 0.09."

"We use a mutation rate (per basepair) of 1.2 × 10−8
and a recombination rate per basepair) of 1.3 × 10−8
and simulate 3,000 1 Mb regions for a total of 3GB of simulated sequence."

ms command in SM, p.16:

ms 202 1 -t 480 -r 519.99948 1000000 -I 4 100 100 1 1 -en 0 3 0.1 -en 0 4 0.1
-es 0.04310345 1 0.97 -ej 0.05387931 5 3 -ej 0.05387931 2 1
-es 0.055 1 ${admix} -es 0.2155172 4 0.94 -en 0.2155172 7 0.2
-es 0.2836207 3 0.95 -ej 0.2836207 8 1 -ej 0.3577586 4 3
-ej 0.4741379 3 1 -ej 0.6 6 1 -ej 0.7844828 7 1

Given the mutation rate and L=1e6 (cf -r), N0 = 10000 (= theta/(4*mu*L))

###

Generation time (yrs.) 29

PopMap = ["CEU", "YRI", "Neanderthal", "Denisovan", "Anc1", "Anc2", "Anc3"]

"""

import msprime, demes

def samples(ARGS):
    Samples, PopLabels = [], []
    for x in ARGS["samples"]:
        x = x.split(":")
        assert len(x)==4, "Sample definition should have 4 `:`-delimited values, DemeName:nDiploids:TimeYa:Label"
        assert all(c.isalpha for c in x[0]), "For sample definition, you must provide the name of the deme to sample in."
        Samples.append( msprime.SampleSet(  num_samples = int(x[1]),
                                            population = x[0],
                                            time = float(x[2]) / ARGS["generation_time"],
                                            ploidy = 2) )
        PopLabels += [x[3]] * int(x[1]) * 2
    return Samples, PopLabels

def history(Args, Pars):
    N0 = 10000 # = 480 / (4*1.2e-8*1e6)
    alpha = 0.05 # cf. p.9
    #
    cmd = "-I 4 100 100 1 1 "
    cmd += "-en 0 3 0.1 "
    cmd += "-en 0 4 0.1 "
    cmd += "-es 0.04310345 1 0.97 "
    cmd += "-ej 0.05387931 5 3 "
    cmd += "-ej 0.05387931 2 1 "
    cmd += "-es 0.055 1 "+str(alpha)+" "
    cmd += "-es 0.2155172 4 0.94 "
    cmd += "-en 0.2155172 7 0.2 "
    cmd += "-es 0.2836207 3 0.95 "
    cmd += "-ej 0.2836207 8 1 "
    cmd += "-ej 0.3577586 4 3 "
    cmd += "-ej 0.4741379 3 1 "
    cmd += "-ej 0.6 6 1 "
    cmd += "-ej 0.7844828 7 1 "
    #
    Demes = ["CEU", "YRI", "Neanderthal", "Denisovan", "Anc1", "Anc2", "UA"]
    #
    Graph = demes.from_ms(cmd, N0 = N0, deme_names = Demes)
    Demography = msprime.Demography.from_demes(Graph)


    Args["generation_time"] = 29
    #Args["mut_rate"] = Pars["mutation_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels





#___
