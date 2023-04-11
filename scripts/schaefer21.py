#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> scrm command
Schaefer et al. (2021)
Base of the model is derived from Gutenkunst et al. (2019) OOA model.

scrm command in main text, p.12:

scrm 456 1 -t 17253.7128713 -r 13802.970297 25000000 -T -I 7 150 150 150 0 0 0 0 -eI q 0 0 0 0 0 2 0 -eI 0.0891112545728 0 0 0 2 0 0 0 -eI 0.0579585395596 0 0 0 0 0 0 2 -n 1 1.68 -n 2 3.74 -n 3 7.29 -n 4 0.231834158238 -n 5 0.231834158238 -n 6 0.231834158238 -n 7 0.0260813428018 -eg 0 2 116.010723 -eg 0 3 160.246047 -m 2 3 2.797460 -m 3 2 2.797460 -ej 0.028985 3 2 -en 0.028985 2 0.287184 -em 0.028985 1 2 7.293140 -em 0.028985 2 1 7.293140 -es 0.0144896348899 3 [1-x] -ej 0.0144896348899 8 7 -es 0.0362240872247 2 [1-x] -ej 0.0362240872247 9 5 -ej 0.0724481744495 6 5 -ej 0.197963 2 1 -en 0.303501 1 1 -ej 0.099616239868 5 4 -ej 0.304282332688 7 4 -ej 0.416577003084 4 1, where [1-x] means one minus the admixture proportion in a given run.

Mutation rate: 0.5e-9 per year => 1.25e-8
Recombination rate: 1cM/Mb
Given the mutation rate (1.5e-8) and L=25e6 (cf -r), N0 = 13803 (= theta/(4*mu*L))

###

Generation time (yrs.) 25


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
    N0 = 13803 # = 17253.7128713 / (4*1.25e-8*25e6)
    # removed these scrm-specific switches of ancient sampling:
    # -eI q 0 0 0 0 0 2 0 -eI 0.0891112545728 0 0 0 2 0 0 0 -eI 0.0579585395596 0 0 0 0 0 0 2
    alpha = 0.02
    cmd = "-I 7 150 150 150 0 0 0 0 -n 1 1.68 -n 2 3.74 -n 3 7.29 -n 4 0.231834158238 -n 5 0.231834158238 \
    -n 6 0.231834158238 -n 7 0.0260813428018 -eg 0 2 116.010723 -eg 0 3 160.246047 -m 2 3 2.797460 \
    -m 3 2 2.797460 -ej 0.028985 3 2 -en 0.028985 2 0.287184 -em 0.028985 1 2 7.293140 \
    -em 0.028985 2 1 7.293140 -es 0.0144896348899 3 "+str(1-alpha)+" -ej 0.0144896348899 8 7 \
    -es 0.0362240872247 2 "+str(1-alpha)+" -ej 0.0362240872247 9 5 \
    -ej 0.0724481744495 6 5 -ej 0.197963 2 1 -en 0.303501 1 1 -ej 0.099616239868 5 4 \
    -ej 0.304282332688 7 4 -ej 0.416577003084 4 1"
    #
    Demes = ["YRI", "CEU", "CHB", "Altai", "Nea_Intro", "Vindija", "Den"]
    #
    Graph = demes.from_ms(cmd, N0 = N0, deme_names = Demes)
    Demography = msprime.Demography.from_demes(Graph)

    Args["generation_time"] = 25
    #Args["mut_rate"] = Pars["mutation_rate"]
    #Args["mut_rate"] = Pars["recombination_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels





#___
