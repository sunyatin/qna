#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> ms command
Bottleneck older than admixture (best fit to DCFS), with admixture rate of 0.05 and 4Nm=0 (best fit from Fig2c)
Yang et al. (2012)

ms command in Appendix B, p. 8 using parameters of Table 1

###

Generation time (yrs.) 25

PopMap = ["NEA", "YRI", "CEU"]

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
    N0 = 10000 # p.4
    # -t 20
    cmd = "-I 3 1 1 10 "
    cmd += "-n 2 100 "
    cmd += "-n 3 100 "
    cmd += "-m 3 2 0 "
    cmd += "-m 2 3 0 "
    cmd += "-es 0.05 3 0.95 "
    cmd += "-ej 0.05 4 1 "
    cmd += "-en 0.1 3 1 "
    cmd += "-ej 0.1125 3 2 "
    cmd += "-en 0.1150 2 1 "
    #cmd += "-en 0.125 3 100 "
    cmd += "-ej 0.3 2 1 "
    #
    Demes = ["NEA", "YRI", "CEU"]
    #
    Graph = demes.from_ms(cmd, N0 = N0, deme_names = Demes)
    Demography = msprime.Demography.from_demes(Graph)

    Args["generation_time"] = 25
    #Args["mut_rate"] = Pars["mutation_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels





#___
