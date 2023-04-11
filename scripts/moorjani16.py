#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> ms command
Simulation (a) Estimation of the date of Neanderthal admixture
Moorjani et al. (2016)

ms command in SM, p.7:

ms 44 1 -r 20000 50000000 -t 30000 -I 6 20 20 1 1 1 1 -en 0 1 1 -en 0 2 1 -en 0 3 1e-10 \
-en 0 4 1e-10 -en 0 5 1e-10 -en 0 6 1e-10 -es t n 2 0.97 -en 0.02500025 7 0.25 -en \
0.02500025 2 1 -ej 0.05 4 3 -ej 0.05 6 5 -en 0.05000025 3 0.25 -en 0.05000025 5 0.25 \
-ej 0.0500025 5 3 -en 0.050005 3 0.25 -ej 0.075 2 1 -en 0.0750025 1 1 -ej 0.1 7 3 -en \
0.1000025 3 0.25 -ej 0.3 3 1 -en 0.3000025 1 1


###

Generation time (yrs.) 28 [estimated from empirical ancient DNA data]

pop1 = west Africans, pop2 = Europeans, pop3-6 = Neanderthals (to simulate branch shortening in ms, we generate
four haploid chromosomes to obtain one diploid ancient genome)

"""

import msprime, demes, numpy

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
    # originally, for model S2a, should be `N0 = 10000`
    N0 = 14000 # contrary to model S2a (10000) used the No from S2b because better estimates with the latter at mu=1.2e-8
    t_n = 2000 # Manuscript p.2: "we set the date of the shared Neanderthal gene flow to 2,000 generations ago (9)"
    cmd = "-I 6 20 20 1 1 1 1 -en 0 1 1 -en 0 2 1 -en 0 3 1e-10 \
-en 0 4 1e-10 -en 0 5 1e-10 -en 0 6 1e-10 -en 0.02500025 2 1 \
-ej 0.05 4 3 -ej 0.05 6 5 -en 0.05000025 3 0.25 -en 0.05000025 5 0.25 \
-ej 0.0500025 5 3 -en 0.050005 3 0.25 -ej 0.075 2 1 -en 0.0750025 1 1 -en 0.1000025 3 0.25 \
-ej 0.3 3 1 -en 0.3000025 1 1"
    cmd += " -es "+str(numpy.true_divide(t_n, 4*N0))+" 2 0.97"
    #cmd += " -en 0.02500025 7 0.25" this command has to be changed to avoid bugs (Note: change has no impact on the coalescent)
    cmd += " -en "+str(numpy.true_divide(t_n, 4*N0))+" 7 0.25"
    cmd += " -ej 0.1 7 3"
    #
    Demes = ["YRI", "CEU", "Nea_1", "Nea_2", "Nea_3", "Nea_4", "Nea_intro"]
    #
    print(Demes)
    Graph = demes.from_ms(cmd, N0 = N0, deme_names = Demes)
    print(Graph)
    Demography = msprime.Demography.from_demes(Graph)

    Args["generation_time"] = 28
    #Args["mut_rate"] = Pars["mutation_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels





#___
