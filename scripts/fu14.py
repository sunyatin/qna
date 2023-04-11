#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> ms command
Simulation 4
Fu et al. (2014)
Note that this command is the same as the one in Moorjani et al. (2016) SM p. 7-8, Simulation (b) Complex demographic scenarios

ms command in SM, p.113:

ms 68 1 -I 11 20 20 20 1 1 1 1 1 1 1 1 -en 0 1 1 -en 0 2 2.41428571428571 -en 0 3 3.23571428571429 -eg 0 2
97.6909397920288 -eg 0 3 123.512032930849 -en 0 4 7.14285714285714e-11 -en 0 5 7.14285714285714e-11 -
en 0 6 7.14285714285714e-11 -en 0 7 7.14285714285714e-11 -en 0 8 7.14285714285714e-11 -en 0 9
7.14285714285714e-11 -en 0 10 7.14285714285714e-11 -en 0 11 7.14285714285714e-11 -ej
0.0321428571428571 5 4 -ej 0.0321428571428571 7 6 -en 0.0321430357142857 4 0.714285714285714 -en
0.0321430357142857 6 0.714285714285714 -ej 0.0321446428571429 6 4 -en 0.0321464285714286 4
0.714285714285714 -ej 0.0357142857142857 2 3 -ej 0.0357142857142857 4 3 -en 0.0357160714285714 3
0.132857142857143 -en 0.0357160714285714 3 0.132857142857143 -es 0.0392857142857143 3 0.97 -en
0.0392875 12 0.178571428571429 -en 0.0392875 3 0.132857142857143 -ej 0.0428571428571429 9 8 -ej
0.0428571428571429 11 10 -en 0.0428573214285714 8 0.178571428571429 -en 0.0428573214285714 10
0.178571428571429 -ej 0.0428589285714286 10 8 -en 0.0428607142857143 8 0.178571428571429 -ej
0.0535714285714286 3 1 -en 0.0535732142857143 1 1 -ej 0.0714285714285714 12 8 -en 0.0714303571428571
8 0.178571428571429 -en 0.107142857142857 1 0.521428571428571 -ej 0.214285714285714 8 1 -en
0.2142875 1 0.714285714285714 -r 28000 50000000 -t 42000 -p 12 -seeds 46 47 48

Given the mutation rate (1.5e-8) and L=50e6 (cf -r), N0 = 14000 (= theta/(4*mu*L))

###

Generation time (yrs.) 29

pop1= African, pop2=present-day Europeans, pop3=present-day Asian, pop4-7= Ust’-Ishim, pop8-11= Neandertal.

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
    N0 = 14000 # = 42000 / (4*1.5e-8*50e6)
    # The admixture will be handled by msprime command after the ms=>msprime conversion
    # Thus commenting this switch:
    # -es 0.0392857142857143 3 0.97
    # -en 0.0392875 12 0.178571428571429
    # -ej 0.0714285714285714 12 8
    cmd = "-I 11 20 20 20 1 1 1 1 1 1 1 1 -en 0 1 1 -en 0 2 2.41428571428571 -en 0 3 3.23571428571429 -eg 0 2 \
97.6909397920288 -eg 0 3 123.512032930849 -en 0 4 7.14285714285714e-11 -en 0 5 7.14285714285714e-11 \
-en 0 6 7.14285714285714e-11 -en 0 7 7.14285714285714e-11 -en 0 8 7.14285714285714e-11 -en 0 9 \
7.14285714285714e-11 -en 0 10 7.14285714285714e-11 -en 0 11 7.14285714285714e-11 -ej \
0.0321428571428571 5 4 -ej 0.0321428571428571 7 6 -en 0.0321430357142857 4 0.714285714285714 -en \
0.0321430357142857 6 0.714285714285714 -ej 0.0321446428571429 6 4 -en 0.0321464285714286 4 \
0.714285714285714 -ej 0.0357142857142857 2 3 -ej 0.0357142857142857 4 3 -en 0.0357160714285714 3 \
0.132857142857143 -en 0.0357160714285714 3 0.132857142857143 \
-en 0.0392875 3 0.132857142857143 -ej 0.0428571428571429 9 8 -ej \
0.0428571428571429 11 10 -en 0.0428573214285714 8 0.178571428571429 -en 0.0428573214285714 10 \
0.178571428571429 -ej 0.0428589285714286 10 8 -en 0.0428607142857143 8 0.178571428571429 -ej \
0.0535714285714286 3 1 -en 0.0535732142857143 1 1 -en 0.0714303571428571 \
8 0.178571428571429 -en 0.107142857142857 1 0.521428571428571 -ej 0.214285714285714 8 1 -en \
0.2142875 1 0.714285714285714"
    #
    Demes = ["YRI", "CEU", "CHB", "UI_1", "UI_2", "UI_3", "UI_4", "Nea_1", "Nea_2", "Nea_3", "Nea_4"]
    #
    Graph = demes.from_ms(cmd, N0 = N0, deme_names = Demes)
    Demography = msprime.Demography.from_demes(Graph)
    # add the Nea admixture edge
    Demography.add_population(name="Nea_Intro", initial_size=0.178571428571429*N0)
    Demography.add_mass_migration(time = 0.0392857142857143*4*N0, source = "CHB", dest = "Nea_Intro", proportion = 0.03)
    Demography.add_mass_migration(time = 0.0714285714285714*4*N0, source = "Nea_Intro", dest = "Nea_1", proportion = 1.0)
    Demography.sort_events()


    Args["generation_time"] = 29
    #Args["mut_rate"] = Pars["mutation_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels





#___
