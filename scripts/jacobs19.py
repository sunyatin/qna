#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> stdpopsim model
Out-of-Africa with archaic admixture into Papuans
Jacobs et al., 2019. https://doi.org/10.1016/j.cell.2019.02.035
Malaspinas et al., 2016. https://doi.org/10.1038/nature18299
Generation time (yrs.) 29

PopMap = [{"pop_0": "YRI",
          "pop_1": "CEU",
          "pop_2": "CHB",
          "pop_3": "Papuan",
          "pop_4": "DenA",
          "pop_5": "NeaA",
          "pop_6": "Den1",
          "pop_7": "Den2",
          "pop_8": "Nea1",
          "pop_9": "Ghost" }]

"""

import stdpopsim, demes, msprime

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
    model = stdpopsim.get_species("HomSap").get_demographic_model("PapuansOutOfAfrica_10J19")
    Demography = model.model

    Args["generation_time"] = model.generation_time
    #Args["mut_rate"] = Pars["mutation_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels

#___
