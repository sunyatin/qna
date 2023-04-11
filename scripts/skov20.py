#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 21:31:58 2022
@author: rtournebize

==> msprime0-formatted demography
Iceland - Archaic (Den + Nea) admixture into ancestral Europeans
Skov et al. (2020)

Demography downloaded from:
    https://github.com/LauritsSkov/ArchaicSimulations/blob/master/SI%203%20Dataset%20-%20Simulation%20script.py

###

Generation time (yrs.) 29

PopMap = ["Africa", "Iceland", "Den_intro", "Den_seq", "Altai", "Nea_intro", "Vindija"]

"""

import msprime as msp
import random

def samples(ARGS):
    Samples, PopLabels = [], []
    for x in ARGS["samples"]:
        x = x.split(":")
        assert len(x)==4, "Sample definition should have 4 `:`-delimited values, DemeName:nDiploids:TimeYa:Label"
        assert all(c.isalpha for c in x[0]), "For sample definition, you must provide the name of the deme to sample in."
        Samples.append( msp.SampleSet(  num_samples = int(x[1]),
                                            population = x[0],
                                            time = float(x[2]) / ARGS["generation_time"],
                                            ploidy = 2) )
        PopLabels += [x[3]] * int(x[1]) * 2
    return Samples, PopLabels

def history(Args, Pars):
    Demography = skov20_original()

    Args["generation_time"] = 29
    #Args["mut_rate"] = Pars["mutation_rate"]

    print("/!\ Warning. Replacing generation_time by the model-specific one: "+str(Args["generation_time"]))
    #print("/!\ Warning. Replacing mut_rate by the one suggested in the parameter file: "+str(Args["mut_rate"])+"\n")

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels

#___


def skov20_original():
    # Generation time, mutation rate and recomination rate
    gen_time = 29.0
    #rec_rate = 1.2e-8
    #mu = 1.2e-8
    #chrom = 19
    #Overall_scenario = 'test'

    # Populations sizes (number of individuals)
    Ne_Africa = 27122
    Ne_Europe = 5000
    Modernhumans = 7000
    Denisovasize = 5000
    Neanderthal_recentsize = 1000
    Neanderthal_latersize = 2000
    Neanderthal_Denisovasize = 5000

    # Split times (years)
    OutOfafrica = 62000/gen_time
    Denisova_split = 350000
    IntrogressingVindijaSplit = 90000
    DenisovaNeanderthal = 420000

    # Admix parameters (years and admixture proportions (in percent))
    Denisova_admix_time = 45000
    DenisovaProportion = 0.08
    admix_into = 1 # (1 for humans, 5 for neanderthals)

    # Number of samples
    #n_ingroup = 2
    #African_samples = 1000
    #samples = [msp.Sample(0, 0)]*African_samples + [msp.Sample(1, 0)]*n_ingroup + [msp.Sample(3, 80000/gen_time)]*n_ingroup + [msp.Sample(4, 120000/gen_time)]*n_ingroup + [msp.Sample(6, 60000/gen_time)]*n_ingroup

    population_configurations = [
    msp.PopulationConfiguration(initial_size = Ne_Africa, metadata = {"name": "Africa"}), #0
    msp.PopulationConfiguration(initial_size = Ne_Europe, metadata = {"name": "Iceland"}), #1
    msp.PopulationConfiguration(initial_size = Denisovasize, metadata = {"name": "Den_intro"}), #2
    msp.PopulationConfiguration(initial_size = Denisovasize, metadata = {"name": "Den_seq"}), #3
    msp.PopulationConfiguration(initial_size = Neanderthal_recentsize, metadata = {"name": "Altai"}), #4
    msp.PopulationConfiguration(initial_size = Neanderthal_recentsize, metadata = {"name": "Nea_intro"}), #5
    msp.PopulationConfiguration(initial_size = Neanderthal_recentsize, metadata = {"name": "Vindija"}), #6
    ]

    demographic_events_dict = {

    # admixture times
    1551.06896552 + random.uniform(0,1)/100000.0: msp.MassMigration(time = 1551.06896552, source = 1, destination = 5,proportion = 0.02),
    Denisova_admix_time/gen_time + random.uniform(0,1)/100000.0: msp.MassMigration(time = Denisova_admix_time/gen_time, source = admix_into, destination = 2 ,proportion = DenisovaProportion),

    # Human parameters
    OutOfafrica - 300 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica - 300, initial_size = 1305, growth_rate = 0, population_id = 1),
    OutOfafrica - 200 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica - 200, initial_size = 5000, growth_rate = 0, population_id = 1),
    OutOfafrica - 100 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica - 100, initial_size = 250, growth_rate = 0, population_id = 1),
    OutOfafrica + random.uniform(0,1)/100000.0: msp.MassMigration(time = OutOfafrica, source = 1, destination = 0, proportion = 1.0),
    OutOfafrica + 0.0001 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica  + 0.0001, initial_size = Modernhumans, growth_rate = 0, population_id = 0),

    # archaic population merges parameters
    IntrogressingVindijaSplit/gen_time + random.uniform(0,1)/100000.0: msp.MassMigration(time = IntrogressingVindijaSplit/gen_time, source = 6, destination = 5, proportion = 1.0),
    130000/gen_time - 0.0001 + random.uniform(0,1)/100000.0: msp.MassMigration(time = 130000/gen_time - 0.0001, source = 5, destination = 4, proportion = 1.0),
    Denisova_split/gen_time - 0.0001 + random.uniform(0,1)/100000.0: msp.MassMigration(time = Denisova_split/gen_time, source = 3, destination = 2, proportion = 1.0),
    DenisovaNeanderthal/gen_time - 0.0001 + random.uniform(0,1)/100000.0: msp.MassMigration(time = DenisovaNeanderthal/gen_time, source = 4, destination = 2, proportion = 1.0),

    # psms pop sizes
    130000/gen_time + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = 130000/gen_time + 0.0001, initial_size = Neanderthal_latersize, growth_rate = 0, population_id = 4),

    # Neanderthal and denisova merges and have a bigger pop size
    DenisovaNeanderthal/gen_time + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = DenisovaNeanderthal/gen_time, initial_size = Neanderthal_Denisovasize, growth_rate = 0, population_id = 2),

    # Humans and archaics are the same population
    580000/gen_time + random.uniform(0,1)/100000.0: msp.MassMigration(time = 580000/gen_time, source = 2, destination = 0,proportion = 1.0),
    580000/gen_time + 0.0001 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = 580000/gen_time + 0.0001, initial_size = 7000, growth_rate = 0, population_id = 0),

    }


    demographic_events = []
    for i, key in enumerate(sorted(demographic_events_dict)):
        demographic_events.append(demographic_events_dict[key])
    Demography = msp.Demography.from_old_style(
        population_configurations,
        demographic_events = demographic_events,
        population_map = None,
    )
    return Demography
