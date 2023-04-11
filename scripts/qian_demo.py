#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 00:25:57 2021
@author: rtournebize

Thrown on Dec 8th 2021, 00:27.
Being humble and yielding. Persisting makes things go well.

Encoding of the parameters:

    *** All time values are in years BP ***

    - n.EPOCH.CHAIN: size change

    - m.EPOCH.CHAIN<: within-chain migration rate changes from right to left deme (forward)
    - m.EPOCH.CHAINA_from_CHAINB: between-chain migration rate change, with p inds in CHAINA coming from (backward) CHAINB

    - j.EPOCH.CHAINA_into_CHAINB: expansion starting (forward) at EPOCH from CHAINA into CHAINB, the value corresponds to the time elapsed between each successive founder event
    - jb.EPOCH.CHAINA_into_CHAINB: adds an instantaneous bottleneck just after the founding, VALUE is the strength of the bottleneck (cf msp1 doc)
    - jn.EPOCH.CHAINA_into_CHAINB: adds a single-generation bottleneck just after the founding, VALUE is the size of the bottleneck (= diploid Nf)
    - jd.EPOCH.CHAINA_into_CHAINB: sets the duration of the bottleneck for jn.EPOCH.CHAINA_into_CHAINB (if absent, assumes that the duration is 1 generation)

    - T.EPOCH: defines the timing in years BP of EPOCH

    - C.CHAIN: defines the chain CHAIN, the value corresponds to the number of demes

    - a.EPOCH.CHAIN__WITHININDEXA_from_CHAIN__WITHININDEXB: admixture event at EPOCH, with VALUE lineages moving from CHAIN__WITHININDEXA into CHAIN_WITHININDEXB (backward)

    - g.EPOCH.CHAIN: [Works only for growth till present!]: exponential growth, the value is the growth rate: g = -1/T * log(No/Nf)

"""

import msprime as msp
import copy

def samples(Demo, ARGS, PAR):
    Samples, PopLabels = [], []
    for x in ARGS["samples_within"]:
        x = x.split(":")
        assert len(x)==5, "Sample definition should have 5 `:`-delimited values, Chain:Index:nDiploids:TimeYa:Label"
        if any(c.isalpha() for c in x[1]):
            Pops = [y.name for y in Demo.populations if y.name.startswith(x[0]+"__")]
            x[1] = int(round(float(PAR[x[1]]) * (len(Pops)-1)))
        elif "." in x[1]:
            Pops = [y.name for y in Demo.populations if y.name.startswith(x[0]+"__")]
            x[1] = int(round(float(x[1]) * (len(Pops)-1)))
        Samples.append( msp.SampleSet(  num_samples = int(x[2]),
                                        population = x[0]+"__"+str(x[1]),
                                        time = float(x[3]) / ARGS["generation_time"],
                                        ploidy = 2) )
        PopLabels += [x[4]] * int(x[2]) * 2
    return Samples, PopLabels


def connector(Demo, source_chain, dest_chain):
    c0 = [x.id for x in Demo.populations if x.name.startswith(source_chain+"__")]
    c1 = [x.id for x in Demo.populations if x.name.startswith(dest_chain+"__")]
    if c0[0] < c1[0]:
        source, dest = c0[-1], c1[0]
    else:
        source, dest = c0[0], c1[-1]
    source, dest = Demo.populations[source].name, Demo.populations[dest].name
    return source, dest


def migration_within(Demo, PAR, epoch, chain):
    Pops = [x for x in Demo.populations if x.name.startswith(chain+"__")]
    for pop in Pops:
        chain, idw = pop.name.split("__")
        idw = int(idw)
        try:
            Demo.add_migration_rate_change(  time = PAR["T."+epoch] / PAR["generation_time"],
                                             source = chain+"__"+str(idw),
                                             dest = chain+"__"+str(idw+1),
                                             rate = PAR["m."+epoch+"."+chain+"<"])
        except: pass
        try:
            Demo.add_migration_rate_change(  time = PAR["T."+epoch] / PAR["generation_time"],
                                             source = chain+"__"+str(idw),
                                             dest = chain+"__"+str(idw-1),
                                             rate = PAR["m."+epoch+"."+chain+">"])
        except: pass


def migration_between(Demo, PAR, epoch, source, dest):
    rate = PAR["m."+epoch+"."+source+"_from_"+dest]
    source, dest = connector(Demo, source, dest)
    Demo.add_migration_rate_change(  time = PAR["T."+epoch] / PAR["generation_time"],
                                     source = source,
                                     dest = dest,
                                     rate = rate)


def resize_within(Demo, PAR, epoch, chain):
    Pops = [x.name for x in Demo.populations if x.name.startswith(chain+"__")]
    for pop in Pops:
        Demo.add_population_parameters_change( time = PAR["T."+epoch] / PAR["generation_time"],
                                               initial_size = PAR["n."+epoch+"."+chain], \
                                               growth_rate = 0., \
                                               population = pop)


def migration(Demo, PAR, x):
    _, epoch, x = x.split(".")
    if "_from_" in x:
        source, dest = x.split("_from_")
        migration_between(Demo, PAR, epoch, source, dest)
    else:
        chain = x.split("<")[0].split(">")[0]
        migration_within(Demo, PAR, epoch, chain)


def expansion(Demo, PAR, epoch, source, dest, lapse):
    # assumes that the expansion proceeds from left to right deme (forward in time)
    Pops = [x.name for x in Demo.populations if x.name.startswith(source+"__")]
    for i, pop in enumerate(Pops):
        if i == 0:
            _, ancestral = connector(Demo, source, dest)
        else:
            ancestral = Pops[i-1]
        ya = PAR["T."+epoch] - lapse * i
        Demo.add_population_split(time = ya / PAR["generation_time"],
                                  derived = [pop],
                                  ancestral = ancestral)
        assert not ("jb."+epoch+"."+source+"_into_"+dest in PAR.keys() and "jn."+epoch+"."+source+"_into_"+dest in PAR.keys()), "Cannot have concomitant jb and jn events"
        if "jb."+epoch+"."+source+"_into_"+dest in PAR.keys():
            Demo.add_instantaneous_bottleneck(time = ya / PAR["generation_time"] - 1.,
                                              population = pop,
                                              strength = PAR["jb."+epoch+"."+source+"_into_"+dest] / PAR["generation_time"])
        if "jn."+epoch+"."+source+"_into_"+dest in PAR.keys():
            if "jd."+epoch+"."+source+"_into_"+dest in PAR.keys():
                duration = PAR["jd."+epoch+"."+source+"_into_"+dest] / PAR["generation_time"]
            else:
                duration = 1.
            Demo.add_population_parameters_change( time = ya / PAR["generation_time"] - duration,
                                                   initial_size = PAR["jn."+epoch+"."+source+"_into_"+dest], \
                                                   growth_rate = 0., \
                                                   population = pop)


def history(ARGS, PAR):
    Keys = list(PAR.keys())
    Demo = msp.Demography()
    # Initialize
    PAR = copy.deepcopy(PAR)
    PAR["T.0"] = 0.
    PAR["generation_time"] = float(ARGS["generation_time"])
    # Metapopulation definition
    for key in [x for x in Keys if x.startswith("C.")]:
        chainName = key.replace("C.", "")
        nDemes = int(PAR[key])
        for i in range(nDemes):
            Demo.add_population(name = chainName+"__"+str(i),
                                default_sampling_time = 0.,
                                growth_rate = PAR["g."+"0."+chainName] if "g."+"0."+chainName in Keys else 0.,
                                initially_active = True,
                                initial_size = PAR["n."+"0."+chainName])
    # Demographic changes
    for key in Keys:
        if key.startswith("C."):
            pass
        elif key.startswith("m."):
            migration(Demo, PAR, key)
        elif key.startswith("n.") and not key.startswith("n.0."):
            _, epoch, chain = key.split(".")
            resize_within(Demo, PAR, epoch, chain)
        elif key.startswith("j."):
            _, epoch, x = key.split(".")
            source, dest = x.split("_into_")
            expansion(Demo, PAR, epoch, source, dest, PAR[key])
        elif key.startswith("a."):
            _, epoch, x = key.split(".")
            source, dest = x.split("_from_")
            Demo.add_mass_migration(time = PAR["T."+epoch] / PAR["generation_time"], source = source, dest = dest, proportion = PAR[key])
        else:
            pass
    # Sort and export
    Demo.sort_events()
    Samples, PopLabels = samples(Demo, ARGS, PAR)
    if False:
        print(Demo.debug())
    return ARGS, Demo, Samples, PopLabels


