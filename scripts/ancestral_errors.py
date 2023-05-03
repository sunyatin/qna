#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 11:18:05 2021
@author: rtournebize

Implement an ancestral polarization error in *.geno*.gz files
   +
Copy the other *.ind, *.snp.gz and (if exists) *.par file.

0.0  231222
0.1  090123  implements the other-file copying to preserve the full dataset structure
0.2  090123  adds a `seed` option
0.3  030523  added option to argparse to print out default values when using --help

"""

import argparse, os, gzip, yaml, shutil
import numpy as np
from numpy.random import Generator, PCG64

__version__ = "0.3"

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', type=str, required=True, help="Prefix of the EIGENSTRAT dataset")
parser.add_argument('-o', '--output', type=str, required=True, help="Prefix of output the EIGENSTRAT dataset")
parser.add_argument('-p', '--proba_polarization_error', type=float, required=True, help="Probability of ancestral polarization error")
parser.add_argument('--ploidy', type=int, required=False, default=2, help="Ploidy")
parser.add_argument('--seed', type=int, required=False, default=None, help="Random seed")
ARGS = vars(parser.parse_args())

print(yaml.dump(ARGS))
p = ARGS["proba_polarization_error"]

assert p>=0 and p<=1, "Probability of polarization error must be contained between 0 and 1."
assert ARGS["input"] != ARGS["output"], "Input and output must have different prefixes."

def reverse_haploid(ARGS, suffix, Reverse):
    FOUT = gzip.open(ARGS["output"]+"."+suffix+".gz", "wt")
    with gzip.open(ARGS["input"]+"."+suffix+".gz", "rt") as F:
        i = 0
        for line in F:
            if Reverse[i]:
                line = np.array(list(line.strip())).astype(np.int8)
                line = 1 - line
                line = "".join(line.astype(str))+"\n"
            else:
                pass
            FOUT.write(line)
            i += 1
    FOUT.close()

if ARGS["seed"] is not None:
    rg0 = Generator(PCG64(ARGS["seed"]))

print("geno")
FOUT = gzip.open(ARGS["output"]+".geno.gz", "wt")
Reverse = []
with gzip.open(ARGS["input"]+".geno.gz", "rt") as F:
    for line in F:
        if ARGS["seed"] is not None:
            x = rg0.uniform()
        else:
            x = np.random.uniform()
        rev = False
        if x <= p:
            line = np.array(list(line.strip())).astype(np.int8)
            line = ARGS["ploidy"] - line
            line = "".join(line.astype(str))+"\n"
            rev = True
        else:
            pass
        Reverse.append(rev)
        FOUT.write(line)
FOUT.close()

Reverse = np.array(Reverse)
print("Reversed "+str(np.sum(Reverse))+" among "+str(len(Reverse))+"  |  "+str(np.round(100*np.true_divide(np.sum(Reverse), len(Reverse)), 2))+"%" )

print("geno1")
if os.path.exists(ARGS["input"]+".geno1.gz"):
    reverse_haploid(ARGS, "geno1", Reverse)

print("geno12")
if os.path.exists(ARGS["input"]+".geno12.gz"):
    reverse_haploid(ARGS, "geno12", Reverse)

for suffix in ["ind", "par", "snp", "snp.gz"]:
    if os.path.exists(ARGS["input"]+"."+suffix):
        shutil.copy2(ARGS["input"]+"."+suffix, ARGS["output"]+"."+suffix)

print("Done.")
