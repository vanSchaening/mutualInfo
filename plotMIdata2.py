# First attempt at visualizing data
# Receives output from DoubleInteractions or IndirectInteractions
# Generates a histogram of deltaMI values and
# Plots the deltaMI for every (mir,TF,t) against correlation(mir, t)

# -------- Argument Parsing ----------------------------------------------------
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--results", dest="MIdata",
                  help="output from MI scripts.")
parser.add_option("--correlation", dest="corrfile",
                  help="matrix of correlations between mRNAs in a tissue.")
parser.add_option("--outdir", dest="outdir",
                  help="output directory")

(files,args) = parser.parse_args()

if not files.outdir.endswith("/"):
    files.outdir = files.outdir+"/"
# ------------------------------------------------------------------------------

def getName(word):
    return word.split("|")[0]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pprint import pprint
import cPickle as pickle
from itertools import izip
import numpy
 
max_p = 0.05
fname = files.MIdata.split("/")[-1]
outfile = files.outdir + fname.strip(".txt") + "_RNA_corr.txt"

# Read data from MIdata and store it in a dict
# {(mir,tf,t):(MI, correlation(tf,t))}
# { mir: {(tf,t):(MI,corr)} }
MI = dict()
f = open(files.MIdata,'r')
for line in f:
    line = line.strip().split()
    try:
        val = float(line[3])
    except ValueError:
        continue
    interaction = tuple(line[0:3])
    if not interaction in MI:
        MI[interaction] = (val,0)
f.close()

# Open files
f = open(files.corrfile, 'r')
o = open(outfile, 'w')

# Get names of mRNAs in correlation file and a list of transcription factors
names = f.readline().strip().split()
tfs = set([ tf for (mir,tf,t) in MI ])

# Find tf-t correlations for each interaction, write results to outfile
for line in f:
    line = line.strip().split()
    if not getName(line[0]) in tfs:
        continue
    for (mir,tf,t) in MI:
        if tf == getName(line[0]):
            mutualinfo = MI[(mir,tf,t)][0]
            corr = line[names.index(t)+1]
            if corr == "nan" or corr == "-2":
                corr = -2
            MI[(mir,tf,t)] = (mutualinfo,float(corr))
            o.write("\t".join([mir,tf,t,str(mutualinfo),corr])+"\n")
o.close()

# -------- XY-SCATTER - MI vs correlation --------------------------------------

# Make a list of MI and a list of correlations to pass to scatter()
# At the same time, store information in outfile
MI = [ values for interaction, values in MI.iteritems()]                  
MI = zip(*MI)


plt.scatter(MI[0],MI[1])
plt.xlabel("Change in MI between TF and target, with respect to miR.")
plt.ylabel("Correlation between TF and target.")
plt.savefig(files.outdir + fname.strip(".txt") + "_vs_rna_correlation.png")

exit()



