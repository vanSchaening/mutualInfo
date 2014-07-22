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
                  help="matrix of correlations between single miRNAs and mRNAs.")
parser.add_option("--outdir", dest="outdir",
                  help="output directory")

(files,args) = parser.parse_args()

if not files.outdir.endswith("/"):
    files.outdir = files.outdir+"/"
# ------------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pprint import pprint
import cPickle as pickle
from itertools import izip
import numpy
 
max_p = 0.05
fname = files.MIdata.split("/")[-1]
outfile = files.outdir + fname.strip(".txt") + "_with_corr.txt"

# Read data from MIdata and store it in a dict
# {(mir,tf,t):(MI, correlation(mir,t))}
MI = dict()
f = open(files.MIdata,'r')
for line in f:
    # Get line and MI value, if it exists
    line = line.strip().split()
    try:
        val = float(line[3])
    except ValueError:
        continue
    if not line[0] in MI:
        MI[line[0]] = dict()
    # Store value
    MI[line[0]][tuple(line[1:3])] = (val,0)
f.close()

# Get names of miRs in correlation file and the id's of the ones that 
# appear in the results
f = open(files.corrfile, 'r')
names = f.readline().strip().split()
names = [ name.strip('"').lstrip('"') for name in names ]
names = [ name.lower() for name in names ]
ids = [ names.index(mir) for mir in MI.keys() ]

# Store correlation data in the dictionary
for line in f:
    line = line.strip().split()
    mrna = line[0].strip('"').lstrip('"')
    for i in ids:
        mir = names[i]
        for (tf, target) in MI[mir]:
            if mrna == target:
                MI[mir][(tf,target)] = (MI[mir][(tf,target)][0],float(line[i+1]))

# -------- XY-SCATTER - MI vs correlation --------------------------------------

# Make a list of MI and a list of correlations to pass to scatter()
# At the same time, store information in outfile
MI_new = [] 
pairs = []
o = open(outfile, 'w')
for mir, interactions in MI.iteritems():
    for (tf, t), (mi, corr) in interactions.iteritems():
        MI_new.append((mi,corr))
        pairs.append((mir,tf,t))
        o.write("\t".join([mir,tf,t,
                           str(mi),str(corr)])+"\n")
o.close()
MI = zip(*MI_new)
MI_new = []

plt.scatter(MI[0],MI[1])
plt.xlabel("Change in MI between TF and target, with respect to miR.")
plt.ylabel("Correlation between miR and target.")
plt.savefig(files.outdir + fname.strip(".txt") + "_vs_correlation.png")

exit()

# -------- PROBABILITY DISTRIBUTION --------------------------------------------

from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D

# Use gaussian kde to approximate the probability distribution
distribution = gaussian_kde([MI[0],MI[1]])
P = [ distribution.evaluate([m,c]) 
      for (m, c) in izip(MI[0],MI[1]) ]

# Plot points in a 3d plot, where the vertical value is the probability
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(MI[0], MI[1], P)
plt.savefig(files.outdir+fname.strip(".txt")+"_vs_correlation_dist.png")



