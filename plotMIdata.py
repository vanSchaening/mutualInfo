import matplotlib.pyplot as plt
from pprint import pprint
import cPickle as pickle

# extract a dictionary from a .pkl
def unpickle(pklfile):
    f = open(pklfile, 'rb')
    pkl = pickle.load(f)
    f.close()
    return pkl

MIdata = "/home/axolotl/Results/MI_double_results.txt"
interactions = "/home/axolotl/Results/double_interaction.pkl"
corrfile = "/home/axolotl/Data/correlations/COAD_singlemiRNA_vs_mRNA_correlation_247samples.xls"

# Read data from MIdata and store it in a dict
# {(mir,tf,t):(MI,mir-t correlation)}
MI = dict()
f = open(MIdata,'r')
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
f = open(corrfile, 'r')
names = f.readline().strip().split()
names = [ name.strip('"').lstrip('"') for name in names ]
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

MI = zip(*[ (mi, correlation) for mir, interactions in MI.iteritems()
            for interaction,(mi,correlation) in interactions.iteritems() ])

plt.scatter(MI[0],MI[1])
plt.xlabel("Change in MI between TF and target, with respect to miR.")
plt.ylabel("Correlation between miR and target.")
plt.savefig("MI_vs_correlation.png")


exit()


# Plot a simple histogram to visualize the distribution
n, bins, patches = plt.hist(MI[0],35)
plt.savefig("MI_histogram.png")

