# Perform hierarchical clustering on the set of interactions
from Pycluster import treecluster, read
import numpy as np
import sys

if len(sys.argv) < 3:
    exit("Usage:\n\tcluster.py <infile> <out_name>\nWhere <out_name> is a name for output files without an extension.\n")

# Open file into record object
f = open(sys.argv[1])
record = read(f)

# Apply a mask to missing values
def maskValue(v):
    if v == 0.0:
        return 0
    return 1
record.mask = np.array([ map(maskValue, row)
                         for row in record.data ])

# Hierarchical clustering
interactions_tree = record.treecluster(method='m',transpose=0)
cancer_tree = record.treecluster(method='m', transpose=1)

record.save(sys.argv[2], interactions_tree, cancer_tree)
f.close()



