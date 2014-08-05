# Getting sorted files from various cancer types, see which miRs are repeatedly expressed
from scipy.stats import scoreatpercentile
import sys
from collections import defaultdict
from pprint import pprint

cutoff = 95

cancers = ([ cancer.strip("/") for cancer in sys.argv[1:] ])
files = dict( (cancer, "/".join([cancer,"sorted_interactions.txt"])) 
              for cancer in cancers )

mirs = dict( (cancer, defaultdict(int)) 
             for cancer in cancers )
for cancer in cancers:
    f = open(files[cancer])
    f.readline()
    score = [ float(line.strip().split()[3])
              for line in f ]
    min_MI = scoreatpercentile(score,85)
    f.seek(0)
    f.readline()
    for line in f:
        line = line.strip().split()
        if float(line[3]) > min_MI:
            mirs[cancer][line[0]] += 1
        else:
            break
    f.close()

pprint(mirs)
