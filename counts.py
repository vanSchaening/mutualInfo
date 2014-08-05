# Receive MI results and count how many times each component
# appears as a tf or a target, and how many times each miR appears
import sys
if len(sys.argv) != 3:
    print ("Usage:\n\tcounts.py <infile> <outfile>\n")

from collections import defaultdict
import operator

# Open the input file and an output file for each role
f = open(sys.argv[1], 'r')
filemap = { 'miR' : open(sys.argv[2]+"_sortedMirs.txt",'w') ,
            'TF' : open(sys.argv[2]+"_sortedTFs.txt",'w') ,
            'target' : open(sys.argv[2]+"_sortedTargets.txt",'w') }

# Skip the header
f.readline()

# Initialize a dictionary of counts
counts = { 'miR' : defaultdict(int),
           'TF' : defaultdict(int),
           'target' : defaultdict(int) }

# Over all interactions, count the number of times a component
# appears in a given role
for line in f:
    line = line.strip().split()
    counts['miR'][line[0]] += 1
    counts['TF'][line[1]] += 1
    counts['target'][line[2]] += 1

# Sort each list by count and write them to file
def writeCountsToFile(counts,f):
    for (name, count) in counts:
        f.write(name + "\t" + str(count) + "\n")
for name, numbers in counts.iteritems():
    counts[name] = sorted(numbers.iteritems(), key=operator.itemgetter(1), reverse=True)
    writeCountsToFile(counts[name],filemap[name])

# Close files
for o in filemap.values():
    o.close()
