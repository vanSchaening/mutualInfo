# Import functions. Argument parsing occurs at import.
# (files, args) = parser.parse_args()
from MI import *


# -------- Function Definitions ------------------------------------------------

# Assumes that the mirdict has already been trimmed
def getIndirectInteractions(mirdict,tfdict):
    indirect = dict()
    for mir, rnas in mirdict.iteritems():
        mir = mir.lower()
        indirect[mir] = dict()
        for rna, targets in tfdict.iteritems():
            indirect[mir][rna] = targets
    return indirect

def flipDictionary(interactions):
    flipped = dict()
    for mir, mrnas in interactions.iteritems():
        for mrna, targets in mrnas.iteritems():
            if not targets:
                continue
            try:
                flipped[mrna].extend([ (mir,target) 
                                       for target in targets ])
            except KeyError:
                flipped[mrna] = [ (mir,target)
                                  for target in targets ]
    return flipped

def getExpressionData(mirfile, rnafile, interactions):
    # Set file pointers to the beginning of the data
    for f in [mirfile,rnafile]:
        f.seek(0)
        f.readline()
    # Initialize dictionary
    expression = dict()
    for tf, pairs in interactions.iteritems():
        expression[tf] = list()
        for (mir, target) in pairs:
            expression[mir], expression[target] = list(), list()
    # Get RNA expression data
    for line in rnafile:
        line = line.strip().split()
        if proteinName(line[0]) in expression:
            expression[proteinName(line[0])] = map(float,line[1:])
    # Get miR expression data
    for line in mirfile:
        line = line.strip().split()
        if line[0] in expression:
            expression[line[0]] = map(float,line[1:])
    return expression

# -------- MAIN / TESTING ------------------------------------------------------

from pprint import pprint

# Get list of miRNAs
mirlist = getMirList(files.mirlist)

# Unpack dictionaries containing {tf:[target]} and {mir:[target]} interactions
tfdict = unpickle(files.rnapkl)
mirdict = trimInteractions(mirlist,
                           unpickle(files.mirpkl))

# Get indirect interactions {mir:{target:[targets]}}, and write to a pkl
indirect_interactions = getIndirectInteractions(mirdict,tfdict)
d = open(files.outdir + "indirect_interactions.pkl",'w')
pickle.dump(indirect_interactions, d)
d.close()
# Flip dictionary into {tf:[(mir,target)]}
indirect_interactions = flipDictionary(indirect_interactions)

# miRNA and mRNA expression matrices, filtered for shared samples
mirfile, rnafile = filterExpressionFiles(files.mirdata, files.rnadata)
mirfile = open(mirfile, 'r')
rnafile= open(rnafile, 'r')

# Get sample names from both expression matrices
mir_samples = getSampleNames(mirfile)
rna_samples = getSampleNames(rnafile)

# Map the sample columns in the mir expression file to those in the
# rna expression file
id_map = mapCorrespondingIndices(mir_samples,rna_samples)

expression = getExpressionData(mirfile,rnafile,indirect_interactions)
expression = dict( (target, exp) 
                   for target, exp in expression.iteritems()
                   if exp )

for tf, interactions in indirect_interactions.iteritems():
    # Find expression outliers and their IDs
    try:
        tfexp = expression[tf]
    except KeyError:
        print "Missing mRNA expression data for", tf, "\n"
        continue
    if not tfexp:
        print "Missing mRNA expression data for", tf, "\n"
        continue
    top, bot = selectOutliers(tfexp)
    top_id, bot_id = getIndices(top, tfexp), getIndices(bot, tfexp)

