# -------- Argument Parsing ----------------------------------------------------
# If you're reading through this code, skip this part. It's just storing file-
# names in the relevant variables.
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--mirlist",dest="mirlist",
                  help="plain text file with the name of a miRNA in each line.")
parser.add_option("--mirpkl",dest="mirpkl",
                  help="pickle with {miR:target} interactions.")
parser.add_option("--rnapkl",dest="rnapkl",
                  help="pickle with {TF:target} interactions for a specific tissue.")
parser.add_option("--mirdata",dest="mirdata",
                  help="expression matrix for individual miRNAs (not families) in a tissue. Obtained from TCGA data, with the samples labelled using TCGA barcodes.")
parser.add_option("--rnadata",dest="rnadata",
                  help="expression matrix for RNA transcripts in a specific tissue. In the same format as <mirdata>.")
parser.add_option("--outdir", dest="outdir",
                  help="Directory where result files will be written.")
(files, args) = parser.parse_args()


# -------- Imports -------------------------------------------------------------
from itertools import izip
import cPickle as pickle
from scipy.stats.kde import gaussian_kde
from scipy.stats import entropy, scoreatpercentile

# -------- Function Definitions ------------------------------------------------

# the two expression files (mir and rna) will rarely have exactly the same
# samples, so we filter both files to contain only the samples shared by both
def filterExpressionFiles(mirfile,rnafile):
    mir, rna = open(mirfile, 'r'), open(rnafile, 'r')
    mir_samples = getSampleNames(mir)
    rna_samples = getSampleNames(rna)
    overlap = set(mir_samples) & set(rna_samples)

    mir_indices = sorted([ mir_samples.index(name) for name in overlap ])
    rna_indices = sorted([ rna_samples.index(name) for name in overlap ])

    mir.seek(0)
    mir.readline()
    rna.seek(0)
    rna.readline()

    mir_file = mirfile.strip(".txt") + (".filtered.txt")
    rna_file = rnafile.strip(".txt") + (".filtered.txt")
    mir_out, rna_out = open(mir_file, 'w'), open(rna_file, 'w')
               
    mir_out.write("\t".join([mir_samples[i] for i in mir_indices])+"\n")
    rna_out.write("\t".join([rna_samples[i] for i in rna_indices])+"\n")

    for f,ids,o in izip([mir, rna],[mir_indices,rna_indices],[mir_out,rna_out]):
        for line in f:
            line = line.strip().split()
            out = [line[0]]
            out.extend([ line[i+1] for i in ids ])
            o.write("\t".join(out)+"\n")
        f.close()
        o.close()
    return mir_file, rna_file

# extract a dictionary from a .pkl
def unpickle(pklfile):
    f = open(pklfile, 'rb')
    pkl = pickle.load(f)
    return pkl

# retrieve list of miR or TF from file
def getMirList(mtffile):
    f = open(mtffile, 'r')
    mtflist = [ line.strip() for line in f ]
    f.close()
    return mtflist

# Get a list of all direct targets of a group of miRs or TFs
def getTargetsList(mtflist, interactions):
    mylist = list(set(mtflist) & set(interactions.keys()))
    return [ target for name in mylist for target in interactions[name] ]

# receive a list of miRs or TFs (mtflist) and a dict of target interactions
# return interactions for members of mtflist
def trimInteractions(mtflist, interactions):
    mylist = list(set(mtflist) & set(interactions.keys()))
    return dict( (mtf,interactions[mtf]) for mtf in mylist )

# Find targets that a miR regulates both directly and indirectly,
# through a TF. Returns a dict { miR : [(tf, target)] }
def getDoubleInteractions(mirlist, mirdict, tfdict):
    doubles = dict()
    for mir, tfs in mirdict.iteritems():
        doubles[mir.lower()] = list()
        for tf in tfs:
            try:
                overlap = set(tfs) & set(tfdict[tf])
            except KeyError:
                continue
            doubles[mir.lower()].extend([ (tf,target) for target in overlap ])
    return dict( (mir,interaction) 
                 for mir, interaction in doubles.iteritems() if interaction )
 
# Get the names of samples in expression matrix. The second line 
# removes the barcode information after the sample number to 
# facilitate comparisons between files.
def getSampleNames(f):
    columns = f.readline().strip().split()
    return [ name[0:15] for name in columns ]

# Find the individuals who have comparatively high and low expression
# of a specific mirna, two lists of values
def selectOutliers(expression):
    top_cutoff = scoreatpercentile(expression, 80)
    bottom_cutoff = scoreatpercentile(expression, 20)
    top = [ val for val in expression
            if val > top_cutoff ]
    bot = [ val for val in expression
            if val < bottom_cutoff ]
    return top, bot

# Find indices of a list of values from a line in expression matrix
def getIndices(values, expr):
    return [ expr.index(val) for val in values ]

# Map the names of targets and transcription factors to lists of their
# expression data
def getExpressionData(mirfile, rnafile, interactions):
    # Set file pointers to beginning of data
    rnafile.seek(0)
    mirfile.seek(0)
    rnafile.readline()
    mirfile.readline()
    # Initialize dictionary with targets and tfs
    targets = dict()
    for mir, pairs in interactions.iteritems():
        targets[mir] = []
        for pair in pairs:
            targets[pair[0]], targets[pair[1]] = [], []
    # Find RNA expression data
    for line in rnafile:
        line = line.strip().split()
        if proteinName(line[0]) in targets:
            targets[proteinName(line[0])] = map(float,line[1:])
    # Find miR expression data
    for line in mirfile:
        line = line.strip().split()
        if line[0] in targets:
            targets[line[0]] = map(float,line[1:])
    return targets

# If the expression values for a compound consists entirely of zeroes, it will 
# be impossible to use gaussian_kde to infer its probability distribution.
# Use this function to filter out those instances
def notExpressed(expr):
    if len(set(expr)) == 1:
        return True
    return False

def proteinName(name):
    return name.split('|')[0]

# Once miRNA outliers are identified, find those individuals expression
# of a given TF or target
def getOutlierCoexpression(outlier_ids,exp,id_map):    
    top_rna = [ exp[id_map[i]] for i in outlier_ids[0] if id_map[i] ]
    bot_rna = [ exp[id_map[i]] for i in outlier_ids[1] if id_map[i] ]
    return top_rna, bot_rna
    
# Map the indexes of samples between different expression matrices.
def mapCorrespondingIndices(names1, names2):
    id_map = [None]*len(names1)
    for name in set(names1) & set(names2):
        id_map[names1.index(name)] = names2.index(name)
    return id_map

# Calculate the mutual information between two vectors
def mutualInformation(X,Y):
    # Use a gaussian kernel estimator to approximate the pdfs
    pX = gaussian_kde(X)
    pY = gaussian_kde(Y)
    # Estimate joint pdf
    pXY = gaussian_kde([X,Y])
    # Use estimated distributions to approx. entropies
    sX = entropy(pX.evaluate(X))
    sY = entropy(pY.evaluate(Y))
    sXY = entropy(pXY.evaluate([X,Y]))
    # Calculate and return mutual information between X and Y
    MI = sX + sY - sXY
    return MI


# -------- MAIN / TESTING ------------------------------------------------------

# List of miRNAs to test
mtflist = getMirList(files.mirlist)

# Unpack dictionaries containing {tf:[target]} and {mir:[target]} interactions
tfdict = unpickle(files.rnapkl)
mirdict = trimInteractions(mtflist,
                           unpickle(files.mirpkl))

# Find double interactions and write them to a pickle
double_interactions = getDoubleInteractions(mtflist,mirdict,tfdict)
d = open(files.outdir + "double_interactions.pkl", 'w')
pickle.dump(double_interactions, d)
d.close()

# miRNA and mRNA expression matrices
mirfile, rnafile = filterExpressionFiles(files.mirdata, files.rnadata)
mirfile = open(mirfile, 'r')
rnafile= open(rnafile, 'r')

# Open output file
outfile = files.outdir + "MI_double_results.txt"
o = open(outfile, 'w')

# Get sample names from both expression matrices
mir_samples = getSampleNames(mirfile)
rna_samples = getSampleNames(rnafile)

# Map the sample columns in the mir expression file to those in the
# rna expression file
id_map = mapCorrespondingIndices(mir_samples,rna_samples)

expression = getExpressionData(mirfile,rnafile,double_interactions)
expression = dict( (target, exp) 
                   for target, exp in expression.iteritems()
                   if exp )

for mir, interactions in double_interactions.iteritems():
    # Find expression outliers and their IDs
    mirexp = expression[mir]
    if not mirexp:
        o.write("\n" + mir + "\nMISSING MIR EXPRESSION DATA.\n")
        continue
    top, bot = selectOutliers(mirexp)
    top_id, bot_id = getIndices(top, mirexp), getIndices(bot, mirexp)

    # For every (mir, tf, t) trio in double_interactions, calculate the
    # change in mutual information between tf and t for high mir expression
    # and low mir expression
    for (tf, t) in interactions:
        if tf == t:
            continue
        o.write("\t".join([mir,tf,t])+"\t")
        try:
            top_tf, bot_tf = getOutlierCoexpression((top_id,bot_id),
                                                    expression[tf],
                                                    id_map)
            top_t, bot_t = getOutlierCoexpression((top_id,bot_id),
                                                  expression[t],
                                                  id_map)
        except KeyError:
            o.write("MISSING RNA EXPRESSION DATA.\n")
            continue
        
        top_MI = mutualInformation(top_tf, top_t)
        bot_MI = mutualInformation(bot_tf, bot_t)
    
        o.write(str(top_MI-bot_MI) + "\n")
o.close()
