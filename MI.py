
# This file contains argument parsing and functions to be used both by 
# DoubleInteractions and IndirectInteractions, and is imported by both of them.

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

# extract a dictionary from a .pkl
def unpickle(pklfile):
    f = open(pklfile, 'rb')
    pkl = pickle.load(f)
    f.close()
    return pkl

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

# Get the names of samples in expression matrix. The second line 
# removes the barcode information after the sample number to 
# facilitate comparisons between files.
def getSampleNames(f):
    columns = f.readline().strip().split()
    return [ name[0:15] for name in columns ]

# Find the individuals who have comparatively high and low expression
# of a specific mirna, two lists of values
def selectOutliers(expression):
    top_cutoff = scoreatpercentile(expression, 70)
    bottom_cutoff = scoreatpercentile(expression, 30)
    top = [ val for val in expression
            if val > top_cutoff ]
    bot = [ val for val in expression
            if val < bottom_cutoff ]
    return top, bot

# Find indices of a list of values from a line in expression matrix
def getIndices(values, expr):
    return [ expr.index(val) for val in values ]

# If the expression values for a compound consists entirely of zeroes, it will 
# be impossible to use gaussian_kde to infer its probability distribution.
# Use this function to filter out those instances
def isExpressed(expr):
    if len(set(expr)) <= 1:
        return False
    return True

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
