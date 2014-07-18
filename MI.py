import cPickle as pickle
from tempfile import NamedTemporaryFile
from itertools import izip
from scipy.stats import entropy, scoreatpercentile 
from scipy.stats.kde import gaussian_kde

# -------- Argument Parsing ----------------------------------------------------
from optparse import OptionParser
parser = OptionParser()
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
parser.add_option("--mirlist",dest="mirlist",
                  help="(optional) plain text file with the name of a miRNA in each line.")
(files, args) = parser.parse_args()

if not files.outdir.endswith("/"):
    files.outdir = files.outdir + "/"

# -------- Functions -----------------------------------------------------------

# extract an object from a pickle
def unpickle(pklfile):
    f = open(pklfile, 'rb')
    pkl = pickle.load(f)
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

    mir_out, rna_out = NamedTemporaryFile(), NamedTemporaryFile()
               
    mir_out.write("\t".join([mir_samples[i] for i in mir_indices])+"\n")
    rna_out.write("\t".join([rna_samples[i] for i in rna_indices])+"\n")

    for f,ids,o in izip([mir, rna],[mir_indices,rna_indices],[mir_out,rna_out]):
        for line in f:
            line = line.strip().split()
            out = [line[0]]
            out.extend([ line[i+1] for i in ids ])
            o.write("\t".join(out)+"\n")
        o.seek(0)
    return mir_out, rna_out

# Get the names of samples in expression matrix. The second line 
# removes the barcode information after the sample number to 
# facilitate comparisons between files.
def getSampleNames(f):
    columns = f.readline().strip().split()
    return [ name[0:15] for name in columns ]

# Find indices of a list of values from a line in expression matrix
def getIndices(values, expr):
    return [ expr.index(val) for val in values ]

# Map the indexes of samples between different expression matrices.
def mapCorrespondingIndices(names1, names2):
    id_map = [None]*len(names1)
    for name in set(names1) & set(names2):
        id_map[names1.index(name)] = names2.index(name)
    return id_map

def proteinName(name):
    return name.split('|')[0]

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

# Is the target expressed?
def isExpressed(expr):
    return not(not expr or len(set(expr))==1)

# Find the individuals who have comparatively high and low expression
# of a specific product, two lists of values
def selectOutliers(expression):
    top_cutoff = scoreatpercentile(expression, 70)
    bottom_cutoff = scoreatpercentile(expression, 30)
    top = [ val for val in expression
            if val > top_cutoff ]
    bot = [ val for val in expression
            if val < bottom_cutoff ]
    return top, bot

# Once miRNA outliers are identified, find those individuals expression
# of a given TF or target
def getOutlierCoexpression(outlier_ids,exp,id_map):    
    top_rna = [ exp[id_map[i]] for i in outlier_ids[0] if id_map[i] ]
    bot_rna = [ exp[id_map[i]] for i in outlier_ids[1] if id_map[i] ]
    return top_rna, bot_rna

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

# Unpack dictionaries of interactions: {mir:[targets]} and {rna:[targets]}
mirdict = unpickle(files.mirpkl)
tfdict = unpickle(files.rnapkl)

# Filter dictionary of miR interactions using mirlist, if provided
if files.mirlist:
    f = open(files.mirlist)
    mirlist = [ line.strip() for line in f ]
    mirdict = dict( (mir.lower(),targets) 
                    for mir, targets in mirdict.iteritems()
                    if mir.lower() in mirlist)
    f.close()

# Filter miRNA and mRNA expression matrices for shared samples
mirfile, rnafile = filterExpressionFiles(files.mirdata, files.rnadata)

# Get sample names from both expression matrices
mir_samples = getSampleNames(mirfile)
rna_samples = getSampleNames(rnafile)

# Map the sample columns in the mir expression file 
# to those in the rna expression 
id_map = mapCorrespondingIndices(mir_samples,rna_samples)
