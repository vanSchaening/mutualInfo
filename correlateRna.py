from scipy.stats import pearsonr
import sys

infile = sys.argv[1]#"/home/axolotl/Data/ExpressionMatrices/COAD_20140416_RNASeqV2.txt"
outfile = sys.argv[2]#"/home/axolotl/Data/mRNAcorrelations/COAD_20140616-mRNA.txt"
pvals = sys.argv[3]#"/home/axolotl/Data/mRNAcorrelations/COAD_20140616-pVals.txt"

if len(sys.argv) != 4:
    print "Usage:\n\t$ python correlateRna.py <input file> <correlation output file> <pvalues output file>\n"

# Open input file twice
f = open(infile, 'r')
sup = open(infile, 'r')

def getExpression(line):
    line = line.strip().split()
    name = line[0].split("|")[0]
    expr = map( float, line[1:] )
    return name, expr

# Set file pointers at the beginning of data
sup.readline()
f.readline()
data_ptr = f.tell()

# Open output files
o = open(outfile,'w')
p = open(pvals, 'w')

# Get RNA names
names = [ line.strip().split()[0] for line in f ]
names = [ name.split("|")[0] for name in names ]
f.seek(data_ptr)

# Write column names to files
o.write("\t".join(names)+"\n")
p.write("\t".join(names)+"\n")

# Get correlation between RNA in infile and write them to outfile,
# Simultaneously write corresponding p values to file.
line = f.readline()
while line != '':
    gene1, expr1 = getExpression(line)
    oline = [ gene1 ]
    pline = [ gene1 ]
    
    olist, plist = [], []

    for expr2 in sup:
        gene2, expr2 = getExpression(expr2)
        corr,prob = pearsonr(expr1,expr2)
        olist.append(str(corr))
        plist.append(str(prob))
  
    # Matrix is symmetrical. Use placeholders on left triangles to
    # to avoid redundant computations
    filler = ["-2"]*(len(names)-len(olist))
    oline.extend(filler)
    pline.extend(filler)

    # Write calculated  correlations
    oline.extend(olist)
    pline.extend(plist)
    
    o.write("\t".join(oline)+"\n")
    p.write("\t".join(pline)+"\n")
    
    # Set sup's pointer back in place
    sup.seek(data_ptr,0)
    line = f.readline()
    data_ptr = f.tell()

f.close()
o.close()
p.close()

