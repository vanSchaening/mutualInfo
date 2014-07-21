from scipy.stats import pearsonr
import sys

infile = sys.argv[1]#"/home/axolotl/Data/ExpressionMatrices/COAD_20140416_RNASeqV2.txt"
outfile = sys.argv[2]#"/home/axolotl/Data/mRNAcorrelations/COAD_20140616-mRNA.txt"
pvals = sys.argv[3]#"/home/axolotl/Data/mRNAcorrelations/COAD_20140616-pVals.txt"

# Open input file twice
f = open(infile, 'r')
sup = open(infile, 'r')

def getExpression(line):
    line = line.strip().split()
    name = line[0].split("|")[0]
    expr = map( float, line[1:] )
    return name, expr

# Set file pointers at the beginning of data
bufsize = len(f.readline()) * 128
sup.readline()
data_ptr = f.tell()

# Open output files
o = open(outfile,'w',bufsize)
p = open(pvals, 'w',bufsize)

# Get RNA names
names = [ line.strip().split()[0] for line in f ]
names = [ name.split("|")[0] for name in names ]
f.seek(data_ptr)

# Write column names to files
o.write("\t".join(names)+"\n")
p.write("\t".join(names)+"\n")

# Get correlation between RNA in infile and write them to outfile,
# Simultaneously write corresponding p values to file.
for expr1 in f:
    gene1, expr1 = getExpression(expr1)
    oline = [ gene1 ]
    pline = [ gene1 ]
    for expr2 in sup:
        gene2, expr2 = getExpression(expr2)
        corr,prob = pearsonr(expr1,expr2)
        oline.append(str(corr))
        pline.append(str(prob))

    o.write("\t".join(oline)+"\n")
    p.write("\t".join(pline)+"\n")

    # Set sup's pointer back in place
    sup.seek(data_ptr)

# I'm currently doing a lot of redundant computation:
    # Each correlation is being computed twice, since the result is 
    # a symmetric matrix. 

f.close()
o.close()

