# Some code I wrote when I was looking at the expression matrices
# I was figuring out the structure of the TCGA barcodes and 
# how much of the data I would be able to use for correlations / MI

p = "/home/axolotl/Data/ExpressionMatrices/"

mir = open(p+"BRCA_20140316_miRNASeq.txt", 'r')
rna1 = open(p+"BRCA_20140316_RNASeq.txt", 'r') 
rna2 = open(p+"BRCA_20140316_RNASeqV2.txt", 'r')

def getPatientIds(f):
    line = f.readline().strip().split()
    return [ ".".join(name.split('.')[1:3]) for name in line]

samples = {'mir':getPatientIds(mir),
           'rna1':getPatientIds(rna1),
           'rna2':getPatientIds(rna2)}

#print samples['mir']

print
print len(set(samples['rna2']) & set(samples['mir'])), "shared samples."
print len(set(samples['rna2'])), "unique samples from RNA-Seq."
print len(set(samples['mir'])), "unique miR samples."
print
print len(samples['rna2']), "total samples from RNA-Seq."
print len(samples['mir']), "total miR samples."
print
