from pprint import pprint

def getSampleNames(f):
    columns = f.readline().strip().split()
    f.seek(0)
    return [ name[0:15] for name in columns ]

def filterExpressionFiles(mirfile,rnafile):
    mir = open(mirfile, 'r')
    rna = open(rnafile, 'r')

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

    mir_out = open(mir_file, 'w')
    rna_out = open(rna_file, 'w')
               
    mir_out.write("\t".join([mir_samples[i] for i in mir_indices])+"\n")
    rna_out.write("\t".join([rna_samples[i] for i in rna_indices])+"\n")

    for line in mir:
        line = line.split()
        out = [line[0]]
        out.extend([ line[i+1] for i in mir_indices ])
        mir_out.write("\t".join(out)+"\n")
    mir.close()
    mir_out.close()

    for line in rna:
        line = line.split()
        out = [line[0]]
        out.extend([ line[i+1] for i in rna_indices ])
        rna_out.write("\t".join(out)+"\n")
    rna.close()
    rna_out.close()

    return mir_file, rna_file

filterExpressionFiles("/home/axolotl/Data/ExpressionMatrices/COAD_20140416_miRNASeq.txt",
                      "/home/axolotl/Data/ExpressionMatrices/COAD_20140416_RNASeqV2.txt")
