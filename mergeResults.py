# Receives a list of directories, each of which contains results from MI
# and a name for an output file without an extension
import sys

def formatDir(dirname):
    if not dirname.endswith("/"):
        return dirname+"/"
    return dirname

out_pos = open(sys.argv[-1]+"_coherent.txt",'w')
out_neg = open(sys.argv[-1]+"_incoherent.txt",'w')

indirs = sys.argv[1:sys.argv.index("--out")]
indirs = map(formatDir, indirs)

header = ["interaction"]
header.extend([indir.strip("/") for indir in indirs])
out_pos.write("\t".join(header)+"\n")
out_neg.write("\t".join(header)+"\n")

interactions = dict()
for directory in indirs:
    f = open(directory+"MI_double_results_RNA_corr.txt",'r')
    f.readline()
    for line in f:
        line = line.strip().split()
        try:
            interactions[(line[0],line[1],line[2])][directory.strip("/")] = (float(line[3]), float(line[4]))
        except KeyError:
            interactions[(line[0],line[1],line[2])] = { directory.strip("/"):(float(line[3]), float(line[4])) }
    f.close()

indirs = [ indir.strip("/") for indir in indirs ]
for (mir,tf,t),results in interactions.iteritems():
    pos_line = [",".join([mir,tf,t])]
    neg_line = [",".join([mir,tf,t])]
    for cancer in indirs:
        if cancer not in results:
            pos_line.append(0)
            neg_line.append(0)
        else:
            mutualinfo, corr = results[cancer]
            if corr > 0:
                pos_line.append(mutualinfo)
                neg_line.append(0)
            else:
                neg_line.append(mutualinfo)
                pos_line.append(0)
    out_pos.write("\t".join(map(str,pos_line))+"\n")
    out_neg.write("\t".join(map(str,neg_line))+"\n")

out_pos.close()
out_neg.close()
