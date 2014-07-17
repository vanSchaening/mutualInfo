import sys

# Get the list of miRs in an expression matrix
def getMirsFromFile(exprfile):
    f = open(exprfile, 'r')
    f.readline()
    mirlist = list()
    for line in f:
        line = line.strip().split()
        if len(set(line[1:])) <= 1:
            continue
        mirlist.append(line[0])
    f.close()
    return mirlist

# ---- Input validation ----
if len(sys.argv) <= 3 or ("--out" not in sys.argv):
    exit("Usage:\n\tgetMirlist.py filename1 [filename2 ... filename n] --out outfile\n\nWill return a list of the mirs expressed in all files.\n")

# ---- Find and open output file ----
outpos = sys.argv.index("--out")
outfile = sys.argv[outpos+1]
o = open(outfile, 'w')

# ---- Obtain list of miRs ----
mirlist = set(getMirsFromFile(sys.argv[1]))
#  For every file after the first one, get a list of miRs, and remove
#  miRs in the original list that aren't shared with the new list
for mirfile in sys.argv[2:outpos]:
    mirlist = mirlist & set(getMirsFromFile(mirfile))

o.write("\n".join(mirlist))
o.close()
