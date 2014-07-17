import sys

if len(sys.argv) == 1:
    exit("Usage:\n\tgetMirlist.py filename1 [filename2 ... filename n]\n\nWill return a list of the mirs expressed in all files.\n")

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

# Obtain list of miRs from the first file
mirlist = getMirsFromFile(sys.argv[1])
# For every subsequent miR, get a list of miRs, and remove
# miRs in mirlist that aren't shared with the new list
for mirfile in sys.argv[2:]:
    mirlist = set(mirlist) & set(getMirsFromFile(mirfile))

o = open("shared_mirs.txt", 'w')
o.write("\n".join(mirlist))
