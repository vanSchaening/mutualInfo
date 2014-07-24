# Get input file from command line and open it
import sys
f = open(sys.argv[1])
# Output files will be saved in the same directory as input file
outdir = sys.argv[1].split("/")
outdir.pop()
outdir = "/".join(outdir)+"/"

# Get list of (MI,corr) points with corr past the threshold
plus_points = []
minus_points = []                
for line in f:
    line = line.strip().split()
    if float(line[4]) >= 0:
        plus_points.append( (float(line[3]),float(line[4])) ) 
    elif float(line[4]) < 0:
        minus_points.append( (float(line[3]), float(line[4])) )

# Change the points into two vectors
plus_points = zip(*plus_points)
minus_points = zip(*minus_points)

# Use plotting library with non-interactive backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Plot points with a positive correlation
plt.figure(0)
plt.scatter(plus_points[0],plus_points[1])
plt.xlabel("Change in MI between TF and target, with respect to miR.")
plt.ylabel("Correlation between TF and target")
plt.savefig(outdir+"positive_RNAcorr.png")

# Plot points with a negative correlation
plt.figure(1)
plt.scatter(minus_points[0],minus_points[1])
plt.xlabel("Change in MI between TF and target, with respect to miR.")
plt.ylabel("Correlation between TF and target")
plt.savefig(outdir+"negative_RNAcorr.png")

f.close()
