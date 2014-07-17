Scripts for calculating the Mutual Information for miRNAs and their targets. `MI.py` contains general functions used by `MI_DoubleInteraction.py` and `MI_IndirectInteraction.py`, which calculate the MI for different cases (see below). Both of them receive the following items:
* `.pkl`s containing interactions between miRNAs and their direct targets, and interactions between gene products
* Matrices containing expression data for both miRNAs and mRNA in a specific type of cancer tissue.
* A file containing the name of a miRNA in each line

**`DoubleInteraction`** By a "double interaction", I mean a (miRNA, transcription factor, target) tuple where the miRNA directly targets both the transcription factor and the target, while the transcription factor also regulates the target (is there a name for this?). For each such tuple, the script finds the change in mutual information of the expression of the transcription factor and the target between samples where the miRNA expression is either above or below a certain threshold.   

**`IndirectInteraction`** (not implemented yet) In this case, we consider (miRNA, transcription factor, target) tuples where the miRNA targets the transcription factor, the transcription factor regulates the targets, but there is no direct interaction between the miRNA and the target. Here we want to find the mutual information of the miRNA and the target between different expression levels of the transcription factor.

### NOTES ###
I implemented `IndirectInteraction` but it currently takes too long to run. In part, it's unavoidable, since it evaluates several times the number of interactions evaluated by `DoubleInteraction`. However, it's also taking a really long time because of redundancies in calculating the Mutual Information. I am calculating the probability distributions and the entropy from within the loop, instead of using the function, as an attempt to solve this problem, but the last time I tried to run it, it still took over an hour for the COAD dataset.


### PENDING ###
* Finishing implementation of 'IndirectInteraction' today (17-July-2014)
* Adding visualization scripts this week
