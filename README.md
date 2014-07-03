Most of the scripts contained in this repository were just exploratory ones that I wrote while getting familiar with the project data. Currently, the important (read: useful) ones are `MI_DoubleInteraction.py` and `MI_IndirectInteraction.py`. Both of them receive the following items:
* `.pkl`s containing interactions between miRNAs and their direct targets, and interactions between gene products
* Matrices containing expression data for both miRNAs and mRNA in a specific type of cancer tissue.
* A file containing the name of a miRNA in each line

*DoubleInteraction* By a "double interaction", I mean a (miRNA, transcription factor, target) tuple where the miRNA directly targets both the transcription factor and the target, while the transcription factor also regulates the target (is there a name for this?). For each such tuple, the script finds the change in mutual information of the expression of the transcription factor and the target between examples where the miRNA expression is either above or below a certain threshold.   

*IndirectInteraction*(not implemented yet)


### PENDING ###
* While the functions in `MI_DoubleInteraction` have already been tested, I still need to implement an argument parser and a main body to have it run for various datasets.
* Many of the functions I wrote will be reused when implementing `MI_IndirectInteraction`, so I'll probably move those to a separate module.
* I still need to come up with a good way to visualize the data produced by the scripts