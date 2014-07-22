# MI is a script that performs actions needed for both double and indirect
# interactions. From it, import:
#     mirdict := {mir:[targets]} interactions
#     tfdict :=  {rna:[targets]} interactions
#     mirfile, rnafile := expression matrices, with corresponding samle names
#     mir_samples, rna_samples
#     id_map := simple map of names in mir_samples to names in rna_samples
# as well as several methods

from MI import *

# -------- MAIN ----------------------------------------------------------------

# ---- Get double interactions ----
double_interactions = dict()
for mir,rnas in mirdict.iteritems():
    mir = mir.lower()
    double_interactions[mir] = list()
    for tf in rnas:
        try:
            overlap = set(rnas) & set(tfdict[tf])
        except KeyError:
            continue
        double_interactions[mir].extend([ (tf,target) for target in overlap ])
double_interactions = dict( (mir, pairs)
                            for mir, pairs in double_interactions.iteritems()
                            if pairs )

# ---- Get expression data ----
expression = getExpressionData(mirfile, rnafile, double_interactions)
expression = dict( (target, exp)
                   for target, exp in expression.iteritems()
                   if isExpressed(exp) )

# ---- Update interactions dictionary, remove mirs with no data ----
double_interactions = dict( (mir, double_interactions[mir])
                            for mir in double_interactions.keys()
                            if mir in expression )

# ---- Open output file ----
o = open(files.outdir + "MI_double_results.txt", 'w')

# ---- Calculate MI for all (miR,tf,target) interactions ----
for mir, interactions in double_interactions.iteritems():
    # Find the outliers in miR expression and their corresponding IDs
    mirexp = expression[mir]
    top_mir, bot_mir = selectOutliers(mirexp)
    if not (isExpressed(top_mir) and isExpressed(bot_mir)):
        continue
    top_id, bot_id = getIndices(top_mir,mirexp), getIndices(bot_mir,mirexp)

    # Calculate the change in conditional MI between the TF and target
    for (tf,target) in interactions:
        # Skip missing/invalid targets
        if tf == target:
            continue
        if tf not in expression or target not in expression:
            continue

        # Get expression of tf and t in outliers
        top_tf, bot_tf = getOutlierCoexpression([top_id,bot_id],
                                                expression[tf],
                                                id_map)
        if not (isExpressed(top_tf) and isExpressed(bot_tf)):
            continue
        top_t, bot_t = getOutlierCoexpression([top_id,bot_id],
                                              expression[target],
                                              id_map)
        if not (isExpressed(top_t) and isExpressed(bot_t)):
            continue

        # Calculate the conditional mutual information
        top_MI = mutualInformation(top_tf, top_t)
        bot_MI = mutualInformation(bot_tf, bot_t)

        # Write names and delta MI to file
        o.write("\t".join([mir, tf, target,
                           str(top_MI-bot_MI)]) + "\n")
o.close()
