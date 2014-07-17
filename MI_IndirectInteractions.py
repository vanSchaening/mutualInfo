from MI import *

def flipInteractions(interactions):
    flipped = dict()
    for mir, pairs in interactions.iteritems():
        for (tf, target) in pairs:
            if tf not in flipped:
#                flipped[tf] = list()
#            flipped[tf].append((mir,target))
                flipped[tf] = dict()
            if mir not in flipped[tf]:
                flipped[tf][mir] = list()
            flipped[tf][mir].append(target)
   #         except KeyError:
    #            flipped[tf][mir] = [target]
    return flipped

# -------- MAIN ----------------------------------------------------------------

# ---- Get indirect interactions ----
indirect_interactions = dict()
for mir, rnas in mirdict.iteritems():
    mir = mir.lower()
    overlap = set(rnas) & set(tfdict.keys())
    indirect_interactions[mir] = list()
    for tf in overlap:
        indirect_interactions[mir].extend([ (tf, target)
                                            for target in tfdict[tf] 
                                            if target not in rnas ])

# ---- Get expression data ----
expression = getExpressionData(mirfile, rnafile, indirect_interactions)
expression = dict( (target, exp)
                   for target, exp in expression.iteritems()
                   if isExpressed(exp))

# ---- Update interaction dictionary, remove parts with missing data ----
indirect_interactions = dict( (mir, indirect_interactions[mir])
                              for mir in indirect_interactions.keys()
                              if mir in expression )
for mir, interactions in indirect_interactions.iteritems():
    indirect_interactions[mir] = [ (tf, target) for tf, target in interactions
                                   if tf in expression and target in expression ]

# ---- Flip the dictionary into {tf:[(mir,target)]} ----
indirect_interactions = flipInteractions(indirect_interactions)

# ---- DEBUG ----
#test = indirect_interactions.keys()[7]
#print indirect_interactions[test]

# ---- Open output file ----
o = open(files.outdir + "MI_indirect_results.txt", 'w')

# ---- Calculate MI for all (miR,tf,target) interactions ----
for tf, interactions in indirect_interactions.iteritems():
    
    tfexp = expression[tf]
    top_tf, bot_tf = selectOutliers(tfexp)
    if not (isExpressed(top_tf) and isExpressed(bot_tf)):
        continue
    top_id, bot_id = getIndices(top_tf,tfexp), getIndices(bot_tf,tfexp)

    for mir, targets in interactions.iteritems():
        top_mir, bot_mir = getOutlierCoexpression([top_id,bot_id],
                                                  expression[mir],
                                                  id_map)
        if not (isExpressed(top_mir) and isExpressed(bot_mir)):
            continue

        # Estimate probability distributions for the outliers' mir expression
        p_top_mir = gaussian_kde(top_mir)
        p_bot_mir = gaussian_kde(bot_mir)

        # Find the entropies for both
        s_top_mir = entropy(p_top_mir.evaluate(top_mir))
        s_bot_mir = entropy(p_bot_mir.evaluate(bot_mir))
        
        for target in targets:
            top_t, bot_t = getOutlierCoexpression([top_id,bot_id],
                                                  expression[target],
                                                  id_map)
            if not (isExpressed(top_t) and isExpressed(bot_t)):
                continue

            # Estimate probability distributions for the outliers' mir expression
            p_top_t = gaussian_kde(top_t)
            p_bot_t = gaussian_kde(bot_t)

            # Find the entropies for both
            s_top_t = entropy(p_top_t.evaluate(top_t))
            s_bot_t = entropy(p_bot_t.evaluate(bot_t))

            # Joint probability and joint entropy
            p_top_joint = gaussian_kde([top_mir,top_t])
            p_bot_joint = gaussian_kde([bot_mir,bot_t])
            s_top_joint = entropy(p_top_joint.evaluate([top_mir,top_t]))
            s_bot_joint = entropy(p_bot_joint.evaluate([bot_mir,bot_t]))

            # Conditional Mutual Information
            top_MI = s_top_t + s_top_mir - s_top_joint
            bot_MI = s_bot_t + s_bot_mir - s_bot_joint
#    for (mir, target) in interactions:
        # Get expression of mir and t in the outliers
#        top_mir, bot_mir = getOutlierCoexpression([top_id,bot_id],
#                                                  expression[mir],
#                                                  id_map)
#        if not (isExpressed(top_mir) and isExpressed(bot_mir)):
#            continue
#        top_t, bot_t = getOutlierCoexpression([top_id,bot_id],
#                                              expression[target],
 #                                             id_map)
  #      if not (isExpressed(top_t) and isExpressed(bot_t)):
   #         continue
    #    # Calculate conditional mutual information

            #at this rate, p_mir is being calculated redundantly. Find
# a way to estimate the kernel density outside the function and instead
# just pass the approximated distribution
#            top_MI = mutualInformation(top_mir, top_t)
 #           bot_MI = mutualInformation(bot_mir, bot_t)

            o.write("\t".join([mir,tf,target,
                               str(top_MI - bot_MI)])+"\n")
o.close()
