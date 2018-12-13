#!/usr/bin/python3
from pymdd.mdd import *
from os import system
from random import shuffle

# Script parameters
numVars = 18
showGraph = False

# Threshold function counterexample (Hosaka et al 1997)
# Corresponding MDD has width >= 2^(sqrt(n)/2) in *any* variable ordering

# Set k <- ceil(sqrt(numVars))
k = 0
while k*k < numVars:
    k += 1
# Set weights and threshold in a particular way
weight = {(i,j): 2**i + 2**(j+k) for i in range(k) for j in range(k)}
varorder = list(weight.keys())
threshold = sum(weight[varorder[i]] for i in range(numVars))/2

# Change varorder randomly
shuffle(varorder)

# MDD parameters
numLayers = numVars
domain = lambda i: (0,1)
rootState = threshold
trFunc = lambda s,d,i: s-d*weight[varorder[i]]
costFunc = lambda s,d,i,ns: 0
isFeas = lambda s,i: s >= 0
mergeFunc = lambda slist,i: min(slist)
name = 'threshold'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down compilation
print('Compiling...')
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
print('Pruning...')
mymdd.prune_all()
print('Reducing...')
mymdd.reduce_bottom_up(mergeFunc)

# Print the contents
#print(mymdd)
print(mymdd.maxWidth)

# Output and display MDD with GraphViz
if showGraph and mymdd.maxWidth <= 50:
    ndf = lambda s,j: '[label=""];'
    adf = lambda l,w,j: '[%slabel="%d"];' % ('style=dotted,' if l == 0 else '',l)
    mymdd.output_to_dot(ndf, adf, MDD._default_arcsortargs)
    system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
    system('gv ' + name + '.ps')
