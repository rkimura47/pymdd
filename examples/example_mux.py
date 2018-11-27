#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Script parameters
# NOTE: n = 4 takes a few minutes to compile.
n = 3
reverseOrdering = False

# 2^n-to-1 Multiplexer
# Point (S_0, ..., S_{n-1}, D_0, ..., D_{(2^n)-1}) is feasible if
# D_{<S_0...S_{n-1}>_2} == 1
#     where <X_0...X_{n-1}>_2 = sum(X_i * 2^i for i = 0, ..., n-1)
# i.e., output is equal to value of data bit D_k where k is the index
# selected by the selector bits S_0, ..., S_{n-1}
#
# This is one example where the variable ordering has a significant effect
# on the MDD width.
# Using variable ordering (S_0, ..., S_{n-1}, D_0, ..., D_{(2^n)-1})
# results in an MDD with width 2^n,
# but the reverse variable ordering (D_{(2^n)-1}, ..., D_0, S_{n-1}, ..., S_0)
# results in an MDD with width 2^(2^n) - 1.

# Multiplexer parameters
numSel = n
numWays = 2**numSel
varorder = [(0, i) for i in range(numSel)] + [(1, i) for i in range(numWays)]
if reverseOrdering:
    varorder.reverse()

# MDD-DP parameters
numLayers = numSel + numWays
domain = lambda i: (0,1)
# state = (sum of selection bits, sum of data bits)
rootState = (0,0)
# transition -> s[varorder[i][0]] += d*(2**varorder[i][1])
trFunc = lambda s,d,i: ( s[0] + (1-varorder[i][0])*d*(2**varorder[i][1]), s[1] + varorder[i][0]*d*(2**varorder[i][1]) )
costFunc = lambda s,d,i,ns: 0
isFeas = lambda s,i: i < numLayers or (s[1] >> s[0])%2 == 1
mergeFunc = lambda slist,i: slist[0]
name = 'mux'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down compilation
if n > 3:
    print('Compiling...')
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
if n > 3:
    print('Pruning...')
mymdd.prune_all()
if n > 3:
    print('Reducing...')
mymdd.reduce_bottom_up(mergeFunc)

# Print the contents
print(mymdd)
print(mymdd.maxWidth)

# Output and display MDD with GraphViz
if mymdd.maxWidth <= 50:
    ndf = lambda s,j: '[label=""];'
    adf = lambda l,w,j: '[%slabel="%d"];' % ('style=dotted,' if l == 0 else '',l)
    mymdd.output_to_dot(ndf, adf, MDD._default_arcsortargs)
    system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
    system('gv ' + name + '.ps')
