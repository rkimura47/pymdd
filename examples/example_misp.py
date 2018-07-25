#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Maximum Independent Set Problem Example from Section 3.5 of Bergman et. al, Decision Diagrams for Optimization
G = [[(1,2), (1,3)], [(2,1), (2,3), (2,4)], [(3,1), (3,2), (3,4)], [(4,2), (4,3), (4,5)], [(5,4)]]
neighbors = [ [v for (u,v) in G[j]] for j in range(5) ]
weight = [3, 4, 2, 2, 7]
numLayers = 5
domain = lambda j: (0,1)
rootState = frozenset([1,2,3,4,5])
def trFunc(s,d,j):
    if d == 1:
        if j+1 in s:
            return s - {j+1} - set(neighbors[j])
        else:
            return None
    else:
        return s - {j+1}
costFunc = lambda s,d,j: d*weight[j]
isFeas = lambda s,j: s is not None
mergeFunc = lambda slist,j: frozenset(set.intersection(*(set(s) for s in slist)))
adjFunc = lambda w,os,ms,j: w
name = 'misp'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down compilation
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
mymdd.reduce_bottom_up(mergeFunc, adjFunc)

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Display MDD with GraphViz
mymdd.output_to_dot(nodeDotFunc=lambda s: '[label="{' + ', '.join(str(i) for i in s)  + '}"];')
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
