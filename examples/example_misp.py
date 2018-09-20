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
costFunc = lambda s,d,j,ns: d*weight[j]
isFeas = lambda s,j: s is not None
maxWidth = lambda j: 2
nodeSelFunc = lambda vlist,j: [min(vlist), max(vlist)]
mergeFunc = lambda slist,j: frozenset(set.union(*(set(s) for s in slist)))
adjFunc = lambda w,os,ms,j: w
name = 'misp'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down exact compilation
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
# Merge nodes
#mymdd.merge_nodes([MDDNode(2,frozenset([5])), MDDNode(2,frozenset([3,4,5]))], lambda l: mergeFunc(l,2))

# Perform DP-based top-down relaxed compilation
#mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas, maxWidth, mergeFunc, adjFunc, nodeSelFunc)

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Output and display MDD with GraphViz
mymdd.output_to_dot(nodeDotFunc=lambda s,j: '[label="{' + ', '.join(str(i) for i in s)  + '}"];')
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
