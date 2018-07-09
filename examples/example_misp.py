#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Maximum Independent Set Problem Example from Section 3.5 of Bergman et. al, Decision Diagrams for Optimization
G = [[(0,1), (0,2)], [(1,0), (1,2), (1,3)], [(2,0), (2,1), (2,3)], [(3,1), (3,2), (3,4)], [(4,3)]]
neighbors = [ [v for (u,v) in G[i]] for i in range(5) ]
weight = [3, 4, 2, 2, 7]
numLayers = 5
domain = lambda i: (0,1)
rootState = frozenset([0,1,2,3,4])
def trFunc(s,d,i):
    if d == 1:
        if i in s:
            return s - {i} - set(neighbors[i])
        else:
            return None
    else:
        return s - {i}
costFunc = lambda s,d,i: d*weight[i]
isFeas = lambda s: s is not None
mergeFunc = lambda slist: frozenset(set.intersection(*(set(s) for s in slist)))
name = 'misp'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down compilation
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
mymdd.reduce_bottom_up(mergeFunc)

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Display MDD with GraphViz
mymdd.output_to_dot(nodeDotFunc=lambda s: '[label="{' + ', '.join(str(i) for i in s)  + '}"];', arcSortArgs=MDD._default_arcsortargs, nodeSortArgs=MDD._default_nodesortargs)
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
