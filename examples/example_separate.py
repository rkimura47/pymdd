#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Compilation by Separation Example from Section 3.11 of Bergman et. al, Decision Diagrams for Optimization
G = [[(1,2), (1,3)], [(2,1), (2,3)], [(3,1), (3,2)]]
neighbors = [ [v for (u,v) in G[j]] for j in range(3) ]
weight = [1, 1, 1]
numLayers = 3
domain = lambda j: (0,1)

rootState = frozenset([1,2,3])
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
maxWidth = lambda j: 100
name = 'smisp'

def counter(k = [-1]):
    k[0] += 1
    return k[0]
# This is admittedly a bit of a hack... a more 'appropriate' option would be
#from uuid import uuid1
#counter = uuid1



# Construct the MDD
mymdd = MDD(name)

mymdd.compile_trivial(numLayers, domain, lambda d,j: costFunc(None,d,j,None), lambda j: counter())
mymdd.filter_and_refine_constraint(trFunc, rootState, isFeas, lambda s,j: counter(), maxWidth)
mymdd.reduce_bottom_up(lambda slist,j: slist[0])

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Output and display MDD with GraphViz
ndf = lambda s,j: 'label=""'
mymdd.output_to_dot(nodeDotFunc=ndf, arcSortArgs=MDD._default_asa)
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
