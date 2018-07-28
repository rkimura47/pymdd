#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Knapsack Example from Section 3.2 of Bergman et al, Decision Diagrams for Optimization
numLayers = 4
capacity = 6
profit = [8, 7, 6, 14]
weight = [3, 3, 4, 6]

domain = lambda i: (0,1)
rootState = capacity
trFunc = lambda s,d,i: s-d*weight[i]
costFunc = lambda s,d,i: d*profit[i]
isFeas = lambda s,i: s >= 0
mergeFunc = lambda slist,i: min(slist)
name = 'knapsack'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down compilation
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
mymdd.reduce_bottom_up(mergeFunc)

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Display MDD with GraphViz
mymdd.output_to_dot()
#mymdd.output_to_dot(arcSortArgs=MDD._default_arcsortargs, nodeSortArgs=MDD._default_nodesortargs)
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
