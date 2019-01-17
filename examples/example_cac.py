#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Canonical Art Costs Example from Hooker, "Decision Diagrams and Dynamic Programming"
feasSolutions = [(0,1,0,1), (0,1,1,0), (0,1,1,1), (1,0,1,1), (1,1,0,0), (1,1,0,1), (1,1,1,0), (1,1,1,1)]
# Non-separable cost
cost = [6, 7, 8, 5, 6, 8, 7, 9]
# Separable cost
#cost = [sum(map(lambda i,j: i*j, (3,5,4,6), s)) for s in feasSolutions]

# MDD Parameters
mergeFunc = lambda slist,i: min(slist)
name = 'cac'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down compilation
pathList = list(zip(cost, feasSolutions))
mymdd.compile_pathlist(pathList)
mymdd.reduce_bottom_up(mergeFunc)

# Print the contents
print(mymdd)

# Output and display MDD with GraphViz
ndf = lambda s,j: 'label=""'
def adf(label, weight, layer): 
    lab = []
    if label == 0:
        lab.append('style=dotted')
    if weight != 0:
        lab.append('label="%d"' % weight)
    return ','.join(lab)
mymdd.output_to_dot(ndf, adf)
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
