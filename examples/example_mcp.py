#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Maximum Cut Problem Example from Section 3.9 of Bergman et. al, Decision Diagrams for Optimization
numVrt = 4
wght = [[0, 1, 2, -2], [1, 0, 3, -1], [2, 3, 0, -1], [-2, -1, -1, 0]]
pos = lambda a: max(a,0)
neg = lambda a: min(a,0)

numLayers = numVrt
domain = lambda k: ['S'] if k == 0 else ['S', 'T']
rootState = tuple(0 for j in range(numVrt))
rootValue = sum(neg(wght[i][j]) for i in range(numVrt) for j in range(numVrt) if i < j)
def trFunc(s,d,k):
    if d == 'S':
        return tuple(s[l] + wght[k][l] if l > k else 0 for l in range(numVrt))
    else:
        return tuple(s[l] - wght[k][l] if l > k else 0 for l in range(numVrt))
def costFunc(s,d,k,ns):
    if k == 0:
        return rootValue
    else:
        if d == 'S':
            return pos(-s[k]) + sum(min(abs(s[l]), abs(wght[k][l])) for l in range(k+1,numVrt) if s[l]*wght[k][l] <= 0)
        else:
            return pos(s[k]) + sum(min(abs(s[l]), abs(wght[k][l])) for l in range(k+1,numVrt) if s[l]*wght[k][l] >= 0)
isFeas = lambda s,k: True
maxWidth = lambda k: 1 if k == 2 else 100
def mergeFunc(slist,k):
    newstate = []
    for l in range(numVrt):
        if l < k:
            newstate.append(0)
        elif all(u[l] >= 0 for u in slist):
            newstate.append(min( u[l] for u in slist ))
        elif all(u[l] <= 0 for u in slist):
            newstate.append(-min( abs(u[l]) for u in slist ))
        else:
            newstate.append(0)
    return tuple(newstate)
adjFunc = lambda w,os,ms,k: w + sum(abs(os[l]) - abs(ms[l]) for l in range(k,numVrt))
name = 'mcp'



# Construct the MDD
mymdd = MDD(name)
# Perform DP-based top-down exact compilation
#mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)
# Perform DP-based top-down relaxed compilation
mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas, maxWidth, mergeFunc, adjFunc)

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Output and display MDD with GraphViz
def adf(l,w,j):
    if l == 'S':
        return '[style=dotted,label="' + str(w) + '"];'
    else:
        return '[label="' + str(w) + '"];'
mymdd.output_to_dot(arcDotFunc=adf)
system('dot -Tps ' + name + '.gv -o ' + name + '.ps')
system('gv ' + name + '.ps')
