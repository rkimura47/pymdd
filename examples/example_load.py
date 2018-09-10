#!/usr/bin/python3
from pymdd.mdd import *
from os import system

# Load MDD from file
jsondir = './pymdd/examples/json/'
jsonfiles = ['knapsack.json', 'mcp.json', 'misp.json', 'smisp.json']
selectfile = jsonfiles[2]

mymdd = MDD()

if selectfile == 'smisp.json':
    from uuid import UUID
    mymdd.loadJSON(jsondir + selectfile, stateLoadFunc=lambda s: eval(s, globals()))
else:
    mymdd.loadJSON(jsondir + selectfile)

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Output and display MDD with GraphViz
mymdd.output_to_dot()
system('dot -Tps ' + mymdd.name + '.gv -o ' + mymdd.name + '.ps')
system('gv ' + mymdd.name + '.ps')
