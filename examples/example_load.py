#!/usr/bin/python3
from argparse import ArgumentParser
from pymdd.mdd import *
from os import system

# Parser
parser = ArgumentParser(description='Loads and displays MDD from JSON file')
parser.add_argument('json_file', help='json file with MDD data')
args = parser.parse_args()

# Load MDD from file
json_file = args.json_file
mymdd = MDD()
mymdd.loadJSON(json_file)
# NOTE: for smisp.json, use the following instead
#from uuid import UUID
#mymdd.loadJSON(json_file, stateLoadFunc=lambda s: eval(s, globals()))

# Print the contents
print(mymdd)
print(mymdd.find_longest_path())

# Output and display MDD with GraphViz
mymdd.output_to_dot()
system('dot -Tps ' + mymdd.name + '.gv -o ' + mymdd.name + '.ps')
system('gv ' + mymdd.name + '.ps')
