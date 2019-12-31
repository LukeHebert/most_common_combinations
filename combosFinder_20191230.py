'''
author: Luke Hebert

date begun: 10/01/2018

powerset_list() definition adapted from:   https://docs.python.org/3/library/itertools.html     retrieved on 10/01/2018
'''

#need chain & combination for the powerset_list function definition
from itertools import chain, combinations
from os import name as OSname
import sys


#powerset_list takes an iterable and finds all possible nonredundant combinations of items in that list (each combination is internally nonredundant too)
def powerset_list(iterable):
    '''#example
       powerset_list([1,2,3]) --> [(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]'''
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(1,5+1)))

slash = '\\' if OSname == 'nt' else '/'

inFile_path = sys.argv[1]

#the dictionary of patients and corresponding lists of their variant-containing genes "vars"
#dictionary format/template:
'''
pt_vars = {'pt1':('A','B','D','E','P'),
         'pt2':('B','D','Q','O','A','E','P'),
         'pt3':('B','C','R'),
         'pt4':('D','L','O','A','E','P'),
         'pt5':('B','C','O')}
'''
pt_vars = {}
with open(inFile_path, 'r') as inFile:
	for i, line in enumerate(inFile):
		line_list = line.strip('\n').strip('\r').split('\t')
		if i == 0:
			genes_index = line_list.index("wnt&shroom_list")
			id_index = line_list.index("#ID")
		else:
			pt_id = line_list[id_index]
			varsList = line_list[genes_index].split(',')
			pt_vars[pt_id] = set(varsList)
'''
#print for analysis
for pt, vars in pt_vars.iteritems():
    print(pt + ': ' + str(vars))
print('\n')    
'''

#make a list with the format [['pt1', [all possible letter combos for pt1 as tuples]], ['pt2', [all possible letter combos for pt2 as tuples]]]
pt_var_combos = {}
for pt, vars in pt_vars.iteritems():
    pt_var_combos[pt] = powerset_list(vars)
'''
#print for analysis
for pt, combos in pt_var_combos.iteritems():
    print(pt + ': ' + str(combos) + '\n')
'''

#make a dictionary with the format {[(combo1): 'pt1', 'pt4', 'pt5', (combo2): 'pt4', (combo3): 'pt2'}
var_combos_pt = {}
for pt, combos in pt_var_combos.iteritems():
    for var_combo in combos:#iterate through every variant combination list 'pt[1]' for every patient item 'pt' in nested list 'pt_var_combos'
        if var_combo not in var_combos_pt.keys():#if the combinations-to-patients-they're-in dictionary doesn't contain this variant combination, assign it
            var_combos_pt[var_combo] = [pt]
        elif var_combo in var_combos_pt.keys():
        	var_combos_pt[var_combo].append(pt)
'''
#print for analysis
for key, val in var_combos_pt.iteritems():
    print(str(key) + ': ' + str(val))
'''

#convert dictionary to a list of the format 
'''[
[(combo1), number_of_pts_with_combo1, [list of pts with combo1]], 
[(combo2), number_of_pts_with_combo2, [list of pts with combo2]]
]'''

var_combos_pt_list = []
for key, val in var_combos_pt.iteritems():
    var_combos_pt_list.append([key, len(val), val])

print('\n')

#sort that list by the number of patients each combo is in, so that most common combos are listed first and least common last
var_combos_pt_list.sort(key = lambda i: i[1], reverse = True)
#print for analysis
combo_count = 0

#output results
outFile_path = inFile_path.replace('.txt', '_combos.txt')
with open(outFile_path, 'w') as outFile:
	outFile.write('\t'.join(['gene combo count', 'gene combo', 'subjects affected count', 'subjects affected IDs', '\n']))
	for combo in var_combos_pt_list:
		combo_count += 1
		numCombo_str = str(len(list(combo[0])))
		combo_str = ','.join(list(combo[0]))
		numPts_str = str(combo[1])
		ptsList_str = ','.join(combo[2])
		outFile.write('\t'.join([numCombo_str, combo_str, numPts_str, ptsList_str, '\n']))
print('Number of combinations found: ' + str(combo_count))

