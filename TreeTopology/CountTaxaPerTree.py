import os
import sys
import re
from p4 import *


var.doRepairDupedTaxonNames = 1


def help():
	
	print('Aaaah!')
	exit()


def get_args():

	if('--help' in sys.argv or '-h' in sys.argv): help()

	input_dir = ''

	if('--input_dir' in sys.argv or '-i' in sys.argv):
		try:
			if('--input_dir' in sys.argv):
				input_dir = sys.argv[sys.argv.index('--input_dir') + 1]
			else:
				input_dir = sys.argv[sys.argv.index('-i') + 1]
		except IndexError:
			help()
	else:
		help()
		
	return input_dir	
	

def nexus_to_newick(fname):

	nexus = False
	
	for line in open(fname):
		if('NEXUS' in line):
			nexus = True
			
	if(nexus):
		newick = ''
		for line in open(fname):
			line = line.split(' ')[-1]
			if(line.startswith('(') or line.startswith('tree1=')):
				newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')
				
		with open(fname, 'w') as o:
			o.write(newick)
			
			
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	tree_corrected.write(tree2correct)
	tree_corrected.close()


def taxa_by_tree(input_dir):

	counts_by_tree = { }
	
	for tre_file in os.listdir(input_dir):
		if(tre_file.endswith('.tre') or tre_file.endswith('.tree')):
			fname = input_dir + '/' + tre_file
		
			var.trees = []
			
			nexus_to_newick(fname)
																
			correct_tree(fname)								
														
			tree_file = fname.split('.tre')[0] + '_temp.tre'
			read(tree_file) 
			tree = var.trees[0] 						
															
			os.remove(tree_file)
			
			counts_by_tree.update({ 'OG5_' + fname.split('OG5_')[-1][:6] : { taxon : len([leaf for leaf in tree.getAllLeafNames(0) if leaf[:10] == taxon]) for taxon in [l[:10] for l in tree.getAllLeafNames(0)] } })
			
	all_taxa = list(dict.fromkeys([taxon for tree in counts_by_tree for taxon in counts_by_tree[tree]]))
	
	with open('taxon_counts_by_tree.csv', 'w') as o:
		o.write('Taxon,' + ','.join([og for og in counts_by_tree]) + '\n')
		for taxon in all_taxa:
			o.write(taxon)
			for og in counts_by_tree:
				if(taxon in counts_by_tree[og]):
					o.write(',' + str(counts_by_tree[og][taxon]))
				else:
					o.write(',0')
			o.write('\n')
			
	
def main():
	
	input_dir = get_args()

	taxa_by_tree(input_dir)
	
	
main()

