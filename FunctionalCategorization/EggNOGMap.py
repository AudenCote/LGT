#!/usr/bin/env python3
# coding=utf-8

#Author: Auden Cote-L'Heureux
#Contact at acotelheureux@smith.edu, audenemil@gmail.com, or on Slack if a Katzlabber.
#Last updated 05/20/21
#Notes: should technically be able to run with python 2 or 3; it says it requires p4 but that's just for selecting sequences from trees -- you can probably modify this script to just use the EggNOG mapping bits.


import os
import sys
import re
from p4 import *
from Bio import SeqIO
import subprocess
from tqdm import tqdm


#Getting the arguments (the input directory of trees/seqs and whether or not you want to slim your GO terms)
def get_args():
	
	input_dir = ''
	slim = False
	
	if('--input_dir' in sys.argv):
		try:
			input_dir = sys.argv[sys.argv.index('--input_dir') + 1]
		except IndexError:
			bad_script_call()
	else:
		bad_script_call()
		
	if('--slim' in sys.argv):
		slim = True
		
	if(input_dir.endswith('/')):
		return input_dir[:-1], slim
	else:
		return input_dir, slim
		

#Utility function just for reading trees
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

	
#Utility function just for reading trees
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	tree_corrected.write(tree2correct)
	tree_corrected.close()	
	
	
#Utility function just for reading trees. See description
def root_by_major_clade(tree):

	###########################################################################
	# This function takes in a tree of interest (called by below functions)   #
	# and returns the tree rooted on the largest bacterial clade. If there is #
	# not a bacterial clade of significant size, it defaults to archaea, then #
	# Opisthokonta, and so on (see code for full hierarchy). This code was	  #
	# shamelessly stolen from Mario.										  #
	###########################################################################
	
	sizes_cladesBaZa = {}			
	sizes_cladesOp = {}				
	sizes_cladesPl = {}				
	sizes_cladesAm = {}				
	sizes_cladesEx = {}				
	sizes_cladesSr = {}				
			
	for node in tree.iterNodesNoRoot():				
		allTaxa_node = tree.getAllLeafNames(node)	
		
		MC_list =[]

		for taxon in allTaxa_node:					
			MC_list.append(taxon[:2])	
															
		if MC_list.count('Ba') + MC_list.count('Za') == (len(MC_list) or len(MC_list)-1):
			sizes_cladesBaZa[tree.node(node)] = len(allTaxa_node)	
		if MC_list.count('Op') == (len(MC_list) or len(MC_list)-1):
			sizes_cladesOp[tree.node(node)] = len(allTaxa_node)
		if MC_list.count('Pl') == (len(MC_list) or len(MC_list)-1):
			sizes_cladesPl[tree.node(node)] = len(allTaxa_node)
		if MC_list.count('Am') == (len(MC_list) or len(MC_list)-1):
			sizes_cladesAm[tree.node(node)] = len(allTaxa_node)
		if MC_list.count('Ex') == (len(MC_list) or len(MC_list)-1):
			sizes_cladesEx[tree.node(node)] = len(allTaxa_node)
		if MC_list.count('Sr') == (len(MC_list) or len(MC_list)-1):
			sizes_cladesSr[tree.node(node)] = len(allTaxa_node)
					
	if sizes_cladesBaZa: 		
		biggestBaZa = max(sizes_cladesBaZa, key=sizes_cladesBaZa.get) 
		tree.reRoot(biggestBaZa)
	else:
		if sizes_cladesOp:
			biggestOp = max(sizes_cladesOp, key=sizes_cladesOp.get)
			tree.reRoot(biggestOp)
		else:
			if sizes_cladesPl:
				biggestPl = max(sizes_cladesPl, key=sizes_cladesPl.get)
				tree.reRoot(biggestPl)
			else:
				if sizes_cladesAm:
					biggestAm = max(sizes_cladesAm, key=sizes_cladesAm.get)
					tree.reRoot(biggestAm)
				else:
					if sizes_cladesEx:
						biggestEx = max(sizes_cladesEx, key=sizes_cladesEx.get)
						tree.reRoot(biggestEx)
					else:
						if sizes_cladesSr:
							biggestSr = max(sizes_cladesSr, key=sizes_cladesSr.get)
							tree.reRoot(biggestSr)
											
	return tree
	

#Utility function just for writing trees
def write_nexus(path, colored_names, newick):
	
	with open(path, 'w') as o:
		ntax = str(len(colored_names))
			
		o.write('#NEXUS\n')	
		o.write('begin taxa;\n')
		o.write('\tdimensions ntax=' + ntax + ';\n')
		o.write('\ttaxlabels\n')
				
		for taxon in colored_names:
			o.write('\t' + taxon + '\n')
				
		o.write(';\nend;\n\n')
				
		o.write('begin trees;\n')
		o.write('\ttree tree_1 = [&R]\n')
		o.write(newick)
		o.write('end;\n\n')
				
		with open('figtree_format.txt', 'r') as ff:
			for line in ff:
				o.write(line)


#This function takes in a list of GO term identifiers and slims them using the generic GO-Slim database, then uniquifies
def go_slim_list(go_terms):

	slim_tuples = []; slim_list = [];
	for term in go_terms:
		try:
			res = subprocess.check_output('map_to_slim.py go-basic.obo goslim_generic.obo --term=' + term, shell = True)
			for line in res.split('\n'):
				if('GO:' in line):
					line = [term.strip() for term in line.split(';')]
					slim_tuples.append((term, line))
					slim_list.extend(line)
		except subprocess.CalledProcessError:
			continue
	
	return slim_tuples, list(dict.fromkeys(slim_list))
	
	
#I never really used this. For reducing GO term lists in a different way, using similarity & the tool ReviGO (recommended by Jaime & Carlos)
def revigo_reduce_list(go_terms):
	
	go_terms = '\n'.join(go_terms)
	os.system('/usr/local/bin/Rscript --vanilla query-revigo.r ' + go_terms)
	

#This function does various things related to grabbing sequences of interest, identified using trees. You probably won't use this function, and it requires p4 which is a PITA
def select_sequences(input_dir):

	for tree_file in os.listdir(input_dir):
		if(tree_file.endswith('.tre')):
			og = 'OG5_' + tree_file.split('OG5_')[1][:6]
			
			preguidance_file = ''
			for pg_file in os.listdir(input_dir):
				if('preguidance' in pg_file and og in pg_file):
					preguidance_file = input_dir + '/' + pg_file
			if(preguidance_file == ''):
				print('\nCould not find a pre-Guidance file for OG ' + og + '\n')
				bad_script_call()
			
			if(not os.path.isdir('../FunctionSearching/' + og)):
				os.system('mkdir ../FunctionSearching/' + og)
				
			os.system('cp ' + input_dir + '/' + tree_file + ' ../FunctionSearching/' + og + '/' + og + '.tre')
			
			fname = '../FunctionSearching/' + og + '/' + og + '.tre'
			var.trees = []
			nexus_to_newick('../FunctionSearching/' + og + '/' + og + '.tre')											
			correct_tree('../FunctionSearching/' + og + '/' + og + '.tre')														
			tree_file = '../FunctionSearching/' + og + '/' + og + '_temp.tre'
			read(tree_file) 
			tree = var.trees[0] 												
			
			tree = root_by_major_clade(tree)
			tree.writeNewick(fname)
				
			var.trees = []	
			read(tree_file) 
			tree = var.trees[0]
		
			os.remove(tree_file)
			
			continue
			
			seqs_already_mapped = []
			if(os.path.isfile('../FunctionSearching/' + og + '/eggnog_results.fa.emapper.annotations.tsv')):
				for line in open('../FunctionSearching/' + og + '/eggnog_results.fa.emapper.annotations.tsv'):
					if(line[0] != '#'):
						seqs_already_mapped.append(line.split('\t')[0].strip())
				os.system('mv ../FunctionSearching/' + og + '/eggnog_results.fa.emapper.annotations.tsv ../FunctionSearching/' + og + '/eggnog_results.fa.emapper.annotations.OLD.tsv')
			
			leaves = list(tree.getAllLeafNames(0))
			unfiltered_seqs = list(SeqIO.parse(preguidance_file, 'fasta'))
			seqs_in_tree = []
			for record in unfiltered_seqs:
				if(record.id in leaves and record.id not in seqs_already_mapped):
					seqs_in_tree.append(record)
					
			if(len(seqs_in_tree) == 0):
				os.system('mv ../FunctionSearching/' + og + '/eggnog_results.fa.emapper.annotations.OLD.tsv ../FunctionSearching/' + og + '/eggnog_results.fa.emapper.annotations.tsv')
			else:
				with open('../FunctionSearching/' + og + '/' + og + '_all_seqs.fasta', 'w') as o:
					for record in seqs_in_tree:
						o.write('>' + record.id + '\n' + str(record.seq) + '\n\n')

	
#This function is basically just a wrapper for querying the EggNOG databases using the EggNOG-mapper tool (needs to be installed) for each input file
def emap_selected_seqs():
	
	for dir in os.listdir('../FunctionSearching'):
		if('OG5_' in dir):
			already_done = False
			for file in os.listdir('../FunctionSearching/' + dir):
				if(file == 'eggnog_results.fa.emapper.annotations.tsv'):
					already_done = True
			if(not already_done):
				os.system('emapper.py -i ../FunctionSearching/' + dir + '/' + dir + '_all_seqs.fasta  --output eggnog_results.fa --output_dir ../FunctionSearching/' + dir + ' -m diamond')
			
				os.system('mv ../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.annotations ../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.annotations.tsv')
				#os.system('rm ../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.seed_orthologs')
			else:
				print(dir)
				
			if(os.path.isfile('../FunctionSearching/' + dir + '/eggnog_results.fa.emapper.annotations.OLD.tsv')):
				this_lines = []
				top_hashed_lines = []; bottom_hashed_lines = []
				for line in open('../FunctionSearching/' + dir + '/eggnog_results.fa.emapper.annotations.tsv'):
					if(line[0] == '#' and len(this_lines) == 0):
						top_hashed_lines.append(line)
					elif(line[1] == '#' and len(this_lines) > 0):
						bottom_hashed_lines.append(line)
					else:
						this_lines.append(line)
						
				for line in open('../FunctionSearching/' + dir + '/eggnog_results.fa.emapper.annotations.OLD.tsv'):
					if(line[0] != '#'):
						this_lines.append(line)
						
				with open('../FunctionSearching/' + dir + '/eggnog_results.fa.emapper.annotations.tsv', 'w') as o:
					for line in top_hashed_lines:
						o.write(line)
					for line in this_lines:
						o.write(line)
					for line in bottom_hashed_lines:
						o.write(line)
						
						
						
						
						