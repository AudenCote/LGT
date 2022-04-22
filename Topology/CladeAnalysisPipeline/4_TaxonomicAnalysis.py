############################################################
			
			 #Dependencies and global variables
			
############################################################


from p4 import *
import os, re
import sys
import csv
from Bio import SeqIO
import subprocess
import operator

var.doRepairDupedTaxonNames = 1


#If something goes wrong in parsing the arguments
def bad_script_call():
	
	print("\nPlease input the path to a spreadsheet with three columns. The first column have the ten digit codes of all the taxa that could be in the trees, the second the genus + species names, and the third should be the taxonomy as returned from NCBI or found in the TAXASELECTION/PhyloTOL master taxon list spreadsheet. If you don't know what this is, ask someone!\n\nExample: python 1_SupplementaryTable.py --taxonomy_file <path to spreadsheet>\n\nAlternatively, provide a similarly formatted file called 'taxonomy.csv' in the Scripts folder and this script should locate it.\n")
	
	exit()
	

def get_args():
	
	taxonomy_file = ''

	if('--taxonomy_file' in sys.argv):
		try:
			taxonomy_file = sys.argv[sys.argv.index('--taxonomy_file') + 1]
		except IndexError:
			bad_script_call()
	else:
		if(os.path.isfile('taxonomy.csv')):
			taxonomy_file = 'taxonomy.csv'
		else:
			bad_script_call()
			
	target_clades = []
	if('--target_clades' in sys.argv):
		i = sys.argv.index('--target_clades') + 1
		while '--' not in sys.argv[i]:
			target_clades.append(sys.argv[i].lower())
			if(i < len(sys.argv) - 1):
				i += 1
			else:
				break
			
	return [taxonomy_file, target_clades]

#Removing any dashes from the tree so that the tip names can be extracted using p4
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	tree_corrected.write(tree2correct)
	tree_corrected.close()
	
	
#Utility function
def taxa_retrieval_from_tree(fname):
	
	#Reading all of the tip names in the tree using p4
	var.trees = []			
													
	correct_tree(fname)																		
	tree_file = fname.split('.tre')[0] + '_temp.tre'
	read(tree_file) 
	tree = var.trees[0] 													
	os.remove(tree_file)
	
	return [taxon[:10] for taxon in list(tree.getAllLeafNames(0))]
	
	
def root_by_major_clade(tree):

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


############################################################
			
			   #Assessing taxonomic patterns
			
############################################################


def get_largest_euk_clade(tree, clade_id, nodes_to_exclude):

	best_node = None
	best_size = 0
	contam_threshold = .125
	absolute_n_threshold = 1
	
	forbidden_nodes = [node for node in nodes_to_exclude]
	for node in nodes_to_exclude:
		for num in tree.getNodeNumsAbove(node):
			forbidden_nodes.append(tree.node(num))
				
	for node in tree.iterNodesNoRoot():
		if(len(tree.getAllLeafNames(node)) > absolute_n_threshold and node not in forbidden_nodes):
			leaves = tree.getAllLeafNames(node)
		
			num = 0.0; dem = 0.0;
			for leaf in leaves:
				flag = False
				if(leaf[:2].lower() == clade_id):
					flag = True
						
				if(flag):
					dem += 1.0
				else:
					dem += 1.0; num += 1.0;
					
			if(dem == 0.0):
				ratio = 0.0
			else:
				ratio = num / dem
			
			if((ratio <= contam_threshold or (dem >= 3 and num <= 1)) and dem > best_size):
				best_node = node
				best_size = dem
						
	return best_node


#This function returns a spreadsheet listing the size of every clade of every given taxonomic grouping that meets given criteria	
def get_all_clade_sizes():
	
	all_clades = [ 'am', 'sr', 'op', 'pl', 'ex', 'ee', 'ba', 'za' ]
			
	size_by_clade_by_og = { }
	nodes_by_clade_by_og = { }
	leaf_names_by_clade_by_og = { }

	for f, file in enumerate(os.listdir('../Subtrees/Fulltrees')):
		if(file.endswith('.tre')):
			print(str(f) + '. ' + file)
			fname = '../Subtrees/Fulltrees/' + file
			og = 'OG5_' + fname.split('OG5_')[1][:6]
			
			size_by_clade_by_og.update({ og : {} })
		
			var.trees = []			
														
			correct_tree(fname)													
			tree_file = fname.split('.tre')[0] + '_temp.tre'
			read(tree_file) 
			try:
				tree = var.trees[0] 	
			except:
				continue					
												
			os.remove(tree_file)
			
		#	tree = root_by_major_clade(tree)
			
			#Storing the actual nodes for getting taxonomy later
			nodes_by_clade = { }
	
			for clade in all_clades:
				nodes_by_clade.update({ clade : [] })
				
				best_node = get_largest_euk_clade(tree, clade, nodes_by_clade[clade])
				if(best_node != None):
					nodes_by_clade[clade].append(best_node)
							
				while best_node != None:
					best_node = get_largest_euk_clade(tree, clade, nodes_by_clade[clade])
					nodes_by_clade[clade].append(best_node)
					
			nodes_by_clade_by_og.update({ og : nodes_by_clade })
			
			leaf_names_by_clade_by_og.update({ og : {} })
		
			for clade in nodes_by_clade:
				leaf_list = []
				for node in nodes_by_clade[clade]:
					if(node != None):
						leaf_list.extend([taxon[:10] for taxon in tree.getAllLeafNames(node)])
				leaf_names_by_clade_by_og[og].update({ clade : leaf_list })
		
			#Getting the number of children of each node of significant size of each major clade
			for clade in nodes_by_clade:
				for og2 in size_by_clade_by_og:
					if(clade not in size_by_clade_by_og[og2]):
						size_by_clade_by_og[og2].update({ clade : [] })
				for node in nodes_by_clade[clade]:
					if(node != None):
						leaves = []
						for leaf in tree.getAllLeafNames(node):
							if(leaf[:2].lower() == clade):
								leaves.append(leaf)
						leaves = list(dict.fromkeys([leaf[:10] for leaf in leaves]))
						size_by_clade_by_og[og][clade].append(len(leaves))
					
	return leaf_names_by_clade_by_og


#Getting the occurence of taxonomic levels in selected eukaryotic clades
def get_euk_clade_taxonomy_counts(taxonomy_file, leaves_by_major_clade, target_clades):

	taxonomy_counts = {}
	#For ordering the taxonomic levels later on
	positions = { }
	max_pos = 0
	
	for clade_tree in os.listdir('../Subtrees/Fulltrees'):
		if(clade_tree.endswith('.tre')):
			og = 'OG5_' + clade_tree.split('OG5_')[1][:6]
			taxonomy_counts.update({ og : {} })
	
	#Setting up a dictionary to hold taxon counts for each OG (not most efficient method but good for dev)
	for clade_tree in os.listdir('../Subtrees/Fulltrees'):
		if(clade_tree.endswith('.tre')):
			og = 'OG5_' + clade_tree.split('OG5_')[1][:6]
									
			#Get all taxa present in subtree
			all_taxa = taxa_retrieval_from_tree('../Subtrees/Fulltrees/' + clade_tree)
			
			rel_taxa = []
			for taxon in all_taxa:
				if(taxon in leaves_by_major_clade[og][taxon[:2].lower()] and taxon[:2].lower() in target_clades):
					rel_taxa.append(taxon)
					
			print(len(rel_taxa), list(dict.fromkeys([tax[:2] for tax in rel_taxa])))
					
			#Getting the taxonomy counts for the relevant euk clade
			for t, taxon in enumerate(rel_taxa):
				for line in open(taxonomy_file):
					if(taxon in line and taxon[:2].lower() in target_clades):
						#Retrieving and formatting the taxonomy
						line = line.split(',')
						tax_list = line[2].strip().split(';')
						
						idx_list = [idx + 1 for idx, val in
									enumerate(tax_list) if 'cellular' in val]
									
						if(len(idx_list) > 0):
							tax_list = [tax_list[i: j] for i, j in
										zip([0] + idx_list, idx_list + 
										([len(tax_list)] if idx_list[-1] != len(tax_list) else []))][0]
										
						if(line[0][0:2] == 'Ex'):
							tax_list = ['eukarya', 'excavata'] + tax_list
						elif(line[0][0:2] == 'Sr'):
							tax_list = ['eukarya', 'sar'] + tax_list
						elif(line[0][0:2] == 'Pl'):
							tax_list = ['eukarya', 'archaeplastida'] + tax_list
						else:
							tax_list = ['eukarya'] + tax_list
							
						print(tax_list)
						print('')
										
						for tax_level in tax_list[::-1]:
							if('cellular' in tax_level.lower() or 'eukaryot' in tax_level.lower() or tax_level.lower() == '' or 'incertae sedis' in tax_level.lower()):
								tax_list.remove(tax_level)
														
						#Appending the taxonomy to the master dictionary
						for l, tax_level in enumerate(tax_list):
							tax_level = tax_level.strip().lower()
							
							if(l not in positions):
								flag = 0
								for pos in positions:
									if(tax_level in positions[pos]):
										flag = 1
								
								if(flag == 0):
									positions.update({ l : { tax_level : 0 }})
							else:
								flag = 0
								for pos in positions:
									if(tax_level in positions[pos]):
										flag = 1
										
								if(flag == 0):
									positions[l].update({ tax_level : 0 })
								
							if(tax_level not in taxonomy_counts[og]):
								for tree in taxonomy_counts:
									taxonomy_counts[tree].update({ tax_level : 0 })
								
							taxonomy_counts[og][tax_level] += 1
															
	#Sorting taxonomic levels by the sum of their counts across all OGs
	for og in taxonomy_counts:
		for tax_level in taxonomy_counts[og]:
			for pos in positions:
				if(tax_level in positions[pos]):
					positions[pos][tax_level] += taxonomy_counts[og][tax_level]
	
	for pos in positions:
		sorted_tuples = sorted(positions[pos].items(), key=lambda item: item[1], reverse = True)
		positions[pos] = sorted_tuples
		
	print(positions[0])
								
	#Writing the output
	with open('../Spreadsheets/euk_clade_taxonomy_counts.csv', 'w') as o:
		o.write('Gene Family,')
		for pos in positions:
			for tax in positions[pos]:
				o.write(tax[0] + ',')
		o.write('\n')
		
		for og in taxonomy_counts:
			o.write(og + ',')
			for pos in positions:
				for tax in positions[pos]:
					o.write(str(taxonomy_counts[og][tax[0]]) + ',')
			o.write('\n')
			
	return positions
			
	
#Function for hypothesizing the depth of the XGT event		
def get_depth_of_transfer(positions):
	
	taxonomy_counts = { }
	clades = []
	
	#Holding some initial data (as created by function get_euk_clade_taxonomy_counts()
	for line in open('../Spreadsheets/euk_clade_taxonomy_counts.csv'):
		if('Gene Family' in line):
			line = line.split(',')[1:]
			for clade in line:
				if(clade != '\n'):
					clades.append(clade.strip())
		else:
			line2 = line.split(',')[1:]; og = line.split(',')[0]
			taxonomy_counts.update({ og : [] })
			for cnt in line2:
				if(cnt != '\n'):
					taxonomy_counts[og].append(cnt.strip())
					
	final_dict = { }
	
	#Looking at the taxonomy of each OG	
	for og in taxonomy_counts:
		#Total taxa in major clade of interest
		rel_tot = int(taxonomy_counts[og][0])
		for c, cnt in enumerate(taxonomy_counts[og]):
			if(int(cnt) != 0 and int(cnt) < rel_tot):
				#Getting the counts of the occurrence of each taxonomic level
				level_counts = { }
				for pos in positions:
					pos_clades = [tax[0] for tax in positions[pos]]
					if(clades[c] in pos_clades):
						for pos_clade in pos_clades:
							try:
								level_counts.update({ pos_clade : int(taxonomy_counts[og][clades.index(pos_clade)]) })
							except ValueError:
								print('\nCould not find taxonomic level ' + pos_clade + ' for OG ' + og + '\n')
					
				#First hypothesis of depth clade
				for v, val in enumerate(taxonomy_counts[og][:c][::-1]):
					if(int(val) == rel_tot):
						depth_clade = clades[c - (v + 1)]
						break
						
				#Update hypothesis of depth clade using (Bayesian) probability of > .98 of the taxon appearing
				total = 0
				for clade in level_counts:
					total += level_counts[clade] + 1
					
				for clade in level_counts:
					if((level_counts[clade] + 1)/total > .98):
						depth_clade = clade
						
				final_dict.update({ og : [level_counts, depth_clade] })						
				break
				
	#Writing output
	with open('../Spreadsheets/depth_hypotheses.csv', 'w') as o:
		for og in final_dict:
			o.write('Gene Family,Hypothesis,')
			for clade in final_dict[og][0]:
				if(final_dict[og][0][clade] > 0):
					o.write(clade + ',')
			o.write('\n' + og + ',' + final_dict[og][1] + ',')
			for clade in final_dict[og][0]:
				if(final_dict[og][0][clade] > 0):
					o.write(str(final_dict[og][0][clade]) + ',')
			o.write('\n\n')
			
			
def counts_by_level_spreadsheet_for_figure():
	try:
		hypothesis_sheet = open('../Spreadsheets/depth_hypotheses.csv', 'r')
	except:
		print('\nThe spreadsheet containing data on the taxonomic depths of the putative XGT in each tree, depth_hypotheses.csv could not be found. Please make sure that this spreadsheet has been generated before running the function counts_by_level_spreadsheet_for_figure()\n')
		exit()
		
	depth_counts = { }
	with hypothesis_sheet as hts_open:
		for line in hts_open:
			if('Hypothesis' not in line and len(line.split(',')) > 2):
				level = line.split(',')[1].strip()
				if(level not in depth_counts):
					depth_counts.update({ level : 1 })
				else:
					depth_counts[level] += 1
				
	with open('../Spreadsheets/depth_counts_for_figure.csv', 'w') as o:
		o.write('Taxonomic Level, Count\n')
		for level in depth_counts:
			o.write(level + ',' + str(depth_counts[level]) + '\n')
				
			
def main():
		
	taxonomy_file, target_clades = get_args()
	
	leaves_by_major_clade = get_all_clade_sizes()

	positions = get_euk_clade_taxonomy_counts(taxonomy_file, leaves_by_major_clade, target_clades)
	
	get_depth_of_transfer(positions)
	
	counts_by_level_spreadsheet_for_figure()
			
			
main()	
			
		