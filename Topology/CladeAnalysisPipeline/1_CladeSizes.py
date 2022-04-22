############################################################
			
			 #Dependencies and global variables
			
############################################################


from p4 import *
import os, re
import sys
import csv
from Bio import SeqIO
import subprocess
from tree_singleton_identification import get_singletons
import operator
from get_clade_totals import get_clade_totals

var.doRepairDupedTaxonNames = 1

euk_major_clades = [ 'am', 'sr', 'op', 'pl', 'ex', 'ee' ]

plastid_ogs = ['OG5_172393', 'OG5_146609', 'OG5_172394', 'OG5_134778', 'OG5_157649', 'OG5_136935', 'OG5_157648', 'OG5_132296', 'OG5_130791', 'OG5_147313', 'OG5_143116', 'OG5_145081', 'OG5_131049', 'OG5_138787', 'OG5_175720', 'OG5_130618', 'OG5_147239', 'OG5_140077', 'OG5_138109', 'OG5_180830', 'OG5_181659', 'OG5_144993', 'OG5_143146', 'OG5_150838', 'OG5_144937', 'OG5_134854', 'OG5_147208', 'OG5_150425']

############################################################
			
			  	#Various utility functions
			
############################################################


#Parse the arguments
def get_args():
	
	prokaryotes_only = False; include_prokaryotes = False; include_bacterial_bin = False;
	if('--prokaryotes_only' in sys.argv): prokaryotes_only = True
	if('--include_prokaryotes' in sys.argv): include_prokaryotes = True
				
	return [prokaryotes_only, include_prokaryotes]
	
	
def is_bacbin(tdc):

	if((tdc[:2] != 'Ba' and tdc[:2] != 'Za' and tdc[4] != 'b') or (tdc[:2] == 'Ba' or tdc[:2] == 'Za')):
		return False
	else:
		return True
		
		
#Removing any dashes from the tree so that the tip names can be extracted using p4
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	tree_corrected.write(tree2correct)
	tree_corrected.close()
	
	
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
			
	
#Utility function
def taxa_retrieval_from_tree(fname):
	
	#Reading all of the tip names in the tree using p4
	var.trees = []			
	
	nexus_to_newick(fname)						
	correct_tree(fname)																		
	tree_file = fname.split('.tre')[0] + '_temp.tre'
	read(tree_file) 
	tree = var.trees[0] 													
	os.remove(tree_file)
	
	return [taxon[:10] for taxon in list(tree.getAllLeafNames(0))]
	
	
############################################################
			
			 #Assessing eukaryotic diversity
			
############################################################


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

		
#Universal function for aid in the below functions
def assess_diversity(dir):
	
	diversity_by_og = { }
	
	for tree in os.listdir(dir):
		if(tree.endswith('.tre')):
			og = 'OG5_' + tree.split('OG5_')[1][:6]
			diversity_by_og.update({ og : { } })
	
	for tree in os.listdir(dir):
		if(tree.endswith('.tre')):
			og = 'OG5_' + tree.split('OG5_')[1][:6]
	
			#Get all taxa present in subtree
			taxa = taxa_retrieval_from_tree(dir + '/' + tree)
			
			seen_species = []
			
			for taxon in taxa:
				if(taxon[:2].lower() in euk_major_clades and taxon.lower() not in seen_species):
					if(taxon[:5].lower() not in diversity_by_og[og]):
						for all in diversity_by_og:
							diversity_by_og[all].update({ taxon[:5].lower() : 0 })
						
					diversity_by_og[og][taxon[:5].lower()] += 1
					
					seen_species.append(taxon.lower())
											
	return diversity_by_og
					

#Assessing the eukaryotic minor clades in the largest relevant eukaryotic clade
def diversity_in_euk_clade():

	diversity_by_og = assess_diversity('../Subtrees/EukClades')
		
	with open('../Spreadsheets/euk_clade_diversity.csv', 'w') as o:
		o.write('Gene Family,')
		for clade in diversity_by_og[next(iter(diversity_by_og))]:
			o.write(clade + ',')
		o.write('\n')
		
		for og in diversity_by_og:
			o.write(og + ',')
			for clade in diversity_by_og[og]:
				o.write(str(diversity_by_og[og][clade]) + ',')
			o.write('\n')
			
	
#Assessing the eukaryotic minor clades present in the entire original tree
def diversity_in_full_tree():
	
	diversity_by_og = assess_diversity('../Subtrees/Fulltrees')
	
	with open('../Spreadsheets/fulltree_diversity.csv', 'w') as o:
		o.write('Gene Family,')
		for clade in diversity_by_og[next(iter(diversity_by_og))]:
			o.write(clade + ',')
		o.write('\n')
		
		for og in diversity_by_og:
			o.write(og + ',')
			for clade in diversity_by_og[og]:
				o.write(str(diversity_by_og[og][clade]) + ',')
			o.write('\n')
			
			
def get_largest_euk_clade(tree, clade_id, nodes_to_exclude, minor):

	nbu_minor_clade_counts, bb_minor_clade_counts, combined_minor_clade_counts = get_clade_totals()

	dig8 = False
	if(len(clade_id) == 7 or len(clade_id) == 8):
		dig8 = True

	best_node = None
	best_size = 0
	contam_threshold = .125
	#Setting this to zero because filtering occurs after the data has been collected; if you would like to filter for best node size raise this value (e.g. absolute_n_threshold = 3)
	absolute_n_threshold = 0
	bacbin_conc_threshold = .33
	
	forbidden_nodes = [node for node in nodes_to_exclude]
	for node in nodes_to_exclude:
		for num in tree.getNodeNumsAbove(node):
			forbidden_nodes.append(tree.node(num))
				
	for node in tree.iterNodesNoRoot():
		if(len(tree.getAllLeafNames(node)) > absolute_n_threshold and node not in forbidden_nodes):
			leaves = tree.getAllLeafNames(node)
		
			num = 0.0; dem = 0.0;
			
			#Default clade-determining mechanism
			for leaf in leaves:
				flag = False
				if(dig8):
					if((leaf[:4] + leaf[5:8]).lower() == clade_id or leaf[:4].lower() == clade_id[:4]):
						flag = True
				if(leaf[:5].lower() == clade_id or leaf[:8].lower() == clade_id or leaf[:2].lower() == clade_id or leaf[:4].lower() == clade_id):
					flag = True
						
				if(flag):
					dem += 1.0
				else:
					dem += 1.0; num += 1.0;
					
			if(dem == 0.0):
				ratio = 0.0
			else:
				ratio = num / dem
				
			#If the taxa of the relevant minor clade (clade_id) are all bacterial bin, then only count the clade if a certain proportion of the taxa in the minor clade are represented in the bacterial bin clade
			bb_leaves = []
			good_bb_conc = True
			if(minor and (dig8 or len(clade_id) == 4) and clade_id[:2].lower() != 'ba' and clade_id[:2].lower() != 'za'):
				if(dig8):
					if(len(clade_id) == 7):
						bb_leaves = [leaf for leaf in leaves if((leaf[:4] + leaf[5:8]).lower() == clade_id and leaf[4] == 'b')]
						code_for_count = clade_id[:4] + 'b_' + clade_id[5:]
					else:
						bb_leaves = [leaf for leaf in leaves if((leaf[:8]).lower() == clade_id and leaf[4] == 'b')]
						code_for_count = clade_id
				else:
					bb_leaves = [leaf for leaf in leaves if(leaf[:4].lower() == clade_id and leaf[4] == 'b')]
					code_for_count = clade_id[:4] + 'b'
			elif(minor and (clade_id[:2] != 'ba' and clade_id[:2] != 'za' and clade_id[4] == 'b')):
				bb_leaves = [leaf for leaf in leaves if(leaf[:5].lower() == clade_id and leaf[4] == 'b')]
				code_for_count = clade_id
			
			if(clade_id[:2].lower() != 'ba' and clade_id[:2].lower() != 'za'):
				if(len(bb_leaves) == dem - num and (len(bb_leaves) < bacbin_conc_threshold * bb_minor_clade_counts[code_for_count] or len(bb_leaves) < 4)):
					good_bb_conc = False
			
			#rel_leaves = [leaf for leaf in leaves if leaf[:4] == 'Am_t' and leaf[6:8] == 'He']
			#if(len(leaves) - len(rel_leaves) < 10 and len(rel_leaves) > 10):
				#print(leaves)
			if((ratio <= contam_threshold or (dem >= 3 and num <= 1)) and dem > best_size and good_bb_conc):
				best_node = node
				best_size = dem
						
	return best_node
			
		
#This function returns a spreadsheet listing the size of every clade of every given taxonomic grouping that meets given criteria	
def get_all_clade_sizes(prokaryotes_only, include_prokaryotes, minor, bacbin, combine, select_taxa):
	
	if(prokaryotes_only):
		if(minor):
			all_clades = ['Ba_ac', 'Ba_ad', 'Ba_aq', 'Ba_ar', 'Ba_ba', 'Ba_bc', 'Ba_bi', 'Ba_ca', 'Ba_cd', 'Ba_ch', 'Ba_cr', 'Ba_cv', 'Ba_cy', 'Ba_de', 'Ba_df', 'Ba_di', 'Ba_el', 'Ba_fb', 'Ba_fc', 'Ba_fn', 'Ba_fu', 'Ba_ge', 'Ba_is', 'Ba_me', 'Ba_ni', 'Ba_pa', 'Ba_pb', 'Ba_pd', 'Ba_pg', 'Ba_pl', 'Ba_sp', 'Ba_sy', 'Ba_te', 'Ba_th', 'Ba_ts', 'Ba_le', 'Ba_cp', 'Ba_fe', 'Ba_nt', 'Ba_cn', 'Ba_ft', 'Za_cr', 'Za_ea', 'Za_eb', 'Za_ec', 'Za_eh', 'Za_em', 'Za_ep', 'Za_et', 'Za_eu', 'Za_ey', 'Za_ko', 'Za_na', 'Za_nh', 'Za_pa', 'Za_th', 'Za_lo', 'Za_as', 'Za_ba', 'Za_ba', 'Za_al']
		else:
			all_clades = [ 'ba', 'za' ]
	elif(include_prokaryotes):
		if(minor):
			all_clades = ['Am_ar', 'Am_di', 'Am_my', 'Am_hi', 'Am_is', 'Am_th', 'Am_tu', 'Am_va', 'Am_uk', 'EE_ap', 'EE_br', 'EE_cr', 'EE_ha', 'EE_is', 'EE_ka', 'EE_ce', 'EE_uk', 'Ex_eu', 'Ex_fo', 'Ex_he', 'Ex_is', 'Ex_ja', 'Ex_ma', 'Ex_ox', 'Ex_pa', 'Op_ch', 'Op_fu', 'Op_ic', 'Op_is', 'Op_me', 'Op_nu', 'Pl_gr', 'Pl_gl', 'Pl_rh', 'Sr_ap', 'Sr_ch', 'Sr_ci', 'Sr_di', 'Sr_is', 'Sr_pe', 'Sr_rh', 'Sr_st', 'Ba_ac', 'Ba_ad', 'Ba_aq', 'Ba_ar', 'Ba_ba', 'Ba_bc', 'Ba_bi', 'Ba_ca', 'Ba_cd', 'Ba_ch', 'Ba_cr', 'Ba_cv', 'Ba_cy', 'Ba_de', 'Ba_df', 'Ba_di', 'Ba_el', 'Ba_fb', 'Ba_fc', 'Ba_fn', 'Ba_fu', 'Ba_ge', 'Ba_is', 'Ba_me', 'Ba_ni', 'Ba_pa', 'Ba_pb', 'Ba_pd', 'Ba_pg', 'Ba_pl', 'Ba_sp', 'Ba_sy', 'Ba_te', 'Ba_th', 'Ba_ts', 'Ba_le', 'Ba_cp', 'Ba_fe', 'Ba_nt', 'Ba_cn', 'Ba_ft', 'Za_cr', 'Za_ea', 'Za_eb', 'Za_ec', 'Za_eh', 'Za_em', 'Za_ep', 'Za_et', 'Za_eu', 'Za_ey', 'Za_ko', 'Za_na', 'Za_nh', 'Za_pa', 'Za_th', 'Za_lo', 'Za_as', 'Za_ba', 'Za_ba', 'Za_al', 'Am_tu_He', 'Am_tu_Hp', 'Sr_ci_Cu', 'Sr_ci_Lx']
		else:
			all_clades = [ 'am', 'sr', 'op', 'pl', 'ex', 'ee', 'ba', 'za' ]
	else:
		if(minor):
			all_clades = ['Am_ar', 'Am_di', 'Am_my', 'Am_hi', 'Am_is', 'Am_th', 'Am_tu', 'Am_va', 'Am_uk', 'EE_ap', 'EE_br', 'EE_cr', 'EE_ha', 'EE_is', 'EE_ka', 'EE_ce', 'EE_uk', 'Ex_eu', 'Ex_fo', 'Ex_he', 'Ex_is', 'Ex_ja', 'Ex_ma', 'Ex_ox', 'Ex_pa', 'Op_ch', 'Op_fu', 'Op_ic', 'Op_is', 'Op_me', 'Op_nu', 'Pl_gr', 'Pl_gl', 'Pl_rh', 'Sr_ap', 'Sr_ch', 'Sr_ci', 'Sr_di', 'Sr_is', 'Sr_pe', 'Sr_rh', 'Sr_st', 'Am_tu_He', 'Am_tu_Hp', 'Sr_ci_Cu', 'Sr_ci_Lx']
		else:
			all_clades = [ 'am', 'sr', 'op', 'pl', 'ex', 'ee' ]
			
	if(bacbin and minor):
		temp = []
		for clade in all_clades:
			if(clade[:2] == 'ba' or clade[:2] == 'za'):
				temp.append(clade)
			else:
				if(len(clade) == 5):
					temp.append(clade[:4] + 'b')
				elif(len(clade) == 8):
					temp.append(clade[:4] + 'b' + clade[5:])
		all_clades = temp
	elif(combine and minor):
		temp = []
		for clade in all_clades:
			if(clade[:2] != 'ba' and clade[:2] != 'za'):
				if(len(clade) == 5):
					temp.append(clade[:4])
				elif(len(clade) == 8):
					temp.append(clade[:4] + clade[5:])
		all_clades = temp
		
	bad_clades = ['am_th', 'ee_ce', 'op_is', 'sr_ch', 'am_hi', 'ee_ka', 'ex_is', 'sr_is', 'am_hb', 'ee_kb', 'ex_ib', 'sr_ib', 'am_va', 'am_vb', 'op_nu', 'op_nb', 'am_h', 'ee_k', 'ex_i', 'sr_i', 'am_v', 'op_n', 'am_uk', 'am_ub', 'am_u', 'ee_uk', 'ee_ub', 'ee_u', 'ex_ox', 'ex_ob', 'ex_o']
			
	all_clades = [code.lower() for code in all_clades if(code.lower() not in bad_clades)]
	#all_clades = ['am_t_he']
	
	bad_singletons = [line.strip().lower() for line in open('bad_singletons.txt')]
		
	size_by_clade_by_og = { }
	nodes_by_clade_by_og = { }
	leaf_names_by_clade_by_og = { }
	
	leaves_to_keep_by_og = { }
	leaves_to_remove_by_og = { }
	
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
			
			tree = root_by_major_clade(tree)
			
			all_leaves_in_tree = []
			for node in tree.iterNodesNoRoot():				
				all_leaves_in_tree += list(tree.getAllLeafNames(node))
				
			all_leaves_in_tree = list(dict.fromkeys(all_leaves_in_tree))
			leaves_to_keep = []
			leaves_to_remove = []

			#Storing the actual nodes for getting taxonomy later
			nodes_by_clade = { }
	
			for clade in all_clades:
				nodes_by_clade.update({ clade : [] })
				
				best_node = get_largest_euk_clade(tree, clade, nodes_by_clade[clade], minor)
				if(best_node != None):
					nodes_by_clade[clade].append(best_node)
							
				while best_node != None:
					best_node = get_largest_euk_clade(tree, clade, nodes_by_clade[clade], minor)
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
				dig8 = False
				if(len(clade) == 7 or len(clade) == 8):
					dig8 = True
				for og2 in size_by_clade_by_og:
					if(clade not in size_by_clade_by_og[og2]):
						size_by_clade_by_og[og2].update({ clade : [] })
				for node in nodes_by_clade[clade]:
					if(node != None):
						leaves = []
						for leaf in tree.getAllLeafNames(node):
							if(dig8):
								if((leaf[:4] + leaf[5:8]).lower() == clade):
									leaves.append(leaf)
							if(leaf[:2].lower() == clade or leaf[:5].lower() == clade or leaf[:8].lower() == clade or leaf[:4].lower() == clade):
								leaves.append(leaf)
						l_wh = list(dict.fromkeys([leaf for leaf in leaves]))
						leaves = list(dict.fromkeys([leaf[:10] for leaf in leaves]))
						
						if(combine and select_taxa):
							if(len(leaves) == 2):
								if(leaves[0] not in bad_singletons or leaves[1] not in bad_singletons):
									leaves_to_keep += l_wh
							else:
								#if(clade == 'am_t_he'):
								#	print(l_wh)
								leaves_to_keep += l_wh
						if(len(leaves) > 0):
							size_by_clade_by_og[og][clade].append(len(leaves))
						
			#Getting a record of singletons in the tree (turn this off if you dont want any ones showing up in the output spreadsheet)
			for singleton in get_singletons(fname):
				if(singleton[:8].lower() in all_clades):
					size_by_clade_by_og[og][singleton[:8].lower()].append(' 1')
				if(singleton[:5].lower() in all_clades):
					size_by_clade_by_og[og][singleton[:5].lower()].append(' 1')
				if(singleton[:2].lower() in all_clades):
					size_by_clade_by_og[og][singleton[:2].lower()].append(' 1')
				if(singleton[:4].lower() in all_clades):
					size_by_clade_by_og[og][singleton[:4].lower()].append(' 1')
				if((singleton[:4] + singleton[5:8]).lower() in all_clades):
					size_by_clade_by_og[og][(singleton[:4] + singleton[5:8]).lower()].append(' 1')
					
			for leaf in all_leaves_in_tree:
				if(leaf not in leaves_to_keep and leaf[:10].lower() not in bad_singletons and not is_bacbin(leaf[:10])):
					leaves_to_keep.append(leaf)
				if(leaf not in leaves_to_keep):
					leaves_to_remove.append(leaf)
					
			#for leaf in leaves_to_keep:
				#if('am_t' in leaf.lower()):
				#	print(leaf)
					
			leaves_to_keep_by_og.update({ og : leaves_to_keep })
			leaves_to_remove_by_og.update({ og : leaves_to_remove })
				
	#Writing the output spreadsheets					
	with open('../Spreadsheets/node_sizes_by_clade.csv', 'w') as o:
		o.write('Gene Family,')
		for clade in size_by_clade_by_og[next(iter(size_by_clade_by_og))]:
			o.write(clade + ',')
		o.write('\n')
		
		for og in size_by_clade_by_og:
			o.write(og + ',')
			for clade in size_by_clade_by_og[og]:
				if(len(size_by_clade_by_og[og][clade]) == 0):
					o.write(' 0')
				for c, clade_size in enumerate(size_by_clade_by_og[og][clade]):
					if(len(size_by_clade_by_og[og][clade]) > 1 and c != len(size_by_clade_by_og[og][clade]) - 1):
						o.write(' ' + str(clade_size) + ';')
					elif(len(size_by_clade_by_og[og][clade]) == 1 or c == len(size_by_clade_by_og[og][clade]) - 1):
						o.write(' ' + str(clade_size))
				o.write(',')
			o.write('\n')
		
	if(combine and select_taxa):	
		with open('../Spreadsheets/removed_from_preguidance.csv', 'w') as o:
			for og in leaves_to_remove_by_og:
				o.write(og + ',' + ','.join(list(dict.fromkeys(leaves_to_remove_by_og[og]))) + '\n')
				
		select_taxa_from_preguidance(leaves_to_keep_by_og)
					
	return leaf_names_by_clade_by_og
	
	
#This function requires three spreadsheets from the above function, one for NBU clades, one for BB clades, and one for both/combined.
#It merges these measures into three rows per tree in a large easy-to-read spreadsheet, which are compared in later functions
def merge_clade_size_reports(nbu, bb, combined):

	tree_types = [nbu, bb, combined]
	
	nbu_minor_clade_counts, bb_minor_clade_counts, combined_minor_clade_counts = get_clade_totals()
	
	nbu_minor_clade_counts.update({ 'pl_gx' : nbu_minor_clade_counts['pl_gr'] + nbu_minor_clade_counts['pl_gl'] })
	nbu_minor_clade_counts.pop('pl_gr'); nbu_minor_clade_counts.pop('pl_gl')
		
	og_groups = { }
	for line in open('og_group_dataframe.csv'):
		og_groups.update({ line.split(',')[0].strip() : line.split(',')[1].strip() })

	clade_size_minimum = 3
	bacbin_clade_minor_clade_concentration_minimum = .6
	
	nbu_data = { }
	bb_data = { }
	combined_data = { }
	
	for t, type in enumerate(tree_types):
		type_data = { }
		type_clades = [] # A list so that we can index in the right order
		for l, line in enumerate(open(type)):
			if(l == 0):
				type_clades = [cell.strip() for cell in line.split(',')[1:] if(cell.strip() != '')]
			else:
				type_data.update({ line.split(',')[0].strip() : { } })
				for c, cell in enumerate([cell.strip() for cell in line.split(',')[1:] if(cell.strip() != '')]):
					type_data[line.split(',')[0].strip()].update({ type_clades[c] : [str(num) for num in cell.strip().split(';') if(int(num) >= clade_size_minimum)] })
					if(len(type_data[line.split(',')[0].strip()][type_clades[c]]) == 0):
						type_data[line.split(',')[0].strip()][type_clades[c]] = [' ']
						
		if t == 0: nbu_data = type_data
		elif t == 1: bb_data = type_data
		elif t == 2: combined_data = type_data
		
	for og in nbu_data:
		nbu_data[og].update({ 'pl_gx' : nbu_data[og]['pl_gr'] + nbu_data[og]['pl_gl'] })
		nbu_data[og].pop('pl_gr'); nbu_data[og].pop('pl_gl')
		
	bad_clades = ['am_th', 'ee_ce', 'op_is', 'sr_ch', 'am_hi', 'ee_ka', 'ex_is', 'sr_is', 'am_hb', 'ee_kb', 'ex_ib', 'sr_ib', 'am_va', 'am_vb', 'op_nu', 'op_nb', 'am_h', 'ee_k', 'ex_i', 'sr_i', 'am_v', 'op_n', 'am_uk', 'am_ub', 'am_u', 'ee_uk', 'ee_ub', 'ee_u', 'ex_ox', 'ex_ob', 'ex_o']
		
	nbu_keys = [clade for clade in sorted(nbu_data[next(iter(nbu_data))].keys()) if(clade[:2] != 'ba' and clade[:2] != 'za' and clade[:5] not in bad_clades)]
	bb_keys = [clade for clade in sorted(bb_data[next(iter(bb_data))].keys()) if(clade[:2] != 'ba' and clade[:2] != 'za' and clade[:5] not in bad_clades)]
	combined_keys = [clade for clade in sorted(combined_data[next(iter(combined_data))].keys()) if(clade[:2] != 'ba' and clade[:2] != 'za' and clade[:4] not in bad_clades)]
				
	master_data = { }
	for og in nbu_data:
		master_data.update({ og : { 'nbu_data' : { clade : nbu_data[og][clade] for clade in nbu_data[og] } } })
		master_data[og].update({ 'bb_data' : { clade : bb_data[og][clade] for clade in bb_data[og] } })
		master_data[og].update({ 'combined_data' : { clade : combined_data[og][clade] for clade in combined_data[og] } })

	with open('../Spreadsheets/master_tree_type_counts.csv', 'w') as o:
		o.write('Gene Family,Putative Recipient.,')
		for clade in combined_keys:
			if(len(clade) == 7):
				for nbu_cld in nbu_minor_clade_counts:
					if(len(nbu_cld) == 8):
						if(nbu_cld[:4] + nbu_cld[5:] == clade):
							o.write(clade + ' (' + str(nbu_minor_clade_counts[nbu_cld]) + '; ' + str(bb_minor_clade_counts[clade[:4] + 'b' + clade[4:]]) + '),')
							break
			else:
				for nbu_cld in nbu_minor_clade_counts:
					if(nbu_cld[:4] == clade and len(nbu_cld) == 5 and nbu_cld not in bad_clades):
						if(nbu_cld[:2] == 'am'):
							print(nbu_cld)
						o.write(clade + ' (' + str(nbu_minor_clade_counts[nbu_cld]) + '; ' + str(bb_minor_clade_counts[clade + 'b']) + '),')
						break
		o.write('\n')
		for og in master_data:
			o.write(og + ',' + og_groups[og] + ',')
			#o.write('\n,,')
			
			for clade in nbu_keys:
				o.write('; '.join(master_data[og]['nbu_data'][clade]) + ',')
			o.write('\n,,')
			
			for clade in bb_keys:
				o.write('; '.join(master_data[og]['bb_data'][clade]) + ',')
			o.write('\n,,')
			
			for clade in combined_keys:
				o.write('; '.join(master_data[og]['combined_data'][clade]) + ',')
			o.write('\n\n')	
			
	
#This function parses the large easy-to-read spreadsheet and pulls out trees that changed significantly after the addition of bacterial bins		
def pull_out_trees_of_interest():

	lines = [line for line in open('../Spreadsheets/master_tree_type_counts.csv')]
	clade_names = [cell.split('(')[0].strip() for cell in lines[0].split(',')[2:]]
	clades_by_og = { line.split(',')[0].strip() + '_' + line.split(',')[1].strip() : [[cell.strip() for cell in l.split(',')[2:]] for l in lines[i:i+3]] for i, line in enumerate(lines) if 'OG5_' in line }
	
	clades_dict = { }
	
	ogs_to_check = []
	for id in clades_by_og:
		check_flag = False
		
		og = id[:10]; recipient = id[11:]
		nbu = []; bb = []; comb = []
		for t, type in enumerate(clades_by_og[id]):
			nummed_list = []
			for string in clades_by_og[id][t]:
				if(string == '' or string == ' ' or string == ';' or string == '; ' or string == ' ;' or string == ' ; '):
					string = '0'
				nummed_list.append([int(val.strip()) for val in string.split(';') if val.strip() != ''])
			if(nbu == []): nbu = nummed_list
			elif(bb == []): bb = nummed_list
			elif(comb == []): comb = nummed_list
			
		clades_dict.update({ id : { 'nbu' : nbu, 'bb' : bb, 'comb' : comb } })
			
		for c, minor_clade in enumerate(comb):		
			if(recipient == 'EGT'):
				recip_clades = ['pl_g', 'pl_r', 'sr_s', 'sr_r', 'sr_d', 'ee_c', 'ee_h']
			elif(recipient == 'Am_Ex'):
				recip_clades = ['am_a', 'am_d', 'am_i', 'am_m', 'am_t', 'am_v', 'am_u', 'ex_e', 'ex_f', 'ex_h', 'ex_j', 'ex_m', 'ex_o', 'ex_p']
			else:
				recip_clades = [clade for clade in clade_names if clade[:2] == recipient.lower()]

			if(minor_clade != nbu[c] and clade_names[c] not in recip_clades and sum(minor_clade) - sum(nbu[c]) > 1):
				#print(og, clade_names[c], recipient, minor_clade, recip_clades)
				check_flag = True
				break
				
		if(check_flag):
			ogs_to_check.append(og)				
	
	if(not os.path.isdir('../Trees_to_check')):
		os.mkdir('../Trees_to_check')
		
	for og in ogs_to_check:
		for tree in os.listdir('../Subtrees/Fulltrees'):
			if(og in tree):
				os.system('cp ../Subtrees/Fulltrees/' + tree + ' ../Trees_to_check/' + tree)
				
	return clades_dict, clade_names
																															

#This function decides which taxa to include in each of the input trees based on the large spreadsheet 
#generated by above function and the sizes of both NBU, BB, and both/combined clades
def select_taxa_from_preguidance(leaves_per_og):

	if(not os.path.isdir('../Subtrees/Preguidance_Filtered')):
		os.mkdir('../Subtrees/Preguidance_Filtered')

	for og in leaves_per_og:
		pg_seqs = []
		for pg_file in os.listdir('../Subtrees/Preguidance'):
			if(og in pg_file):
				pg_seqs = list(SeqIO.parse('../Subtrees/Preguidance/' + pg_file, 'fasta'))
				break
				
		seqs_to_write = []
		if(len(pg_seqs) > 0):	
			for seq in pg_seqs:
				if(seq.id.replace('-', '') in leaves_per_og[og] or seq.id in leaves_per_og[og]):
					seqs_to_write.append(seq)
						
			with open('../Subtrees/Preguidance_Filtered/' + og + '_preguidance_filtered.fasta', 'w') as o:
				for seq in seqs_to_write:
					o.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
		else:
			print('\nNo pre-Guidance file could be found for OG ' + og + '\n')


#This function gets information on the abundance/concentration of each minor clade in input set of trees
#and is mainly used to generate the main heat-map figure (mother scorpion) of taxon/clade presence per OG
def minor_clade_counts():

	# major_clade_members = { }
# 	minor_clade_members = { }
# 	
# 	for file in os.listdir('../../../../ReadyToGo_Files'):
# 		if('ds_store' not in file.lower() and ((file[:2] == 'Ba' or file[:2] == 'Za') or (file[:2] != 'Ba' and file[:2] != 'Za' and file[4] != 'b'))):
# 			if(file[:2] not in major_clade_members):
# 				major_clade_members.update({ file[:2] : [] })
# 			if(file[:5] not in minor_clade_members):
# 				minor_clade_members.update({ file[:5] : [] })
# 				
# 			if(file[:10] not in major_clade_members[file[:2]]):
# 				major_clade_members[file[:2]].append(file[:10])
# 			if(file[:10] not in minor_clade_members[file[:5]]):
# 				minor_clade_members[file[:5]].append(file[:10])
# 			
# 	for f, file in enumerate(os.listdir('../../../../allOG5Files')):
# 		print(str(f) + '. ' + file)
# 		if('ds_store' not in file.lower()):
# 			taxa = [record.id[:10] for record in list(SeqIO.parse('../../../../allOG5Files/' + file, 'fasta'))]
# 			
# 			for taxon in taxa:
# 				if(taxon[:2] not in major_clade_members):
# 					major_clade_members.update({ taxon[:2] : [] })
# 				if(taxon[:5] not in minor_clade_members):
# 					minor_clade_members.update({ taxon[:5] : [] })
# 				
# 				if(taxon[:10] not in major_clade_members[taxon[:2]]):
# 					major_clade_members[taxon[:2]].append(taxon[:10])
# 				if(taxon[:10] not in minor_clade_members[taxon[:5]]):
# 					minor_clade_members[taxon[:5]].append(taxon[:10])
# 					
# 	print('')
# 				
# 	major_clade_totals = { clade : len(major_clade_members[clade]) for clade in major_clade_members }
# 	minor_clade_totals = { clade : len(minor_clade_members[clade]) for clade in minor_clade_members }
	
	major_clade_totals = {'Ba': 710, 'EE': 53, 'Sr': 464, 'Am': 244, 'Za': 115, 'Ex': 46, 'Pl': 81, 'Op': 124}
	minor_clade_totals = {'Za_ey': 1, 'Za_et': 3, 'Za_eu': 2, 'Za_ep': 4, 'Am_ar': 10, 'Za_em': 21, 'Za_eh': 28, 'Za_ea': 3, 'Za_eb': 6, 'Za_ec': 6, 'Sr_st': 98, 'Ba_el': 5, 'Am_va': 2, 'Ba_me': 1, 'Op_ch': 3, 'Ba_di': 2, 'Ex_ja': 7, 'Za_th': 4, 'Ba_df': 6, 'Ba_de': 10, 'Sr_rh': 91, 'EE_ce': 4, 'Ba_pl': 22, 'Za_cr': 25, 'Op_fu': 51, 'Ba_pd': 49, 'Am_th': 1, 'Ba_pb': 30, 'Ba_pa': 59, 'Am_tu': 165, 'Op_nu': 1, 'Pl_rh': 22, 'EE_cr': 15, 'Ba_cr': 2, 'Ba_cp': 9, 'Ex_is': 1, 'Ba_cv': 22, 'Ba_cy': 122, 'Ba_ca': 8, 'Ba_cd': 14, 'Ba_ch': 28, 'Ba_cn': 2, 'Ex_pa': 3, 'Am_my': 19, 'EE_uk': 1, 'Ba_sp': 23, 'Za_ba': 2, 'Ba_sy': 12, 'Sr_ch': 2, 'Sr_ci': 206, 'Za_ko': 1, 'Am_uk': 1, 'Ex_he': 6, 'Ba_bi': 2, 'Sr_pe': 3, 'Op_me': 65, 'Ba_nt': 2, 'EE_is': 6, 'EE_ap': 6, 'Ba_ni': 3, 'Ba_aq': 14, 'Ba_ar': 1, 'Ba_pg': 47, 'Ba_is': 1, 'Ba_ac': 39, 'Ba_ad': 20, 'Ex_ox': 1, 'EE_br': 4, 'Op_is': 1, 'Sr_ap': 24, 'Ba_bc': 20, 'Ex_fo': 6, 'Am_di': 33, 'Ex_eu': 19, 'Ba_le': 3, 'Ba_ba': 37, 'EE_ha': 17, 'Ba_ts': 11, 'Ba_ge': 3, 'Za_na': 1, 'Ex_ma': 3, 'Ba_th': 14, 'Am_is': 13, 'Sr_di': 40, 'Ba_te': 7, 'Ba_fn': 2, 'Pl_gr': 54, 'Ba_fe': 2, 'Za_al': 1, 'Ba_fc': 22, 'Ba_fb': 20, 'Za_as': 7, 'Op_ic': 3, 'Ba_fu': 12, 'Ba_ft': 2, 'Pl_gl': 5}
	
	ogs_grouped = {}
	for row in open('og_groups.csv'):
		row = row.split(',')
		if(row[1].strip() != ''):
			if(row[1].strip() not in ogs_grouped):
				ogs_grouped.update({ row[1].strip() : [] })
				
			ogs_grouped[row[1].strip()].append(row[0].strip())
							
	minor_clade_master = { 'Am_ar' : 0, 'Am_di' : 0, 'Am_my' : 0, 'Am_hi' : 0, 'Am_is' : 0, 'Am_tu' : 0, 'EE_ap' : 0, 'EE_br' : 0, 'EE_cr' : 0, 'EE_ha' : 0, 'EE_is' : 0, 'EE_ka' : 0, 'EE_ce' : 0, 'EE_uk' : 0, 'Ex_eu' : 0, 'Ex_fo' : 0, 'Ex_he' : 0, 'Ex_is' : 0, 'Ex_ja' : 0, 'Ex_ma' : 0, 'Ex_ox' : 0, 'Ex_pa' : 0, 'Op_ch' : 0, 'Op_fu' : 0, 'Op_ic' : 0,'Op_me' : 0, 'Pl_gl' : 0, 'Pl_gr' : 0, 'Pl_rh' : 0, 'Sr_ap' : 0, 'Sr_ci' : 0, 'Sr_di' : 0, 'Sr_pe' : 0, 'Sr_rh' : 0, 'Sr_st' : 0, 'Ba_ac' : 0, 'Ba_ad' : 0, 'Ba_aq' : 0, 'Ba_ar' : 0, 'Ba_ba' : 0, 'Ba_bc' : 0, 'Ba_bi' : 0, 'Ba_ca' : 0, 'Ba_cd' : 0, 'Ba_ch' : 0, 'Ba_cr' : 0, 'Ba_cv' : 0, 'Ba_cy' : 0, 'Ba_de' : 0, 'Ba_df' : 0, 'Ba_di' : 0, 'Ba_el' : 0, 'Ba_fb' : 0, 'Ba_fc' : 0, 'Ba_fn' : 0, 'Ba_fu' : 0, 'Ba_ge' : 0, 'Ba_is' : 0, 'Ba_me' : 0, 'Ba_ni' : 0, 'Ba_pa' : 0, 'Ba_pb' : 0, 'Ba_pd' : 0, 'Ba_pg' : 0, 'Ba_pl' : 0, 'Ba_sp' : 0, 'Ba_sy' : 0, 'Ba_te' : 0, 'Ba_th' : 0, 'Ba_ts' : 0, 'Ba_le' : 0, 'Ba_cp' : 0, 'Ba_fe' : 0, 'Ba_nt' : 0, 'Ba_cn' : 0, 'Ba_ft' : 0, 'Za_cr' : 0, 'Za_ea' : 0, 'Za_eb' : 0, 'Za_ec' : 0, 'Za_eh' : 0, 'Za_em' : 0, 'Za_ep' : 0, 'Za_et' : 0, 'Za_eu' : 0, 'Za_ey' : 0, 'Za_ko' : 0, 'Za_na' : 0, 'Za_nh' : 0, 'Za_pa' : 0, 'Za_th' : 0, 'Za_lo' : 0, 'Za_as' : 0, 'Za_ba' : 0, 'Za_ba' : 0, 'Za_al' : 0 }
	
	master_counts = minor_clade_master.copy()
	
	special_bacs = { 'Pl_gb' : 'Pl_gr' }
	
	major_clade_master = { 'Am' : 0, 'Sr' : 0, 'Op' : 0, 'EGT' : 0, 'Ex' : 0, 'EE' : 0, 'Am_Ex' : 0, 'Complex' : 0}
	
	mcs = list(minor_clade_master)
	
	counts_by_og = { }
	unique_counts_by_og = { }
	
	majors_by_og = { }
	
	euk_leaves_by_og = { }
	prok_leaves_by_og = { }
	
	singletons_to_remove = [line.strip() for line in open('singletons_to_remove.txt')]

	for f, file in enumerate(os.listdir('../Subtrees/Fulltrees')):
		if(file.endswith('.tre')):
			fname = '../Subtrees/Fulltrees/' + file
			og = 'OG5_' + fname.split('OG5_')[1][:6]
			
			if(og in plastid_ogs):
				continue
			
			gr = ''
			for group in ogs_grouped:
				if(og in ogs_grouped[group]):
					gr = group
			
			if(gr != 'EGT' and gr != 'Complex' and gr != 'Am_Ex'):
				recips_list = [gr]
			elif(gr == 'EGT'):
				recips_list = ['Sr', 'Pl', 'EE']
			elif(gr == 'Am_Ex'):
				recips_list = ['Am', 'Ex']
			else:
				continue
			
			print(str(f) + '. ' + og)
			
			seen_taxa = []
			
			counts_by_og.update({ og : minor_clade_master.copy() })
			unique_counts_by_og.update({ og : minor_clade_master.copy() })
			majors_by_og.update({ og : major_clade_master.copy() })
		
			var.trees = []			
			
			nexus_to_newick(fname)						
			correct_tree(fname)													
			tree_file = fname.split('.tre')[0] + '_temp.tre'
			read(tree_file) 
			try:
				tree = var.trees[0] 	
			except:
				continue					
												
			os.remove(tree_file)
							
			recip_clade = False
			if(not recip_clade):
				leaves = [leaf for leaf in tree.getAllLeafNames(0) if leaf not in singletons_to_remove and leaf[:10] not in singletons_to_remove and leaf[:5] not in singletons_to_remove]
				euk_leaves_by_og.update({ og : [leaf for leaf in leaves if(leaf[:2] != 'Ba' and leaf[:2] != 'Za')] })
				prok_leaves_by_og.update({ og : [leaf for leaf in leaves if(leaf[:2] == 'Ba' or leaf[:2] == 'Za')] })
			else:
				leaves = []
				
				for node in tree.iterNodesNoRoot(): 	
					if node.getNChildren() == 0:						
						taxon_of_interest = node.name[:10]
						full_leaf_name = node.name
																		
						num_sister_taxa = 0; other_rat = 0
						while num_sister_taxa < 3:
							sister_taxa = tree.getAllLeafNames(node.parent)
						
							num_sister_taxa = len(sister_taxa)
							
							n = 0
 							for sister_taxon in sister_taxa:
 								if(sister_taxon[:2] not in recips_list):
 									n += 1
 										
 								other_rat = float(n) / float(len(sister_taxa))
 						 
							node = node.parent
						
						if(other_rat < .25 or taxon_of_interest[:2] == 'Ba' or taxon_of_interest[:2] == 'Za'):
							#print(og, gr, taxon_of_interest, other_rat, recips_list, len(sister_taxa))						
							leaves.append(full_leaf_name)
									
				euk_leaves_by_og.update({ og : [leaf for leaf in leaves if(leaf[:2] != 'Ba' and leaf[:2] != 'Za')] })
				prok_leaves_by_og.update({ og : [leaf for leaf in tree.getAllLeafNames(0) if(leaf[:2] == 'Ba' or leaf[:2] == 'Za')] })
														
			for leaf in leaves:
				if((leaf[:2] != 'Ba' and leaf[:2] != 'Za' and leaf[4] != 'b') or (leaf[:2] == 'Ba' or leaf[:2] == 'Za')):
					key = leaf[:5]
				else:
					if(leaf[:5] in special_bacs):
						key = special_bacs[leaf[:5]]
					else:
						for clade in minor_clade_master:
							if(leaf[:4] in clade):
								key = clade
								break
								
				if(key in minor_clade_master):				
					counts_by_og[og][key] += 1
					if(leaf[:10] not in seen_taxa):
						unique_counts_by_og[og][key] += 1
						seen_taxa.append(leaf[:10])
					master_counts[key] += 1
					if(leaf[:2] == 'Pl'):
						majors_by_og[og]['EGT'] += 1
					elif(leaf[:2] != 'Ba' and leaf[:2] != 'Za'):
						majors_by_og[og][leaf[:2]] += 1
							
					if(leaf[:2] == 'Am' or leaf[:2] == 'Ex'):
						majors_by_og[og]['Am_Ex'] += 1
						
						
	# minor_clade_codes = [clade for clade in minor_clade_master]
# 	with open('../Spreadsheets/absolute_counts.csv', 'w') as o:
# 		o.write('Gene Family,')
# 		for clade in minor_clade_codes:
# 			o.write(clade + ',')
# 		o.write('\n')
# 		for og in counts_by_og:
# 			o.write(og + ',')
# 			for clade in minor_clade_codes:
# 				o.write(str(counts_by_og[og][clade]) + ',')
# 			o.write('\n')
# 			
# 	exit()
# # 														
# 	with open('../Spreadsheets/minor_clade_counts.csv', 'w') as o:	
# 		o.write('Gene Family,Group,N Recips,N Minor Clades,Maj Rat,Minor Clade,OutOfDom,OutOfMaj,OutOfMin\n')
# 		for og in counts_by_og:
# 			gr = ''
# 			for group in ogs_grouped:
# 				if(og in ogs_grouped[group]):
# 					gr = group
# 			if(gr != '' and sum(majors_by_og[og].values()) > 0):
				
# 				#print(og, gr, majors_by_og[og])
						
# 				maj_rat = float(majors_by_og[og][gr])/float(sum(majors_by_og[og].values()))
# 				#n_mins = len(dict.fromkeys([leaf[:5] for leaf in prok_leaves_by_og[og]]))
# 				n_mins = float(len([leaf for leaf in prok_leaves_by_og[og] if leaf[:5] == 'Ba_cy'])) / float(len(prok_leaves_by_og[og]))
# 				#n_mins = len([leaf for leaf in prok_leaves_by_og[og] if leaf[:5] == 'Ba_cy'])
# 				for mc in mcs:
# 					if(master_counts[mc] >= 10):
# 						if(mc[:2] != 'Ba' and mc[:2] != 'Za'):
# 							OutOfDom = float(counts_by_og[og][mc]) / float(len(euk_leaves_by_og[og]))
# 						else:
# 							OutOfDom = float(counts_by_og[og][mc]) / float(len(prok_leaves_by_og[og]))
# 						OutOfMaj = float(unique_counts_by_og[og][mc])/float(major_clade_totals[mc[:2]])
# 						OutOfMin =	float(unique_counts_by_og[og][mc])/float(minor_clade_totals[mc])
# 						o.write(og + ',' + gr + ',' + str(majors_by_og[og][gr]) + ',' + str(n_mins) + ',' + str(maj_rat) + ',' + mc + ',' + str(OutOfDom) + ',' + str(OutOfMaj) + ',' + str(OutOfMin) + '\n')	
						
# 	exit()
												
	# euk_clades = ['Ex', 'Am']
# 	all_rel_euks = list(dict.fromkeys([leaf[:10] for og in euk_leaves_by_og for leaf in euk_leaves_by_og[og] if leaf[:2].lower() in [clade.lower() for clade in euk_clades]]))
# 	all_ex = [euk for euk in all_rel_euks if euk[:2] == 'Ex' and euk[:10] in open('excavata.tre').readlines()[0]]
# 	all_am = [euk for euk in all_rel_euks if euk[:2] == 'Am' and euk[:10] in open('amoebozoa.tre').readlines()[0]]
# 	with open('../Spreadsheets/am_ex_matrix.csv', 'w') as o:
# 		o.write(',' + ','.join(all_ex) + '\n')
# 		cocurr = { am : { ex : 0 for ex in all_ex } for am in all_am }
# 		for og in euk_leaves_by_og:
# 			for leaf in euk_leaves_by_og[og]:
# 				if(leaf[:10] in cocurr and leaf[:10] in open('amoebozoa.tre').readlines()[0]):
# 					for ex in cocurr[leaf[:10]]:
# 						if(ex in [l[:10] for l in euk_leaves_by_og[og]] and ex[:10] in open('excavata.tre').readlines()[0]):
# 							cocurr[leaf[:10]][ex] += 1
# 							
# 		for am in all_am:
# 			o.write(am + ',' + ','.join([str(cocurr[am][ex]) for ex in all_ex]) + '\n')
						
	#euk_clades = ['Pl', 'Sr', 'EE', 'Ex']; recip = 'EGT'
	euk_clades = ['Op']; recip = 'Op'
	presence_by_og = { }
	all_rel_euks = list(dict.fromkeys([leaf[:10] for og in euk_leaves_by_og for leaf in euk_leaves_by_og[og] if leaf[:2].lower() in [clade.lower() for clade in euk_clades]]))
	
	#all_rel_proks = list(dict.fromkeys([leaf[:10] for og in prok_leaves_by_og for leaf in prok_leaves_by_og[og] if leaf[:10] in [line.split('\t')[0].strip() for line in open('proteo_cyano_lineages.tsv')]]))
		
	with open('../Spreadsheets/' + recip + '_taxon_presence_long.csv', 'w') as o:
		o.write('Gene.Family,label,val,Num.Fungi,Num.Met\n')
		#o.write('Gene.Family,label,val,Prop.Pl\n')
		for og in euk_leaves_by_og:
			flag = 0
			for euk_clade in euk_clades:
				if(og in ogs_grouped[recip]):
					flag = 1
					break
			if(flag == 1):
				num_fungi = len(list([leaf for leaf in euk_leaves_by_og[og] if leaf[:5].lower() == 'op_fu']))
				num_met = len(list([leaf for leaf in euk_leaves_by_og[og] if leaf[:5].lower() == 'op_me']))
				#prop_pl = float(len(list([leaf for leaf in euk_leaves_by_og[og] if leaf[:2].lower() == 'pl']))) / float(len(list([leaf for leaf in euk_leaves_by_og[og]])))
				for euk in all_rel_euks:
					if(euk in [leaf[:10] for leaf in euk_leaves_by_og[og]]):
						val = '1'
					else:
						val = '0'
					
					o.write(og + ',' + euk + ',' + val + ',' + str(num_fungi) + ',' + str(num_met) + '\n')
					#o.write(og + ',' + euk + ',' + val + ',' + str(prop_pl) + '\n')
				#for prok in all_rel_proks:
				#	if(prok in [leaf[:10] for leaf in prok_leaves_by_og[og]]):
				#		val = '1'
					#else:
					#	val = '0'
					#
					#o.write(og + ',' + prok + ',' + val + ',' + str(prop_pl) + '\n')

	exit()
					
	with open('../Spreadsheets/' + recip + '_taxon_presence_wide.csv', 'w') as o:
		o.write('Gene.Family,' + ','.join(all_rel_euks) + '\n')
		for og in euk_leaves_by_og:
			o.write(og + ',')
			flag = 0
			for euk_clade in euk_clades:
				if(og in ogs_grouped[recip]):
					flag = 1
					break
			if(flag == 1):
				#num_fungi = len(list([leaf for leaf in euk_leaves_by_og[og] if leaf[:5].lower() == 'op_fu']))
				for euk in all_rel_euks:
					if(euk in [leaf[:10] for leaf in euk_leaves_by_og[og]]):
						val = '1'
					else:
						val = '0'
					
					o.write(val + ',')
			o.write('\n')
				
	exit()				
			

############################################################
			
			   	#Main wrapper
			
############################################################
			
			
#Calling all of the above functions
def main():

	minor = True #Report by minor clades instead of major clades (columns)
	bacbin = False #Report only bacterial bin minor clades (columns)
	combine = True #Report bacterial bin and euk minor clade counts SUMMED (values)
	select_taxa = False
	
	#What is not an option (above) yet is counting only bacterial bins at the major clade level

	if(not os.path.isdir('../Spreadsheets')): os.system('mkdir ../Spreadsheets')

	#prokaryotes_only, include_prokaryotes = get_args()
	
	#leaves_by_major_clade = get_all_clade_sizes(prokaryotes_only, include_prokaryotes, minor, bacbin, combine, select_taxa)
	
	nbu = '../Spreadsheets/node_sizes_by_clade_NBU.csv'
	bb = '../Spreadsheets/node_sizes_by_clade_BB.csv'
	combined = '../Spreadsheets/node_sizes_by_clade_combined.csv'
	
	#merge_clade_size_reports(nbu, bb, combined)
	
	#clades_dict, clade_names = pull_out_trees_of_interest()
	
	#choose_taxa_to_use_per_tree(clades_dict, clade_names)
	
	#minor_clade_counts()
	
	#diversity_in_euk_clade()
	
	diversity_in_full_tree()	
	
	
main()
