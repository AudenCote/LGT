#Dependencies
import dendropy
from p4 import *
import os, re
import sys
import csv
from Bio import SeqIO

var.doRepairDupedTaxonNames = 1

				
# Some trees have '-', which causes problems in p4
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 

	tree_corrected.write(tree2correct)
	tree_corrected.close()
	
#Utility function for later
def remove_bacbin_taxa(taxa):
	updated_taxa = []
	for taxon in taxa:
		if(taxon[:2] == 'Ba' or taxon[:2] == 'Za'):
			updated_taxa.append(taxon)
		elif(taxon[:2] != 'Ba' and taxon[:2] != 'Za' and taxon[4] != 'b'):
			updated_taxa.append(taxon)

	return updated_taxa
	
	
def get_singletons(fname):

	singleton_orfs = []
	
	t = fname.split('/')[-1]	
											
	if('OG5' in t):					
		OG5 = t.split('.')[1].replace('_postguidance', '')
		
		var.trees = []								
						
		correct_tree(fname)								
													
		tree_file = fname.split('.tre')[0] + '_temp.tre'	
		read(tree_file) 							
		tree = var.trees[0] 						
													
		os.remove(fname.split('.tre')[0] + '_temp.tre')
		
		for node in tree.iterNodesNoRoot(): 	
			if node.getNChildren() == 0: 							
				taxon_of_interest = node.name[:10]
				full_leaf_name = node.name
				if(taxon_of_interest[:2] != 'Ba' and taxon_of_interest[:2] != 'Za'):
					append_taxon = True
										
					num_sister_taxa = 0
					while num_sister_taxa < 4:
						sisterTaxa = tree.getAllLeafNames(node.parent)	
					
						#Removing bacterial bins from sisterTaxa
						#updated_sister_taxa = remove_bacbin_taxa(sisterTaxa)
						updated_sister_taxa = sisterTaxa
						num_sister_taxa = len(updated_sister_taxa)
						
						node = node.parent
					
					num_other_euks = 0
					n = 0
 					for sister_taxon in updated_sister_taxa:
 						if(taxon_of_interest[:4] + '_' + taxon_of_interest[6:8] == 'Am_t_he' or taxon_of_interest[:4] + '_' + taxon_of_interest[6:8] == 'Am_t_he'):
 							if(sister_taxon[:4] + '_'  + sister_taxon[6:8] == taxon_of_interest[:4] + '_' + taxon_of_interest[6:8] and sister_taxon[:10] != taxon_of_interest):
 								n += 1
 						else:
 							if(sister_taxon[:4] == taxon_of_interest[:4] and sister_taxon[:10] != taxon_of_interest):
 								n += 1
 				
 					#if(float(n) / float(len(updated_sister_taxa)) >= .25):
 					if(n > 0):
 						append_taxon = False
 						
 					if(append_taxon): 
 						singleton_orfs.append(full_leaf_name)
	
	return list(dict.fromkeys(singleton_orfs))
	