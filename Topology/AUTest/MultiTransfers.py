import os
import re
import sys
#from p4 import *

#var.doRepairDupedTaxonNames = 1

ogs_to_test = [line.strip() for line in open('OGList.txt')]


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
			if(line.startswith('(') or line.startswith('tree1=')):
				newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')
				
		with open(fname, 'w') as o:
			o.write(newick)
			
			
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
	
	
def get_largest_euk_clade(tree, nodes_to_exclude, recips):

	best_node = None
	best_size = 0
	contam_threshold = .125
	absolute_n_threshold = 0
	
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
				
				if(leaf[:2] in recips or leaf[:4] in recips):
					dem += 1.0
				else:
					dem += 1.0; num += 1.0;
					
			if(dem == 0.0):
				ratio = 0.0
			else:
				ratio = num / dem
			
			if((ratio <= contam_threshold or (dem >= 3 and num <= 1)) and dem > best_size and dem - num > 1):
				best_node = node
				best_size = dem
						
	return best_node
	
	
def count_nodes(in_dir):

	if(not os.path.isdir('nmTrees')):
		os.mkdir('nmTrees')


	og_groups = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('og_groups.csv') }
	
	num_clades_by_og = { }
	for f, file in enumerate(os.listdir(in_dir)):
		if(file.endswith('.tre')):
			print(str(f) + '. ' + file)
			fname = in_dir + '/' + file
			og = 'OG5_' + fname.split('OG5_')[1][:6]
			try:
				gr = og_groups[og]
			except:
				continue
			
			try:
				if(gr == 'EGT'):
					recips = ['Sr', 'EE', 'Pl', 'Ex_e']
				elif(gr == 'Am_Ex'):
					recips = ['Sr', 'EE_b', 'Am_a', 'Ex']
				elif(gr == 'Complex'):
					continue
				else:
					recips = [gr]
			except:
				continue
		
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
			
			actual_nodes = []

			best_node = get_largest_euk_clade(tree, actual_nodes, recips)
			if(best_node != None):
				actual_nodes.append(best_node)
							
			while best_node != None:
				best_node = get_largest_euk_clade(tree, actual_nodes, recips)
				actual_nodes.append(best_node)
					
			num_clades_by_og.update({ og : len([node for node in actual_nodes if node != None]) })
			
	with open('euk_clade_counts_by_og.csv', 'w') as o:
		o.write('OG,N.Clades\n')
		for og in num_clades_by_og:
			o.write(og + ',' + str(num_clades_by_og[og]) + '\n')
			
			tre_f = [file for file in os.listdir(in_dir) if og in file and file.endswith('.tre')][0]
			if(num_clades_by_og[og] > 1):
				os.system('cp ' + in_dir + '/' + tre_f + ' nmTrees/' + tre_f)
				
				
def constrain_trees(in_dir):

	if(not os.path.isdir('ConstraintTrees')):
		os.mkdir('ConstraintTrees')

	for file in os.listdir(in_dir):
		if('OG5_' + file.split('OG5_')[-1][:6] in ogs_to_test):
			print(file)
		if(file.endswith('.tre') and 'OG5_' + file.split('OG5_')[-1][:6] in ogs_to_test):
			fname = in_dir + '/' + file
			og = 'OG5_' + fname.split('OG5_')[1][:6]
		
			for line in open(fname):
				if(line.startswith('(') or line.startswith('tree1=')):
					original_newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')
			
			var.trees = []			
			
			nexus_to_newick(fname)								
			correct_tree(fname)													
			tree_file = fname.split('.tre')[0] + '_temp.tre'
			read(tree_file) 
			#tre:
			tree = var.trees[0] 	
			#except:
				#continue					
												
			os.remove(tree_file)
			
			leaves = tree.getAllLeafNames(0)
			
			with open('ConstraintTrees/' + og + '_constraints.treels', 'w') as o:
				o.write(original_newick + '\n\n' + '((' + ','.join([leaf.replace('LKH', '-LKH').replace('--LKH', '-LKH') for leaf in leaves if leaf[:2] != 'Ba' and leaf[:2] != 'Za']) + '),' + ','.join([leaf for leaf in leaves if leaf[:2] == 'Ba' or leaf[:2] == 'Za']) + ');')
				
				
def AU_test(in_dir):

	if(not os.path.isdir('AUTesting')):
		os.mkdir('AUTesting')

	for file in os.listdir('ConstraintTrees'):
		if(file.endswith('.treels')):
			og = file.split('_constraints.treels')[0]
			
			if(not os.path.isdir('AUTesting/' + og)):
				os.mkdir('AUTesting/' + og)
				
			for align in os.listdir(in_dir):
				if(og in align and 'gapTrim' in align and '.tre' not in align):
					os.system('cp ' + in_dir + '/' + align + ' AUTesting/' + og + '/' + align)
					break
					
			for treels in os.listdir('ConstraintTrees'):
				if(og in treels and 'constraint' in treels):
					os.system('cp ConstraintTrees/' + treels + ' AUTesting/' + og + '/' + treels)
					break
					
			os.system('./AUTest_IQTree.sh ' + og + ' ' + align + ' ' + treels)


def summarize_tests():

	with open('AUTestPVals.csv', 'w') as o:
		o.write('OG,p-Value\n')
		for og_dir in os.listdir('AUTesting'):
			if(os.path.isdir('AUTesting/' + og_dir)):
				pval = None
				for file in os.listdir('AUTesting/' + og_dir):
					if(file.endswith('.iqtree')):
						lines = [line for line in open('AUTesting/' + og_dir + '/' + file)]
						for l, line in enumerate(open('AUTesting/' + og_dir + '/' + file)):
							if('p-AU' in line):
								pval = lines[l + 3].split(' ')[-3].strip()
								break

				if(pval != None):
					o.write(og_dir + ',' + pval + '\n')
				else:
					o.write(og_dir + ',FAILED\n')


			
			
def main():

	#in_dir = sys.argv[1]

	#count_nodes(in_dir)
	
	#constrain_trees(in_dir)
	
	#AU_test(in_dir)

	summarize_tests()
	
	
main()

