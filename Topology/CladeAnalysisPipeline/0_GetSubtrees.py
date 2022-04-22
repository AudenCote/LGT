#Dependencies
from p4 import *
import os, re
import sys
import csv
from Bio import SeqIO
import subprocess

var.doRepairDupedTaxonNames = 1

euk_major_clades = [ 'am', 'sr', 'op', 'pl', 'ex', 'ee' ]

def bad_script_call():

	print('\nPlease input a .tre file or a directory or .tre files and a eukaryotic major clade of interest\n\npython find_subtrees.py --input_dir <folder of .tre files> --euk_clade Am\n')

	exit()


def get_args():

	trees_handle = ''
	directory = False
	euk_clade = ''
	
	try:
		if(sys.argv[1] == '--input_file'):
			trees_handle = sys.argv[2]
		elif(sys.argv[1] == '--input_dir'):
			trees_handle = sys.argv[2]
			directory = True
		else:
			bad_script_call()
	except IndexError:
		bad_script_call()

	try:
		if(sys.argv[3] == '--euk_clade' and sys.argv[4].lower() in euk_major_clades):
			euk_clade = sys.argv[4].lower()
		else:
			euk_clade = 'file'
			#bad_script_call()
	except IndexError:
		euk_clade = 'file'
		#bad_script_call()

	
	return [ trees_handle, directory, euk_clade ]
	
	
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
		
	if(tree2correct[0].startswith('(') or tree2correct.startswith('\t(')):
		newick = tree2correct
		newick = newick.split(',')
		
		to_join = []
		for leaf in newick:
			if(leaf.startswith('(') and ')' in leaf):
				leaf = leaf.replace('(', '').split(')')
				del leaf[1]
				leaf = ')'.join(leaf)
				to_join.append(leaf)
			else:
				to_join.append(leaf)
				
		corrected = ','.join(to_join)
		
	tree_corrected.write(corrected)
	tree_corrected.close()	
	
	
def get_largest_euk_clade(tree, euk_clades, nodes_to_exclude):

	print(euk_clades)

	best_node = None
	best_size = 0
	contam_threshold = .125
	#Setting this to zero because filtering occurs after the data has been collected; if you would like to filter for best node size raise this value (e.g. absolute_n_threshold = 3)
	absolute_n_threshold = 0
	
	forbidden_nodes = [node for node in nodes_to_exclude]
	for node in nodes_to_exclude:
		for num in tree.getNodeNumsAbove(node):
			forbidden_nodes.append(tree.node(num))
		
	for node in tree.iterNodesNoRoot():
		if(len(tree.getAllLeafNames(node)) > absolute_n_threshold and node not in forbidden_nodes):
			leaves = tree.getAllLeafNames(node)
							
			num = 0.0; dem = 0.0;
		
			for leaf in leaves:
				if(leaf[:2] in euk_clades):
					dem += 1.0
				else:
					dem += 1.0; num += 1.0;
					
			if(dem == 0.0):
				ratio = 0.0
			else:
				ratio = num / dem
	
			if(ratio <= contam_threshold and dem - num > best_size):
				best_node = node
				best_size = dem
						
	return best_node
	

def get_prok_clade_counts(tree, starting_node):

	num_proks_threshold = 5
	max_minor_clade = False
	num_proks = 0
	
	i = 0
	while num_proks <= num_proks_threshold and max_minor_clade == False:
		if(i > 0): starting_node = starting_node.parent
		
		proks = [leaf for leaf in tree.getAllLeafNames(starting_node) if(leaf[:2].lower() == 'ba' or leaf[:2].lower() == 'za')]
		num_proks = len(proks)
		
		if(num_proks > num_proks_threshold):
			minor_clade_counts = {}
			for prok in proks:
				if(prok in minor_clade_counts):
					minor_clade_counts[prok] += 1
				else:
					minor_clade_counts.update({ prok : 1 })
				
			counts = [ minor_clade_counts[key] for key in minor_clade_counts ]; filtered = list(dict.fromkeys(counts));
			if(len(filtered) == len(counts)):
				max_minor_clade = True
				
		i += 1
			
	return starting_node, minor_clade_counts, num_proks
	
	
def get_alignments(seq_names, fname):

	alignments = []
	
	og = fname.split('/')[-1].split('.')[1][:10]
			
	n = 0
	for file in os.listdir('./' + '/'.join(fname.split('/')[:-1])):
		if('postguidance.fas' in file and '95gapTrimmed' not in file and og in file and '.tre' not in file):
			n += 1
			if(n > 1):
				print('\nMore than one post-guidance alignment file found for gene family OG5_' + str(file.split('OG5_')[-1][:6]) + '\n')
				exit()
				
			records = [(record.id, record.seq) for record in list(SeqIO.parse('./' + '/'.join(fname.split('/')[:-1]) + '/' + file, 'fasta')) if record.id in seq_names]
			
			#Writing the updated post-Guidance alignments
			with open('../Subtrees/Postguidance/' + og + '_subtree_postguidance_Original.fasta', 'w') as o:
				for record in records:
					o.write('>' + record[0] + '\n')
					o.write(str(record[1]) + '\n')
				
			for record in records:
				alignments.append(record[0] + '\t' + record[1])
		
		#Storing the relevant pre-Guidance information
		if('preguidance' in file and og in file):
			records = [(record.id, record.seq) for record in list(SeqIO.parse('./' + '/'.join(fname.split('/')[:-1]) + '/' + file, 'fasta')) if record.id in seq_names]
			
			with open('../Subtrees/Preguidance/' + og + '_subtree_preguidance_Original.fasta', 'w') as o:
			
				for record in records:
					o.write('>' + record[0] + '\n')
					o.write(str(record[1]) + '\n')
			
	if(n == 0):
		print('\nNo guidance files found for gene family OG5_' + file.split('OG5_')[-1][:6] + '\n')
		exit()
				
	return alignments
		
				
def write_nexus(subtree, alignments, path):

	taxa = subtree.getAllLeafNames(0)
	
	newick = subtree.writeNewick(toString = True).replace('subTreeRoot', '')
	
	nchar = len(alignments[0].split('\t')[-1])
	
	with open(path, 'w') as o:
		ntax = str(len(taxa))
			
		o.write('#NEXUS\n')	
		o.write('begin taxa;\n')
		o.write('\tdimensions ntax=' + ntax + ';\n')
		o.write('\ttaxlabels\n')
				
		for taxon in taxa:
			o.write('\t' + taxon + '\n')
				
		o.write(';\nend;\n\n')
		
		o.write('begin data;\n')
		o.write('\tdimensions ntax=' + str(ntax) + ' nchar=' + str(nchar) + ';\n')
		o.write('\tformat datatype=protein gap=-;\n')
		o.write('\tmatrix\n')
		
		for alignment in alignments:
			o.write('\t\t' + str(alignment) + '\n')
		
		o.write('\t;\nend;\n\n')
		
		
		o.write('begin trees;\n')
		o.write('\ttree tree_1 = [&R]\n')
		o.write(newick)
		o.write('end;\n\n')
				
		with open('figtree_format.txt', 'r') as ff:
			for line in ff:
				o.write(line)
	
		
def get_subtree(fname, euk_clade, directory):

	if(euk_clade != 'EGT' and euk_clade != 'Complex' and euk_clade != 'Am_Ex'):
		euk_clades = [euk_clade]
	elif(euk_clade == 'EGT'):
		euk_clades = ['Sr', 'Pl', 'EE']
	elif(euk_clade == 'Am_Ex'):
		euk_clades = ['Am', 'Ex']
	else:
		return None
		
	print(fname)
		
	var.trees = []
	
	nexus_to_newick(fname)
														
	correct_tree(fname)								
												
	tree_file = fname.split('.tre')[0] + '_temp.tre'
	read(tree_file) 
	tree = var.trees[0] 						
												
	os.remove(tree_file)
	
	clade_node = get_largest_euk_clade(tree, euk_clades, [])
	
	print(clade_node)
	
	try:
		euk_clade_tree = tree.dupeSubTree(clade_node, up = True)
	except:
		return None
	
	if(directory):
		euk_clade_tree.writeNewick('../Subtrees/EukClades/' + fname.split('/')[-1].split('.tre')[0] + '_euk_clade.tre')
	else:
		euk_clade_tree.writeNewick('./' + fname.split('.tre')[0] + '_euk_clade.tre')
	
	#Getting the ba/za minor clade in most abundance closest to the euk clade
	prok_concentration_threshold = .125
	prok_concentration = 0
	
	i = 0
	while i <= 3 and prok_concentration <= prok_concentration_threshold:
	
		starting_node, prok_minor_clade_counts, num_proks = get_prok_clade_counts(tree, clade_node)
				
		num_irrelevant_taxa = len([ leaf for leaf in tree.getAllLeafNames(starting_node) if leaf[:2].lower() not in euk_clades ])
		prok_concentration = float(num_proks) / float(num_proks + num_irrelevant_taxa)
		
		i += 1
						
	if(prok_concentration <= prok_concentration_threshold):
		print('\nA suitable concentration of Bacteria and Archaea could not be found near any ' + euk_clade.upper() + ' clade\n')
		return None
		
	#Going back nodes until the majority minor clade disappears
	prok_minor_clade = ''
	inital_max = 0
	for mc in prok_minor_clade_counts:
		if(prok_minor_clade_counts[mc] > inital_max):
			initial_max = prok_minor_clade_counts[mc]
			prok_minor_clade = mc
	
	original = starting_node
	starting_node = starting_nodet
	runout = False	
	if(starting_node == None):
		starting_node = original
		runout = True
	rel_clade_count = 0
	num_proks = 0
	
	print('1')
	
	for leaf in tree.getAllLeafNames(starting_node):
		if(leaf[:5].lower() == prok_minor_clade):
			rel_clade_count += 1
		if(leaf[:2] == 'Ba' or leaf[:2] == 'Za'):
			num_proks += 1
	
	while rel_clade_count > initial_max and not runout: #or num_proks < .3 * len([leaf for leaf in tree.getAllLeafNames(starting_node)]) or num_proks < 10 and not runout:
		
		initial_max = rel_clade_count
	
		starting_node = starting_node.parent
		if(starting_node == None):
			break
	
		rel_clade_count = 0
		for leaf in tree.getAllLeafNames(starting_node):
			if(leaf[:5].lower() == prok_minor_clade):
				rel_clade_count += 1
			if(leaf[:2] == 'Ba' or leaf[:2] == 'Za'):
				num_proks += 1
				
	

	try:
		subtree = tree.dupeSubTree(starting_node, up = True)
	except:
		return None
		
	
		
	subtree.writeNewick('../Subtrees/Subtrees/' + fname.split('/')[-1].split('.tre')[0] + '_subtree.tre')
	
	return None
				
	alignments = get_alignments([leaf for leaf in subtree.getAllLeafNames(0)], fname)
	
	if(directory):
		write_nexus(subtree, alignments, '../Subtrees/Subtrees/' + fname.split('/')[-1].split('.tre')[0] + '_best_subtree_Original.tre')
	else:
		write_nexus(subtree, alignments, './' + fname.split('.tre')[0] + '_best_subtree.tre')
		
		
def main():
	
	trees_handle, directory, euk_clade = get_args()
			
	if(not directory):
		if(euk_clade == 'file'):
			from_file = [line.split(',')[1].strip() for line in open('og_groups.csv') if file.split('OG5')[-1][:6] in line][0]
									
			get_subtree(fname, from_file, directory)
		else:
			get_subtree(fname, euk_clade, directory)
	else:
		if(not os.path.isdir('../Subtrees')):
			os.system('mkdir ../Subtrees')
		if(not os.path.isdir('../Subtrees/Subtrees')):
			os.system('mkdir ../Subtrees/Subtrees')
		if(not os.path.isdir('../Subtrees/EukClades')):
			os.system('mkdir ../Subtrees/EukClades')
		if(not os.path.isdir('../Subtrees/Fulltrees')):
			os.system('mkdir ../Subtrees/Fulltrees')
		if(not os.path.isdir('../Subtrees/Preguidance')):
			os.system('mkdir ../Subtrees/Preguidance')
		#if(not os.path.isdir('../Subtrees/Preguidance_subtrees')):
		#	os.system('mkdir ../Subtrees/Preguidance_subtrees')
		if(not os.path.isdir('../Subtrees/Postguidance')):
			os.system('mkdir ../Subtrees/Postguidance')
	
		for f, file in enumerate(os.listdir(trees_handle)):
			if('OG5_' in file and '.tre' in file):
				fname = trees_handle + '/' + file
				os.system('cp ' + fname + ' ../Subtrees/Fulltrees/' + file.split('.tre')[0] + '_fulltree_Original.tre')
								
				if(euk_clade == 'file'):
					try:
						from_file = [line.split(',')[1].strip() for line in open('og_groups.csv') if file.split('OG5')[-1][:6] in line][0]
					except IndexError:
						continue
					get_subtree(fname, from_file, directory)
				else:
					get_subtree(fname, euk_clade, directory)
													
	
main()
	