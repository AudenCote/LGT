import os
import sys
import re
from p4 import *
from Bio import SeqIO
import subprocess
import matplotlib


var.doRepairDupedTaxonNames = 1


def bad_script_call():
	
	print('\nPlease input a directory of trees and pre-Guidance files, e.g.\n\n\tpython FunctionalCategorization.py --input_dir ../Trees\n')
	exit()


def get_args():
	
	input_dir = ''
	
	if('--input_dir' in sys.argv):
		try:
			input_dir = sys.argv[sys.argv.index('--input_dir') + 1]
		except IndexError:
			bad_script_call()
	else:
		bad_script_call()
		
	if(input_dir.endswith('/')):
		return input_dir[:-1]
	else:
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
	
	
def revigo_reduce_list(go_terms):
	
	go_terms = '\n'.join(go_terms)
	os.system('/usr/local/bin/Rscript --vanilla query-revigo.r ' + go_terms)
	

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
				
			os.system('mkdir ../FunctionSearching/' + og)
			os.system('cp ' + input_dir + '/' + tree_file + ' ../FunctionSearching/' + og + '/' + og + '.tre')
			
			fname = '../FunctionSearching/' + og + '/' + og + '.tre'
			var.trees = []
			nexus_to_newick('../FunctionSearching/' + og + '/' + og + '.tre')											
			correct_tree('../FunctionSearching/' + og + '/' + og + '.tre')														
			tree_file = '../FunctionSearching/' + og + '/' + og + '_temp.tre'
			read(tree_file) 
			tree = var.trees[0] 												
			os.remove(tree_file)
			
			leaves = list(tree.getAllLeafNames(0))
			unfiltered_seqs = list(SeqIO.parse(preguidance_file, 'fasta'))
			seqs_in_tree = []
			for record in unfiltered_seqs:
				if(record.id in leaves):
					seqs_in_tree.append(record)
			
			with open('../FunctionSearching/' + og + '/' + og + '_all_seqs.fasta', 'w') as o:
				for record in seqs_in_tree:
					o.write('>' + record.id + '\n' + str(record.seq) + '\n\n')

	
def emap_selected_seqs():
	
	for dir in os.listdir('../FunctionSearching'):
		if('OG5_' in dir):
			os.system('emapper.py -i ../FunctionSearching/' + dir + '/' + dir + '_all_seqs.fasta  --output eggnog_results.fa --output_dir ../FunctionSearching/' + dir + ' -m diamond')
			
			os.system('mv ../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.annotations ../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.annotations.tsv')
			os.system('rm ../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.seed_orthologs')
				
				
def master_GO_presence():
	
	func_by_seq = { }
	
	for dir in os.listdir('../FunctionSearching'):
		if('OG5_' in dir):
			for line in open('../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.annotations.tsv'):
				if(line[0] != '#'):
					line = line.split('\t')
					go_str = ';'.join(line[12].split(',')).replace('"', '').strip()					
					if(go_str != '-'):
						func_by_seq.update({ line[0].strip() : { 'OG5' : dir, 'GO_list' : go_str.split(';') } })
						
	all_gos = []
	for seq in func_by_seq:
		for term in func_by_seq[seq]['GO_list']:
			if(term not in all_gos):
				all_gos.append(term)
							
	slim_tuples, slim_list = go_slim_list(all_gos)
	
	print('\n' + str(len(all_gos)) + ' Gene Ontology terms slimmed to ' + str(len(slim_list)) + ' terms\n')
		
	for seq in func_by_seq:
		func_by_seq[seq].update({ 'slimmed_list' : { term : 0 for term in slim_list } })
		used_slims = [t for tup in slim_tuples for t in tup[1] if tup[0] in func_by_seq[seq]['GO_list']]
		for term in func_by_seq[seq]['slimmed_list']:		
			
			print(used_slims)
			
			if term in used_slims:
				func_by_seq[seq]['slimmed_list'][term] = 1
				
	with open('../master_go_presence.csv', 'w') as o:
		o.write('Sequence,OG5,')
		for term in slim_list:
			o.write(term + ',')
		o.write('\n')
		for seq in func_by_seq:
			o.write(seq + ',' + func_by_seq[seq]['OG5'] + ',')
			for term in slim_list:
				o.write(str(func_by_seq[seq]['slimmed_list'][term]) + ',')
			o.write('\n')
			
			
def map_unique_vectors_to_trees():

	colors = ['[&!color=' + col + ']' for col in list(matplotlib.colors.cnames.values())]

	#colors = ['[&!color=#0000ff]', '[&!color=#ab2121]', '[&!color=#7b25aa]', '[&!color=#006300]', '[&!color=#ffa100]', '[&!color=#000000]', '[&!color=#ff6288]']
	
	term_list = []

	vectors_by_tip = { }
	str_vectors_by_og = { }
	vectors_by_og = { }
	for l, line in enumerate(open('../master_go_presence.csv')):
		if(l != 0):
			line = line.split(',')
			og = line[1].strip(); tip = line[0].strip()
			if(og not in vectors_by_og):
				vectors_by_og.update({ og : [] })
				str_vectors_by_og.update({ og : { } })
					
			line_vec = [int(cell.strip()) for cell in line[2:] if cell.strip() != '']
			vectors_by_tip.update({ tip : [og, line_vec] })
			if(str(line_vec) not in str_vectors_by_og[og]):
				str_vectors_by_og[og].update({ str(line_vec) : len(vectors_by_og[og]) })
				vectors_by_og[og].append(line_vec)
		elif(l == 0):
			for cell in line.split(','):
				if(cell.strip() != '' and cell.strip() != 'Sequence' and cell.strip() != 'OG5'):
					term_list.append(cell.strip())
					
	for og in str_vectors_by_og:
		newick = ''
		colored_names = []
		for line in open('../FunctionSearching/' + og + '/' + og + '.tre', 'r'):
			temp = line.split(' ')[-1]
			if(temp.startswith('(') or temp.startswith('\t(')):
				newick = temp.split('\t')[-1].replace("'", '').replace('\\', '')
				line = line.split(',')
				for chunk in line:
					chunk = chunk.split('(')[-1].split(')')[0].split(':')[0]
									
					if(chunk in vectors_by_tip):
						#print(vectors_by_tip[chunk])
						colored_names.append(chunk + colors[str_vectors_by_og[og][str(vectors_by_tip[chunk][1])]])
					else:
						colored_names.append(chunk)
					
		write_nexus('../FunctionSearching/' + og + '/' + og + '_colored_by_vector.tre', colored_names, newick)					
		
	## REDUCED VERSION TAKING INTO ACCOUNT ESPECIALLY SIMILAR VECTORS ##
	
	replacement_pairs = { }
	reduced_vectors_by_og = { }
	for og in vectors_by_og:
		reduced_vectors_by_og.update({ og : [] }); replacement_pairs.update({ og : { } })
		dont_use = []
		for v1 in vectors_by_og[og]:
			if(v1 not in reduced_vectors_by_og[og] and v1 not in dont_use):
				hit = False
				matches = [v1]
				for v2 in vectors_by_og[og]:
					if(v1 != v2):
						consensus_list = []
						for i, val in enumerate(v1):
							if(val == v2[i]):
								consensus_list.append(1)
							else:
								consensus_list.append(0)
						
						if(sum(consensus_list) >= .9 * len(consensus_list)):
							matches.append(v2)
				
				reduced_vectors_by_og[og].append(v1)
				if(len(matches) > 0):
					for match in matches:
						replacement_pairs[og].update({ str(match) : str(v1) })
						dont_use.append(match)
					
											
	temp = reduced_vectors_by_og.copy(); reduced_vectors_by_og = { og : { } for og in temp}
	for og in temp:
		for i, vec in enumerate(temp[og]):
			reduced_vectors_by_og[og].update({ str(vec) : i })
									
	for og in vectors_by_og:
		newick = ''
		colored_names = []
		for line in open('../FunctionSearching/' + og + '/' + og + '.tre', 'r'):
			temp = line.split(' ')[-1]
			if(temp.startswith('(') or temp.startswith('\t(')):
				newick = temp.split('\t')[-1].replace("'", '').replace('\\', '')
				line = line.split(',')
				for chunk in line:
					chunk = chunk.split('(')[-1].split(')')[0].split(':')[0]
									
					if(chunk in vectors_by_tip):
						if(str(vectors_by_tip[chunk][1]) in replacement_pairs[og]):
							colored_names.append(chunk + colors[reduced_vectors_by_og[og][replacement_pairs[og][str(vectors_by_tip[chunk][1])]]])
						else:
							colored_names.append(chunk + colors[reduced_vectors_by_og[og][str(vectors_by_tip[chunk][1])]])
					else:
						colored_names.append(chunk)
						
		write_nexus('../FunctionSearching/' + og + '/' + og + '_colored_by_reduced_vectors.tre', colored_names, newick)
		
	with open('../unique_go_vectors.csv', 'w') as o:
		o.write('Gene Family,')
		for term in term_list:
			o.write(term + ',')
		o.write('\n')
		for og in vectors_by_og:
			for v, vec in enumerate(vectors_by_og[og]):
				o.write(og + '_' + str(v) + ',')
				for val in vec:
					o.write(str(val) + ',')
				o.write('\n')
			
		
	with open('../unique_go_vectors_reduced_by_similarity.csv', 'w') as o:
		o.write('Gene Family,')
		for term in term_list:
			o.write(term + ',')
		o.write('\n')
		for og in reduced_vectors_by_og:
			for v, vec in enumerate(reduced_vectors_by_og[og]):
				o.write(og + '_' + str(v) + ',')
				vec = vec.split('[')[-1].split(']')[0].split(',')
				for val in vec:
					o.write(str(val) + ',')
				o.write('\n')
		
	
def main():
	
	#input_dir = get_args()
	
	if(not os.path.isdir('../FunctionSearching')):
		os.system('mkdir ../FunctionSearching')
	
	#select_sequences(input_dir)
	
	#emap_selected_seqs()
	
	#master_GO_presence()
	
	map_unique_vectors_to_trees()
		
	#revigo_reduce_list(go_terms)

	
main()









