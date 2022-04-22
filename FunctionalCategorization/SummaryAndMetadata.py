#!/usr/bin/env python3
# coding=utf-8

#Author: Auden Cote-L'Heureux
#Contact at acotelheureux@smith.edu, audenemil@gmail.com, or on Slack if a Katzlabber.
#This script does a couple of things. It gets a bunch of GO term metadata from the EggNOG-mapper output, which is most of the script. 
#It also has a function that is a wrapper to create dataframes and call the R script to map functional metrics onto trees.
#Last update 05/20/21
#Can run with python 2 or 3

import os
import sys
import re
from Bio import SeqIO
import subprocess
import matplotlib
from tqdm import tqdm

					
def map_all_functional_metrics_to_tree_R_script():

	column_indices = [12, 13, 23, 11, 14, 15]
	file_prefixes = ['go', 'ec', 'pfam', 'preferred-name', 'ko', 'pathway']

	#Mapping all terms from EggNOG
	for og_idx, og in enumerate(os.listdir('../FunctionSearching')):
		if('ds_store' not in og.lower() and 1 == 0):
			for file in os.listdir('../FunctionSearching/' + og):
				if(file == 'eggnog_results.fa.emapper.annotations.tsv'):
					if(not os.path.isdir('../FunctionSearching/' + og + '/tree_mapping_dataframes')):
						os.mkdir('../FunctionSearching/' + og + '/tree_mapping_dataframes')
					if(not os.path.isdir('../FunctionSearching/' + og + '/mapped_figures')):
						os.mkdir('../FunctionSearching/' + og + '/mapped_figures')
						 
					print(str(og_idx) + '. ' + og)
					
					rel_lines = [[cell.strip() for cell in line.split('\t') if cell.strip() != ''] for line in open('../FunctionSearching/' + og + '/' + file) if line[0] != '#']
					
					info_dict = { line[0] : [[term for term in line[i].split(',') if float(len([l for l in rel_lines if term in l[i]])) / float(len(rel_lines)) >= 0] for i in column_indices] for line in rel_lines }
					
					for i, idx in enumerate(column_indices):
						all_terms = list(dict.fromkeys([term for seq in info_dict for term in info_dict[seq][i] if term != '-']))
						if(len(all_terms) >= 1):
							df_handle = '../FunctionSearching/' + og + '/tree_mapping_dataframes/'  + file_prefixes[i] + '.csv'
							with open(df_handle, 'w') as o:
								o.write('Seq,Term,Val\n')
								for seq in info_dict:
									for term in all_terms:
																										
										if(term in info_dict[seq][i]):
											o.write(seq + ',' + term + ',1\n')
										else:
											o.write(seq + ',' + term + ',0\n')
											
	#Mapping the unique PFams from HMMer
	for unique_file in tqdm(os.listdir('../Unique_Pfams')):
		if('PFam' in unique_file and unique_file.startswith('OG5_')):
			info_dict = { }
			
			for l, line in enumerate(open('../Unique_Pfams/' + unique_file)):
				if(l != 0):
					if(line.split('\t')[0].strip() not in info_dict):
						info_dict.update({ line.split('\t')[0].strip() : [] })
	
					info_dict[line.split('\t')[0].strip()].append(line.split('\t')[1].strip())
					
			all_pfams = list(dict.fromkeys([domain for seq in info_dict for domain in info_dict[seq]]))
			
			with open('../FunctionSearching/' + unique_file[:10] + '/tree_mapping_dataframes/unique_pfams.csv', 'w') as o:
				o.write('Seq,Term,Val\n')
				for seq in info_dict:
					for term in all_pfams:
						if(term in info_dict[seq]):
							o.write(seq + ',' + term + ',1\n')
						else:
							o.write(seq + ',' + term + ',0\n')
							
	#Mapping the GOs returned by PFam and the combined GOs (EggNOG + PFam)
	for go_file in os.listdir('../GOTermsFromPFams'):
		if(go_file.startswith('OG5_')):
			info_dict = { }
			
			os.system('rm ../FunctionSearching/' +	go_file.split('.tsv')[0] + '/tree_mapping_dataframes/gos_from_pfam.csv')
			os.system('rm ../FunctionSearching/' +	go_file.split('.tsv')[0] + '/tree_mapping_dataframes/combined_gos.csv')
			
			for line in open('../FunctionSearching/' + go_file.split('.tsv')[0] + '/eggnog_results.fa.emapper.annotations.tsv'):
				if(not line.startswith('#')):
					info_dict.update({ line.split('\t')[0].strip() : [] })
			
			for l, line in enumerate(open('../GOTermsFromPFams/' + go_file)):
				if(l != 0):
					info_dict[line.split('\t')[0].strip()] = [term.strip() for term in line.split('\t')[1].split(',') if term.strip() != '' and term.strip() != '-']
					
			all_terms = list(dict.fromkeys([term for seq in info_dict for term in info_dict[seq]]))
			
			if(len(all_terms) > 1):
				with open('../FunctionSearching/' +	go_file.split('.tsv')[0] + '/tree_mapping_dataframes/gos_from_pfam.csv', 'w') as o:	
					o.write('Seq,Term,Val\n')
					for seq in info_dict:
						for term in all_terms:
							if(term in info_dict[seq]):
								o.write(seq + ',' + term + ',1\n')
							else:
								o.write(seq + ',' + term + ',0\n')	
							
			for line in open('../FunctionSearching/' + go_file.split('.tsv')[0] + '/eggnog_results.fa.emapper.annotations.tsv'):
				if(not line.startswith('#')):
					for term in [term.strip() for term in line.split('\t')[12].split(',') if term.strip() != '' and term.strip() != '-']:
						info_dict[line.split('\t')[0].strip()].append(term)
									
			all_terms = list(dict.fromkeys([term for seq in info_dict for term in info_dict[seq]]))
			
			print(all_terms)
						
			with open('../FunctionSearching/' +	go_file.split('.tsv')[0] + '/tree_mapping_dataframes/combined_gos.csv', 'w') as o:
				o.write('Seq,Term,Val\n')
				for seq in info_dict:
					for term in all_terms:
						if(term in info_dict[seq]):
							o.write(seq + ',' + term + ',1\n')
						else:
							o.write(seq + ',' + term + ',0\n')
			
	os.system('/Library/Frameworks/R.framework/Versions/4.0/Resources/Rscript MapFunc2Tree.r')
	
	
def master_GO_presence():
	
	func_by_seq = { }
	
	for dir in os.listdir('../FunctionSearching'):
		if('OG5_' in dir):
			for line in open('../FunctionSearching/' + dir + '/' + 'eggnog_results.fa.emapper.annotations.tsv'):
				if(line[0] != '#'):
					line = line.split('\t')
					go_str = ';'.join(line[12].split(',')).replace('"', '').strip()					
					func_by_seq.update({ line[0].strip() : { 'OG5' : dir, 'GO_list' : [term.strip() for term in go_str.split(';') if term.strip() != '' and term.strip() != '-'] } })
			
			if(os.path.isfile('../GOTermsFromPFams/' + dir + '.tsv')):
				for l, line in enumerate(open('../GOTermsFromPFams/' + dir + '.tsv')):
					if(l != 0 and line.split('\t')[0].strip() in func_by_seq):
						func_by_seq[line.split('\t')[0].strip()]['GO_list'] += [term.strip() for term in line.split('\t')[1].split(',') if term.strip() != '' and term.strip() != '-']
						
	all_gos = []
	for seq in func_by_seq:
		for term in func_by_seq[seq]['GO_list']:
			if(term not in all_gos):
				all_gos.append(term)
				
	for seq in func_by_seq:
		func_by_seq[seq].update({ 'Original_GO_presence' : { term : 0 for term in all_gos } })
		for term in func_by_seq[seq]['GO_list']:
			func_by_seq[seq]['Original_GO_presence'][term] = 1
				
	with open('../master_go_presence.csv', 'w') as o:
		o.write('Sequence,OG5,')
		for term in all_gos:
			o.write(term + ',')
		o.write('\n')
		for seq in func_by_seq:
			o.write(seq + ',' + func_by_seq[seq]['OG5'] + ',')
			for term in all_gos:
				o.write(str(func_by_seq[seq]['Original_GO_presence'][term]) + ',')
			o.write('\n')
	
	
def vector_counts_by_tree():

	vectors_by_og = { }
	all_terms = []
	
	for l, line in enumerate(open('../master_go_presence.csv')):
		line = [cell.strip() for cell in line.split(',') if cell.strip() != '']
		if(line[1] != 'OG5_145229'):
			continue
		if(l == 0):
			all_terms = line[2:]
		else:
			if(line[1] not in vectors_by_og):
				vectors_by_og.update({ line[1] : { 'pro' : [], 'euk' : [], 'euk_counts' : { }, 'pro_counts' : { }, 'total_euks' : 0, 'total_pros' : 0 , 'euks_hit_bac' : 0, 'euks_hit_bac_w_go' : 0} })
				
			dom = 'pro'
			if(line[0][:2] != 'Ba' and line[0][:2] != 'Za'):
				dom = 'euk'
				if(tuple(line[2:]) not in vectors_by_og[line[1]]['euk_counts']):
					vectors_by_og[line[1]]['euk_counts'].update({ tuple(line[2:]) : 0 })
				vectors_by_og[line[1]]['euk_counts'][tuple(line[2:])] += 1
			else:
				if(tuple(line[2:]) not in vectors_by_og[line[1]]['pro_counts']):
					vectors_by_og[line[1]]['pro_counts'].update({ tuple(line[2:]) : 0 })
				vectors_by_og[line[1]]['pro_counts'][tuple(line[2:])] += 1
				
			if(line[2:] not in vectors_by_og[line[1]][dom]):
				vectors_by_og[line[1]][dom].append(line[2:])
				
			if(vectors_by_og[line[1]]['total_pros'] == 0 and vectors_by_og[line[1]]['total_euks'] == 0):
				for seq in open('../FunctionSearching/' + line[1] + '/eggnog_results.fa.emapper.annotations.tsv'):
					if(seq[:2] == 'Ba' or seq[:2] == 'Za'):
						vectors_by_og[line[1]]['total_pros'] += 1
					elif(seq[0] != '#'):
						vectors_by_og[line[1]]['total_euks'] += 1
						if('Bacteria' in seq.split('\t')[4]):
							vectors_by_og[line[1]]['euks_hit_bac'] += 1
							if('GO' in seq.split('\t')[12]):
								vectors_by_og[line[1]]['euks_hit_bac_w_go'] += 1
							
				
	vector_counts_by_og = { }
				
	for og in vectors_by_og:
		vector_counts_by_og.update({ og : { 'pro' : 0, 'euk' : 0, 'both' : 0, 'pro_props' : [], 'euk_props' : [], 'both_props' : [], 'prop_terms_shared_by_pros' : 0, 'prop_terms_shared_by_euks' : 0, 'prop_terms_shared_by_all' : 0} })
		
		bothed = []
		
		pro_model = []
		pro_sum = []
		if(len(vectors_by_og[og]['pro']) > 0):
			pro_model = [1 for val in next(iter(vectors_by_og[og]['pro']))]
			pro_sum = [0 for val in pro_model]
			all_model = [1 for val in pro_model]
			all_sum = [0 for val in pro_model]
			
		euk_model = []
		euk_sum = []
		if(len(vectors_by_og[og]['euk']) > 0):
			euk_model = [1 for val in next(iter(vectors_by_og[og]['euk']))]
			euk_sum = [0 for val in euk_model]
			all_model = [1 for val in euk_model]
			all_sum = [0 for val in euk_model]
				
		for vec in vectors_by_og[og]['pro']:
			if(vec not in vectors_by_og[og]['euk']):
				vector_counts_by_og[og]['pro'] += 1
				vector_counts_by_og[og]['pro_props'].append(round(float(vectors_by_og[og]['pro_counts'][tuple(vec)]) / float(vectors_by_og[og]['total_pros']), 3))
			elif(vec not in bothed):
				bothed.append(vec)
				vector_counts_by_og[og]['both'] += 1
				vector_counts_by_og[og]['both_props'].append(round(float(vectors_by_og[og]['pro_counts'][tuple(vec)] + vectors_by_og[og]['euk_counts'][tuple(vec)]) / float(vectors_by_og[og]['total_pros'] + vectors_by_og[og]['total_euks']), 3))
				
			for i, val in enumerate(pro_model):
				pro_model[i] *= int(vec[i])
				pro_sum[i] += int(vec[i])
				if(pro_sum[i] == 1):
					print(all_terms[i])

				all_model[i] *= int(vec[i])
				all_sum[i] += int(vec[i])
			
		for vec in vectors_by_og[og]['euk']:
			if(vec not in vectors_by_og[og]['pro']):
				vector_counts_by_og[og]['euk'] += 1
				vector_counts_by_og[og]['euk_props'].append(round(float(vectors_by_og[og]['euk_counts'][tuple(vec)]) / float(vectors_by_og[og]['total_euks']), 3))
			elif(vec not in bothed):
				bothed.append(vec)
				vector_counts_by_og[og]['both'] += 1
				vector_counts_by_og[og]['both_props'].append(round(float(vectors_by_og[og]['pro_counts'][tuple(vec)] + vectors_by_og[og]['euk_counts'][tuple(vec)]) / float(vectors_by_og[og]['total_pros'] + vectors_by_og[og]['total_euks']), 3))
				
			for i, val in enumerate(euk_model):
				euk_model[i] *= int(vec[i])
				euk_sum[i] += int(vec[i])
				
				all_model[i] *= int(vec[i])
				all_sum[i] += int(vec[i])
		
		pro_model = [val for i, val in enumerate(pro_model) if pro_sum[i] > 0]
		if(len(pro_model) > 0):
			vector_counts_by_og[og]['prop_terms_shared_by_pros'] = float(sum(pro_model))/float(len(pro_model))
		else:
			vector_counts_by_og[og]['prop_terms_shared_by_pros'] = ''
			
		euk_model = [val for i, val in enumerate(euk_model) if euk_sum[i] > 0]
		if(len(euk_model) > 0):
			vector_counts_by_og[og]['prop_terms_shared_by_euks'] = float(sum(euk_model))/float(len(euk_model))
		else:
			vector_counts_by_og[og]['prop_terms_shared_by_euks'] = ''
			
		all_model = [val for i, val in enumerate(all_model) if all_sum[i] > 0]
		if(len(all_model) > 0):
			vector_counts_by_og[og]['prop_terms_shared_by_all'] = float(sum(all_model))/float(len(all_model))
		else:
			vector_counts_by_og[og]['prop_terms_shared_by_all'] = ''
				
	with open('../vector_type_count_summary.csv', 'w') as o:
		o.write('Gene Family,Proportion of Euks that Hit Bacteria,and GO,Prokaryote Only,Prokaryote Proportions,Prop Terms Shared by Prokaryotes,Eukaryote Only,Eukaryote Proportions,Prop Terms Shared by Eukaryotes,Both,Both Proportions,Prop Terms shared by All Tips\n')
		for og in vector_counts_by_og:
			o.write(og + ',' + str(float(vectors_by_og[og]['euks_hit_bac'])/float(vectors_by_og[og]['total_euks'])) + ',' + str(float(vectors_by_og[og]['euks_hit_bac_w_go'])/float(vectors_by_og[og]['total_euks'])) + ',' + str(vector_counts_by_og[og]['pro']) + ',' + ' '.join([str(val * 100) + '%' for val in sorted(vector_counts_by_og[og]['pro_props'], reverse = True)]) + ',' + str(vector_counts_by_og[og]['prop_terms_shared_by_pros']) + ',' + str(vector_counts_by_og[og]['euk']) + ',' + ' '.join([str(val * 100) + '%' for val in sorted(vector_counts_by_og[og]['euk_props'], reverse = True)]) + ',' + str(vector_counts_by_og[og]['prop_terms_shared_by_euks']) + ',' + str(vector_counts_by_og[og]['both']) + ',' + ' '.join([str(val * 100) + '%' for val in sorted(vector_counts_by_og[og]['both_props'], reverse = True)]) + ',' + str(vector_counts_by_og[og]['prop_terms_shared_by_all']) + '\n')
			
			
def go_term_metadata():

	og_groups = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('og_groups.csv') }
	
	prev_good_op = ['OG5_130801', 'OG5_132746', 'OG5_133998', 'OG5_138333', 'OG5_130277', 'OG5_131043', 'OG5_131262', 'OG5_131754', 'OG5_130469', 'OG5_132932', 'OG5_135659', 'OG5_135811', 'OG5_137943', 'OG5_155629', 'OG5_146297', 'OG5_153460']

	all_trees = [og_dir for og_dir in os.listdir('../FunctionSearching')]
	trees_with_original = []
	trees_with_pfam = []
	trees_only_pfam = []
	trees_only_original = []
	trees_with_any_go_terms = []
	
	above_10_pfam = [];
	above_10_combined = [];
	above_10_original = [];
	
	go_terms_pfam = []
	go_terms_original = []
	
	added_go_terms = []
	redundant_go_terms = []
	go_terms_unique_to_original = []
	
	num_terms_pfam = { }
	num_terms_original = { }
	
	for og_dir in tqdm(os.listdir('../FunctionSearching')):
	
		if('ds_store' in og_dir.lower()):
			continue
		try:
			recip = og_groups[og_dir]
		except KeyError:
			continue
		if(recip == 'EGT'):
			recips = ['Pl', 'Sr', 'EE']
		elif(recip == 'Am_Ex'):
			recips = ['Am', 'Ex', 'Sr', 'EE']
		elif(len(recip) == 2):
			recips = [recip]
		else:
			recips = []
			continue
			
		if(recip != 'Am_Ex'):
			continue
			
		# if(os.path.isfile('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/pfam.csv')):
# 			trees_with_any_go_terms.append(og_dir)
#  		continue
	
		if(os.path.isfile('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/go.csv')):
			trees_with_original.append(og_dir)
			
			with_go = []; tot = []
			for seq in open('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/go.csv'):
				if(seq.split(',')[0].strip() not in num_terms_original):
					num_terms_original.update({ seq.split(',')[0].strip() : [] })
				
				if(seq[:2] in recips and seq.split(',')[2].strip() == '1'):
					with_go.append(seq.split(',')[0].strip())
					num_terms_original[seq.split(',')[0].strip()].append(seq.split(',')[1].strip())
				if(seq[:2] in recips):
					tot.append(seq.split(',')[0].strip())
			if(float(len((list(dict.fromkeys(with_go)))))/float(len(list(dict.fromkeys(tot)))) >= .1):
				above_10_original.append(og_dir)
			
		actual_pfam_terms = []
		if(os.path.isfile('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/gos_from_pfam.csv')):
			trees_with_pfam.append(og_dir)
			if(og_dir not in trees_with_original):
				trees_only_pfam.append(og_dir)
				
			with_go = []; tot = []
			for seq in open('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/gos_from_pfam.csv'):
				if(seq.split(',')[0].strip() not in num_terms_pfam):
					num_terms_pfam.update({ seq.split(',')[0].strip() : [] })
					
				if(seq[:2] in recips and seq.split(',')[2].strip() == '1'):
					with_go.append(seq.split(',')[0].strip())
					num_terms_pfam[seq.split(',')[0].strip()].append(seq.split(',')[1].strip())
				if(seq[:2] in recips):
					tot.append(seq.split(',')[0].strip())
			if(float(len((list(dict.fromkeys(with_go)))))/float(len(list(dict.fromkeys(tot)))) >= .1):
				above_10_pfam.append(og_dir)
			
		if(os.path.isfile('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/combined_gos.csv')):
			try:
				if(len([ln for ln in open('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/combined_gos.csv')]) == 1):	
					continue
			except:
				pass
			trees_with_any_go_terms.append(og_dir)
			if(og_dir not in trees_with_pfam):
				trees_only_original.append(og_dir)
				
			with_go = []; tot = []
			for seq in open('../FunctionSearching/' + og_dir + '/tree_mapping_dataframes/combined_gos.csv'):
				if(seq[:2] in recips and seq.split(',')[2].strip() == '1'):
					with_go.append(seq.split(',')[0].strip())
				if(seq[:2] in recips):
					tot.append(seq.split(',')[0].strip())
			if(float(len((list(dict.fromkeys(with_go)))))/float(len(list(dict.fromkeys(tot)))) >= .1):
				above_10_combined.append(og_dir)
				
		original_terms = [term for seq in num_terms_original for term in num_terms_original[seq] if og_dir in seq]
		pfam_terms = [term for seq in num_terms_pfam for term in num_terms_pfam[seq] if og_dir in seq]
		
		if(og_dir not in trees_only_pfam and og_dir not in trees_only_original):
			added_go_terms.append(len(list(dict.fromkeys([seq for seq in pfam_terms if seq not in original_terms]))))
			redundant_go_terms.append(len(list(dict.fromkeys([seq for seq in pfam_terms if seq in original_terms]))))
			go_terms_unique_to_original.append(len(list(dict.fromkeys([seq for seq in original_terms if seq not in pfam_terms]))))
			
	average_added_go_terms = float(sum(added_go_terms)) / float(len(added_go_terms))
	average_redundant_go_terms = float(sum(redundant_go_terms)) / float(len(redundant_go_terms))
	average_go_terms_unique_to_original = float(sum(go_terms_unique_to_original)) / float(len(go_terms_unique_to_original))
		
	original_avg_terms = float(sum([len(num_terms_original[seq]) for seq in num_terms_original]))/float(len([seq for seq in num_terms_original if len(num_terms_original[seq]) > 0]))
	pfam_avg_terms = float(sum([len(num_terms_pfam[seq]) for seq in num_terms_pfam]))/float(len([seq for seq in num_terms_pfam if len(num_terms_pfam[seq]) > 0]))
	
	for og in trees_with_original:
		if(og in trees_with_pfam):
			print(og)
							
	with open('../pfam_addition_metadata.txt', 'w') as o:
		o.write('Total number of trees: ' + str(len(all_trees)) + '\n')
		o.write('Total number of trees that originally returned GO terms: ' + str(len(trees_with_original)) + '\n')
		o.write('Total number of trees with PFam-associated GO terms: ' + str(len(trees_with_pfam)) + '\n')
		o.write('Total number of trees with ONLY PFam-associated GO terms: ' + str(len(trees_only_pfam)) + '\n')
		o.write('Total number of trees with ONLY originally returned GO terms: ' + str(len(trees_only_original)) + '\n\n')
		
		o.write('Number of Opisthokont trees with at least 10% of recipient sequences that returned GO terms originally: ' + str(len(above_10_original)) + '\n')
		o.write('Number of Opisthokont trees with at least 10% of recipient sequences that returned PFam-associated GO terms: ' + str(len(above_10_pfam)) + '\n')
		o.write('Number of Opisthokont trees with at least 10% of recipient sequences that returned GO terms, combined: ' + str(len(above_10_combined)) + '\n\n')
		
		o.write('Total number of sequences that returned GO terms from EggNOG: ' + str(len([seq for seq in num_terms_original if len(num_terms_original[seq]) > 0])) + '\n')
		o.write('Total number of sequences that returned GO terms from PFam: ' + str(len([seq for seq in num_terms_pfam if len(num_terms_pfam[seq]) > 0])) + '\n')
		o.write('Average number of GO terms returned by EggNOG : ' + str(original_avg_terms) + '\n')
		o.write('Average number of GO terms returned by PFam : ' + str(pfam_avg_terms) + '\n\n')
		
		o.write('Average number of GO terms added to tree by PFam: ' + str(average_added_go_terms) + '\n')
		o.write('Average number of GO terms redundant in tree and returned by PFam: ' + str(average_redundant_go_terms) + '\n')
		o.write('Average number of GO terms returned only by EggNOG in tree: ' + str(average_go_terms_unique_to_original) + '\n')
	
					
	
def main():
		
	if(not os.path.isdir('../FunctionSearching')):
		os.system('mkdir ../FunctionSearching')
	
	map_all_functional_metrics_to_tree_R_script()
	
	master_GO_presence()
		
	vector_counts_by_tree()
	
	go_term_metadata()
		
	
main()
