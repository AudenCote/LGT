import os
import sys
import re
from Bio import SeqIO
from Bio.Blast.NCBIWWW import qblast
import numpy as np
from subprocess import check_output


rtg_files_handle = '/Volumes/GoogleDrive/.shortcut-targets-by-id/19d0BR7kBq_DD_jz7mKcKd2U4lFy_2zcD/Katzlab_Fulltimers_Only/ReadyToGo_NTD_All'


def get_args():

	fixed_level = -1
	if('--fixed_level' in sys.argv):
		fixed_level = int(sys.argv[sys.argv.index('--fixed_level') + 1])
		
	return fixed_level
	
	
def is_bacbin(tdc):

	if((tdc[:2] != 'Ba' and tdc[:2] != 'Za' and tdc[4] != 'b') or (tdc[:2] == 'Ba' or tdc[:2] == 'Za')):
		return False
	else:
		return True
	

def cluster(id, cluster_size, seqs, to_cluster_handle, clustered_handle, representative_seqs_handle, restrict_by_major_clade):
	
	with open(to_cluster_handle, 'w') as o:		
		for seq in seqs:
			o.write('>' + seq + '\n' + seqs[seq] + '\n')
						
	os.system('vsearch --cluster_smallmem ' + to_cluster_handle + ' --usersort --uc ' + clustered_handle + ' --id ' + str(id))
		
	clusters_to_keep = []
	for line in open(clustered_handle, 'r'):
		line = line.split('\t')
		if(line[0] == 'C'):
			if(int(line[2]) >= cluster_size):
				clusters_to_keep.append([line[1], line[2], line[8][:2]])
				
	final_clusters = []
	if(restrict_by_major_clade):
		major_clades = []
		for cluster in sorted(clusters_to_keep, key = lambda x : x[1], reverse = True):
			if(cluster[2] not in major_clades):
				print(cluster)
				major_clades.append(cluster[2])
				final_clusters.append(cluster[0])
	else:
		final_clusters = [c[0] for c in clusters_to_keep]
				
	representative_seqs = { }
	for line in open(clustered_handle, 'r'):
		line = line.split('\t')
		if(line[0] == 'H' and line[1] in final_clusters):
			if(line[9].split('\n')[0] not in representative_seqs):
				representative_seqs.update({ line[9].split('\n')[0] : seqs[line[9].split('\n')[0]] })
		elif(line[0] == 'S' and line[1] in final_clusters and line[8].split('\n')[0] not in representative_seqs):
			representative_seqs.update({ line[8].split('\n')[0] : seqs[line[8].split('\n')[0]] })
					
	with open(representative_seqs_handle, 'w') as o:
		for seq in representative_seqs:
			o.write('>' + seq + '\n' + representative_seqs[seq] + '\n\n')	
	
def select_taxa():

	for pg_file in os.listdir('../Subtrees/Preguidance_Filtered'):
		og = 'OG5_' + pg_file.split('OG5_')[-1][:6]
		
		all_seqs = { }; length_tuples = [];
		records = [record for record in list(SeqIO.parse('../Subtrees/Preguidance_Filtered/' + pg_file, 'fasta'))]
		
		for record in records:
			all_seqs.update({ record.id : str(record.seq) })
			length_tuples.append((record.id, len(str(record.seq))))
						
		final_sorted = { }; transcriptomic = [];
		for seq in sorted(length_tuples, reverse=True, key=lambda x: x[1]):
			if(not is_bacbin(seq[0][:10])):
				if('Cov' in seq[0][10:]):
					coverage = int(seq[0][10:].split('Cov')[-1].split('_')[0].split('.')[0])
					transcriptomic.append((seq[0], coverage))
				else:
					final_sorted.update({ seq[0] : all_seqs[seq[0]] })
				
		for cov_record in sorted(transcriptomic, reverse=True, key=lambda x: x[1]):
			final_sorted.update({ cov_record[0] : all_seqs[cov_record[0]] })
		
		cluster(.8, 3, final_sorted, '../Subtrees/Clustering/' + og + '_tree_seqs_to_cluster.fasta', '../Subtrees/Clustering/' + og + '_tree_seqs_clustered.uc', '../Subtrees/BLAST/' + og + '_to_blast_from_subtree.fasta', True)
		
def blast_seqs_from_tree():

	for tb_file in os.listdir('../Subtrees/BLAST'):
		if(tb_file.endswith('_to_blast_from_subtree.fasta')):
			og = tb_file[:10]
			records = list(SeqIO.parse('../Subtrees/BLAST/' + tb_file, 'fasta'))
			
			blast_hits = []
			for r, record in enumerate(records):
				print('\nBLAST-ing representative sequence ' + record.id + ' from cluster ' + str(r) + '\n')
				blast_hits.append(qblast(program = 'blastp', sequence = str(record.seq), database = 'nr', hitlist_size = 100))
 						
			with open('../Subtrees/BLAST/' + og + '_blast_output.xml', "w") as o:
				for bh in blast_hits:
					o.write(bh.read())
					bh.close()
					
					
#This function queries Entrez Search with the genus and species name and returns the taxonomy for each name if available, used for the taxonomy_by_id_cov() function
def get_taxonomy(taxa):

	taxonomies = {}
	
	for taxon in taxa:	
		print('\nFetching taxonomy for taxon ' + taxon + '\n')
		if(taxon != ''):
			try:
				output = str(check_output('esearch -db taxonomy -query "' + taxon + '" | efetch -format xml', shell = True, executable='/bin/bash'))
				#sleep(2)
		
				taxonomy_strings = []
		
				for line in output.split('<'):
					if('lineage' in line.lower() and 'lineageex' not in line.lower()):
						if('cellular organisms' in line.split('>')[1].split('\n')[0]):
							taxonomy_strings.append(line.split('>')[1].split('\n')[0])
						
			except:
				continue

			#If multiple taxonomies were returned
			if(len(taxonomy_strings) > 1):
				final_string = 'MTR'
			elif(len(taxonomy_strings) == 0):
				final_string = 'NTR'
			else:
				final_string = taxonomy_strings[0]
															
			taxonomies.update({ taxon : final_string })
    
	return taxonomies
					
					
def get_hit_info(bh_file):
	
	hit_info = []
	if(bh_file.endswith('_blast_output.xml')):
		lines = [line for line in open('../Subtrees/BLAST/' + bh_file)]
		query_len = 0;
		for l, line in enumerate(lines):
			#Written like this in case more information from iteration is needed
			if('<Iteration>' in line):
				i = l
				while '<Iteration_hits>' not in lines[i]:
					if('<Iteration_query-len>' in lines[i]):
						query_len = int(lines[i].split('>')[1].split('<')[0])
					
					i += 1
					
			#For each hit
			if('<Hit>' in line):
				id = 0; query_start = 0; query_end = 0; taxon_name = ''; align_len = 0; seq = '';
				i = l
				while '</Hit>' not in lines[i]:
					if('<Hit_def>' in lines[i]):
						try:
							taxon_name = lines[i].split('[')[1].split(']')[0]
							if(len(taxon_name.split(' ')) > 1):
								if(taxon_name.split(' ')[1].strip() == 'sp' or taxon_name.split(' ')[1].strip() == 'sp.'):
									taxon_name = taxon_name.split(' ')[0]	
								else:
									taxon_name = taxon_name.split(' ')[0] + ' ' + taxon_name.split(' ')[1]
						except IndexError:
							taxon_name = 'IndexError'
					elif('<Hsp_query-from>' in lines[i]):
						query_start = int(lines[i].split('>')[1].split('<')[0])
					elif('<Hsp_query-to>' in lines[i]):
						query_end = int(lines[i].split('>')[1].split('<')[0])
					elif('<Hsp_identity>' in lines[i]):
						id = int(lines[i].split('>')[1].split('<')[0])
					elif('<Hsp_align-len>' in lines[i]):
						align_len = int(lines[i].split('>')[1].split('<')[0])
					elif('<Hsp_hseq>' in lines[i]):
						seq = lines[i].split('>')[1].split('<')[0]
							
					i += 1
				
				if(taxon_name != 'IndexError'):
					#Record information about the hit					
					hit_info.append({ 'hit_name' : taxon_name, 'coverage' : float(query_end - query_start - 1) / float(query_len), 'identity' : float(id) / float(align_len), 'sequence' : str(seq) })
				
	return hit_info
						
						
def taxonomy_by_id_cov(bh_file, fixed_level):
	if(not os.path.isdir('../Subtrees/TaxIDCovSheets')):
		os.system('mkdir ../Subtrees/TaxIDCovSheets')

	og = bh_file [:10]
	
	num_bars = 20
	mod_num = 100/num_bars
	
	hit_info = get_hit_info(bh_file)
	
	taxa = list(dict.fromkeys([hit['hit_name'] for hit in hit_info]))
	
	if(not os.path.isfile('../blast_hit_taxonomy.csv')):
		taxonomies = get_taxonomy(taxa)
	
		with open('../blast_hit_taxonomy.csv', 'w') as o:
			for tax in taxonomies:
				o.write(tax + ',' + taxonomies[tax] + '\n')
	else:
		taxonomies = { line.split(',')[0] : line.split(',')[1].split('\n')[0] for line in open('../blast_hit_taxonomy.csv') if('MTR' not in line.split(',')[1].split('\n')[0] and 'NTR' not in line.split(',')[1].split('\n')[0])}
		
	for hit in hit_info[::-1]:
		if(hit['hit_name'] not in taxonomies):
			hit_info.remove(hit)
	
	ids = [hit['identity'] for hit in hit_info]
	covs = [hit['coverage'] for hit in hit_info]
	
	id_percentiles = { }
	cov_percentiles = { }
	hits_per_id_percentile = { }
	hits_per_cov_percentile = { }
		
	for i in range(101):
		if(i % mod_num == 0):
			id_percentiles.update({ i : np.percentile(ids, i) }); cov_percentiles.update({ i : np.percentile(covs, i) })
			hits_per_id_percentile.update({ i : [] })
			hits_per_cov_percentile.update({ i : [] })
	
	for hit in hit_info:
		for i in range((100 - mod_num) + 1):
			if(i % mod_num == 0):
				if(hit['identity'] > id_percentiles[i] and hit['identity'] <= id_percentiles[i + mod_num]):
					hits_per_id_percentile[i].append(hit)
					
				if(hit['coverage'] > cov_percentiles[i] and hit['coverage'] <= cov_percentiles[i + mod_num]):
					hits_per_cov_percentile[i].append(hit)
	
	levels_per_id_perc = { }
	levels_per_cov_perc = { }
			
	for i in range(101 - mod_num):
		if(i % mod_num == 0):
			perc_tax_id = [taxonomies[hit['hit_name']].strip().split(';') for hit in hits_per_id_percentile[i]]
			perc_tax_cov = [taxonomies[hit['hit_name']].strip().split(';') for hit in hits_per_cov_percentile[i]]	
			perc_tax_both = [perc_tax_id, perc_tax_cov]
			for p, perc_tax in enumerate(perc_tax_both):
				for t, tax_list in enumerate(perc_tax):
					for tax_level in perc_tax_both[p][t][::-1]:
						if('cellular' in tax_level.lower() or 'eukaryot' in tax_level.lower() or tax_level.lower() == '' or 'incertae sedis' in tax_level.lower() or 'unclassified' in tax_level.lower() or 'sample' in tax_level.lower()):
							perc_tax_both[p][t].remove(tax_level)
			
			for p, perc_tax in enumerate(perc_tax_both):
				consensus_levels = { }
				for t, tax_list in enumerate(perc_tax):
					for l, level in enumerate(tax_list):
						if(l not in consensus_levels):
							consensus_levels.update({ l : {} })
						if(level not in consensus_levels[l]):
							consensus_levels[l].update({ level : 0 })
							
						consensus_levels[l][level] += 1
						
				level_of_difference = -1
				ratio_list = []; level_list = []
				for l, level in enumerate(consensus_levels):
					total = 0
					for count in consensus_levels[level]:
						total += consensus_levels[level][count] + 1
					
					ratios = []
					for count in consensus_levels[level]:
						ratio = float(consensus_levels[level][count] + 1)/float(total)
						if(ratio > .5 and ratio <= .95 and fixed_level == -1):
							ratios = [float(consensus_levels[level][count] + 1)/float(total) for count in consensus_levels[level]]
							level_list = consensus_levels[level]
							level_of_difference = l
							break
							
					if(l == fixed_level):
						ratios = [float(consensus_levels[level][count] + 1)/float(total) for count in consensus_levels[level]]
						level_list = consensus_levels[level]
						level_of_difference = l
					
					if(level_of_difference > -1):
						ratio_list = ratios
						break
											
				if(p == 0):
					levels_per_id_perc.update({ i : [level_of_difference, ratio_list, level_list] })
				elif(p == 1):
					levels_per_cov_perc.update({ i : [level_of_difference, ratio_list, level_list] })
	
	inclusive_df_id = { }
	inclusive_df_cov = { }
	for perc in levels_per_id_perc:
		inclusive_df_id.update({ perc : [] })
		inclusive_df_cov.update({ perc : [] })
					
	for p1 in levels_per_id_perc:
		for level in levels_per_id_perc[p1][2]:
			for p2 in inclusive_df_id:
				if(level not in [tup[0] for tup in inclusive_df_id[p2]]):
					inclusive_df_id[p2].append((level, 0))
					
	for p1 in levels_per_cov_perc:
		for level in levels_per_cov_perc[p1][2]:
			for p2 in inclusive_df_cov:
				if(level not in [tup[0] for tup in inclusive_df_cov[p2]]):
					inclusive_df_cov[p2].append((level, 0))
	
	for p1 in inclusive_df_id:
		for p2 in levels_per_id_perc:
			if(p1 == p2):
				for level in levels_per_id_perc[p1][2]:
					if(level in [tup[0] for tup in inclusive_df_id[p1]]):
						for tup in inclusive_df_id[p1][::-1]:
							if(tup[0] == level):
								inclusive_df_id[p1].remove(tup)
					inclusive_df_id[p1].append((level, levels_per_id_perc[p1][2][level]))
					
	for p1 in inclusive_df_cov:
		for p2 in levels_per_cov_perc:
			if(p1 == p2):
				for level in levels_per_cov_perc[p1][2]:
					if(level in [tup[0] for tup in inclusive_df_cov[p1]]):
						for tup in inclusive_df_cov[p1][::-1]:
							if(tup[0] == level):
								inclusive_df_cov[p1].remove(tup)
					inclusive_df_cov[p1].append((level, levels_per_cov_perc[p1][2][level]))
					
	if(fixed_level != -1):
		id_handle = '../Subtrees/TaxIDCovSheets/' + og + '_taxonomic_level_' + str(levels_per_cov_perc[0][0]) + '_by_id.csv'
		cov_handle = '../Subtrees/TaxIDCovSheets/' + og + '_taxonomic_level_' + str(levels_per_cov_perc[0][0]) + '_by_cov.csv'
	else:
		id_handle = '../Subtrees/TaxIDCovSheets/' + og + '_taxonomic_level_by_id.csv'
		cov_handle = '../Subtrees/TaxIDCovSheets/' + og + '_taxonomic_level_by_cov.csv'
					
	with open(id_handle, 'w') as o:
		o.write('Percentile,ID,Level_Num,Tax_Level,Count\n')
		for percentile in inclusive_df_id:
			for level in sorted(inclusive_df_id[percentile], key=lambda x: x[0]):
				if(float(sum([lev[1] for pe in inclusive_df_id for lev in inclusive_df_id[pe] if(lev[0] == level[0])]))/float(len(inclusive_df_id[percentile])) > .5):
					o.write(str(percentile) + ',' + str(id_percentiles[percentile]) + ',' + str(levels_per_id_perc[percentile][0]) + ',' + level[0] + ',' + str(level[1]) + '\n')
			
	with open(cov_handle, 'w') as o:
		o.write('Percentile,Coverage,Level_Num,Tax_Level,Count\n')
		for percentile in inclusive_df_cov:
			for level in sorted(inclusive_df_cov[percentile], key=lambda x: x[0]):
				if(float(sum([lev[1] for pe in inclusive_df_cov for lev in inclusive_df_cov[pe] if(lev[0] == level[0])]))/float(len(inclusive_df_cov[percentile])) > .5):
					o.write(str(percentile) + ',' + str(cov_percentiles[percentile]) + ',' + str(levels_per_cov_perc[percentile][0]) + ',' + level[0] + ',' + str(level[1]) + '\n')


def cluster_hits():
	
	for bh_file in os.listdir('../Subtrees/BLAST'):
		if(bh_file.endswith('blast_output.xml')):
			og = bh_file[:10]
		
			bh_seqs = [(hit['hit_name'], hit['sequence'], hit['identity']) for hit in get_hit_info(bh_file) if(float(hit['identity']) >= .8)]
			
			final_sorted = { }
			included_specs = []
			maj_counts = { 'Ba' : 0, 'Za' : 0, 'Op' : 0, 'Ex' : 0, 'Am' : 0, 'Pl' : 0, 'Sr' : 0, 'EE' : 0 }
			for s, seq in enumerate(sorted(bh_seqs, reverse=True, key=lambda x: x[2])):
				if(not ('_'.join(seq[0].split(' ')) + '_' + str(s)).startswith('_')):
					if(seq[0] not in included_specs and 'protein' not in seq[0]):
						taxonomy = get_taxonomy([seq[0]])[seq[0]].lower()
						if(taxonomy != 'mtr' and taxonomy != 'ntr'):
						
							if('bacter' in taxonomy):
								maj = 'Ba'
							elif('archaea' in taxonomy):
								maj = 'Za'
							elif('opistho' in taxonomy):
								maj = 'Op'
							elif('excavat' in taxonomy):
								maj = 'Ex'
							elif('amoeb' in taxonomy):
								maj = 'Am'
							elif('archaeplastid' in taxonomy):
								maj = 'Pl'
							elif('sar' in taxonomy):
								maj = 'Sr'
							else:
								maj = 'EE'
															
							if(maj_counts[maj] <= 9):
								maj_counts[maj] += 1
								print(maj + '_nw_' + '_'.join(seq[0].split(' ')) + '_' + str(s))
								final_sorted.update({ maj + '_nw_' + '_'.join(seq[0].split(' ')) + '_' + str(s) : seq[1].replace('-', '') })
							
							included_specs.append(seq[0])
						
			rep_seqs_handle = '../Subtrees/Preguidance_Filtered/' + og + '_preguidance_with_representative_seqs.fasta'
			with open(rep_seqs_handle, 'w') as o:
				for seq in final_sorted:
					o.write('>' + seq + '\n' + final_sorted[seq] + '\n')
			#to_cluster_handle = '../Subtrees/Clustering/' + og + '_blast_hits_to_cluster.fasta'		
			#cluster(.99, 1, final_sorted, to_cluster_handle, '../Subtrees/Clustering/' + og + '_blast_hits_clustered.uc', rep_seqs_handle, False)
			
			rep_seqs = list(SeqIO.parse(rep_seqs_handle, 'fasta'))
			rep_ids = [record.id for record in rep_seqs]
			try:
				prev_seqs = list(SeqIO.parse('../Subtrees/Preguidance_Filtered/' + og + '_preguidance_filtered.fasta', 'fasta'))
			except:
				print('\nA pre-Guidance file could not be found for OG ' + og + '\n')
				continue
			
			seqs_to_write = prev_seqs
			seqs_to_write.extend(rep_seqs)
					
			with open(rep_seqs_handle, 'w') as o:
				for r, record in enumerate(seqs_to_write):
					if(record.id in rep_ids):
						o.write('>' + record.id[:6] + ("%04d" % r) + '_' + record.id[6:] + '_' + og + '\n' + str(record.seq) + '\n\n')
					else:
						o.write('>' + record.id + '\n' + str(record.seq) + '\n\n')
						
						
def update_taxon_list():

	already_written = []
	with open('../taxa_inclusive.txt', 'w') as o1:
		with open('../taxa_inclusive_database.txt', 'w') as o2:
			for file in os.listdir('../Subtrees/Preguidance_Filtered/'):
				if('_preguidance_with_representative_seqs' in file):
					for seq in list(SeqIO.parse('../Subtrees/Preguidance_Filtered/' + file, 'fasta')):
						if(seq.id[:10] not in already_written):
							o1.write(seq.id[:10] + '\n')
							o2.write('nw,euk,' + seq.id[:10] + '\n')
							already_written.append(seq.id[:10])
			
						

def main():

	fixed_level = get_args()

	if(not os.path.isdir('../Subtrees')):
		print('\nCould not find the "Subtrees" directory -- make sure that you run the script 0_GetSubtrees.py before running this one (1_ReferenceSeqs.py\n')
		exit()
		
	if(not os.path.isdir('../Subtrees/Clustering')):
		os.system('mkdir ../Subtrees/Clustering')
		
	if(not os.path.isdir('../Subtrees/BLAST')):
		os.system('mkdir ../Subtrees/BLAST')
	
	select_taxa()
	
	blast_seqs_from_tree()
	
	for bh_file in os.listdir('../Subtrees/BLAST'):
		if(bh_file.endswith('blast_output.xml')):
			break
			#taxonomy_by_id_cov(bh_file, fixed_level)
	
	cluster_hits()
	
	update_taxon_list()
		
	
main()
	