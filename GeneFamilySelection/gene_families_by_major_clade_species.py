#Dependencies
import os
import sys
import itertools
#For processing .fasta files
from Bio import SeqIO

#Name of directory containing ReadyToGo files to be processed
rtg_files_handle = '../../ReadyToGo_Files'
og_files_handle = '../../allOG5Files'
#Name of output spreadsheet
output_handle = "../species_per_major_clade_by_gene_family.csv"


def get_data_types():

	genomic = []
	transcriptomic = []
	
	for rtg_file in os.listdir(rtg_files_handle):
		if('ds_store' not in rtg_file.lower()):
			if('genom' in rtg_file.lower() or 'wgs' in rtg_file.lower() or 'wga' in rtg_file.lower()):
				genomic.append(rtg_file[:10])
			else:
				transcriptomic.append(rtg_file[:10])
				
	for og_file in os.listdir(og_files_handle):
		if('Store' not in og_file):
			for record in list(SeqIO.parse(og_files_handle + '/' + 	og_file, "fasta")):
				genomic.append(record.id[:10])
				
	return [list(dict.fromkeys(genomic)), list(dict.fromkeys(transcriptomic))]
	


def generate_dataframe():#genomic, transcriptomic):
	all_gene_families = dict.fromkeys(list(itertools.chain.from_iterable(['OG5_' + record.id.split('OG5_')[-1][:6] for record in list(SeqIO.parse(rtg_files_handle + '/' + 	rtg_file, "fasta"))] for rtg_file in os.listdir(rtg_files_handle) if('Store' not in rtg_file))))

	for gene_family in all_gene_families:	
		#all_gene_families[gene_family] = { 'am_tra' : 0, 'ba_tra' : 0, 'ee_tra' : 0, 'ex_tra' : 0, 'op_tra' : 0, 'pl_tra' : 0, 'sr_tra' : 0, 'za_tra' : 0, 'am_gen' : 0, 'ba_gen' : 0, 'ee_gen' : 0, 'ex_gen' : 0, 'op_gen' : 0, 'pl_gen' : 0, 'sr_gen' : 0, 'za_gen' : 0 }
		all_gene_families[gene_family] = { 'am' : 0, 'ba' : 0, 'ee' : 0, 'ex' : 0, 'op' : 0, 'pl' : 0, 'sr' : 0, 'za' : 0 }
	
	for rtg_file in os.listdir(rtg_files_handle):
		if('Store' not in rtg_file and (rtg_file[4] != 'b' or (rtg_file[4] == 'b' and (rtg_file[:2] == 'Ba' or rtg_file[:2] == 'Za')))):
					
			clade_code = rtg_file[0:2]
		
			#Collecting all of the gene family numbers in the current ReadyToGo file
			#ogs_in_file = list(dict.fromkeys([record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta"))]))
			ogs_in_file = ['OG5_' + record.id.split('OG5_')[-1][:6] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta"))]
			#Iterating all of the entries (and gene families) present in the file
			for og in ogs_in_file:
				try:
					#if(rtg_file[:10] in genomic):
						#all_gene_families[og][clade_code.lower() + '_gen'] += 1
					all_gene_families[og][clade_code.lower()] += 1
					#elif(rtg_file[:10] in transcriptomic):
						#all_gene_families[og][clade_code.lower() + '_tra'] += 1
						#all_gene_families[og][clade_code.lower()] += 1
				except KeyError:
					continue
					
	for og_file in os.listdir(og_files_handle):
		if('Store' not in og_file):
			taxa_in_file = list(dict.fromkeys([record.id[:10] for record in list(SeqIO.parse(og_files_handle + '/' + og_file, "fasta"))]))
			
			for taxon in taxa_in_file:
				try:
					#all_gene_families[og_file][taxon[:2].lower() + '_gen'] += 1
					all_gene_families[og_file][taxon[:2].lower()] += 1
				except KeyError:
					continue

	with open(output_handle, 'w') as o:
		for gene_family in all_gene_families:
			o.write(gene_family + ',')	
			for clade in all_gene_families[gene_family]:
				o.write(clade + ',')
			break
		o.write('\n')
				
		for gene_family in all_gene_families:
			o.write(gene_family + ',')	
			for clade in all_gene_families[gene_family]:
				o.write(str(all_gene_families[gene_family][clade]) + ',')
			o.write('\n')
			
			
def main():

	#genomic, transcriptomic = get_data_types()
		
	generate_dataframe()#genomic, transcriptomic)
	
	
main()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	