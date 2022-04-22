import os
import sys
import itertools
#For processing .fasta files
from Bio import SeqIO


#Name of directory containing ReadyToGo files to be processed
rtg_files_handle = '../../ReadyToGo_Files'
og_files_handle = '../../allOG5Files'
#Name of output spreadsheet
output_handle = '../minor_clade_presence.csv'


def get_species_by_genus(bacbin):
	
	pass


def generate_dataframe(bacbin, relative):

	species_by_genus = get_species_by_genus(bacbin)

	genera = {}					
	
	if(bacbin):
		taxa_per_og = dict.fromkeys(list(itertools.chain.from_iterable([record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta")) if 'OG5_' in record.id[-10:]] for rtg_file in os.listdir(rtg_files_handle) if((rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b') or rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'))))
	elif(not bacbin):
		taxa_per_og = dict.fromkeys(list(itertools.chain.from_iterable([record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta")) if 'OG5_' in record.id[-10:]] for rtg_file in os.listdir(rtg_files_handle) if('ds_store' not in rtg_file.lower()))))
	
	for gene_family in taxa_per_og:
		taxa_per_og[gene_family] = clade_codes.copy()
			
	if(bacbin):
		for rtg_file in os.listdir(rtg_files_handle):
			if((rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b') or rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'):
		
				if(rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'):
					clade_code = rtg_file[:2].lower()
				elif(rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b'):
					clade_code = rtg_file[:5].lower()
					
				all_ogs = [record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta"))]
				ogs_in_file = list(dict.fromkeys(all_ogs))

				#For each gene family present in the cell (file)
				for og in ogs_in_file:
					taxa_per_og[og][clade_code] += 1
	
	elif(not bacbin):
		with open('crap.txt', 'w') as o:
			for rtg_file in os.listdir(rtg_files_handle):
				if('ds_store' not in rtg_file.lower()):
					clade_code = rtg_file[:5].lower()
					
					all_ogs = [record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta"))]
					ogs_in_file = list(dict.fromkeys(all_ogs))
	
					#For each gene family present in the cell (file)
					for og in ogs_in_file:
						if('OG5' in og):
							taxa_per_og[og][clade_code] += 1
																	
	with open(output_handle, 'w') as o:	
		o.write(',')
		for clade in taxa_per_og['OG5_128717']:
			o.write(clade + ',')
		o.write('\n')
		for og in taxa_per_og:
			o.write(og + ',')
			for clade in taxa_per_og[og]:
				if(bacbin and relative):
					if(clade != 'ba' and clade != 'za'):
						o.write(str(float(taxa_per_og[og][clade]) / float(minor_clade_counts[clade])) + ',')
					else:
						o.write(str(taxa_per_og[og][clade]) + ',')
				elif(not bacbin and relative):
					if(clade[:2] != 'ba' and clade[:2] != 'za'):
						o.write(str(float(taxa_per_og[og][clade]) / float(minor_clade_counts[clade])) + ',')
					else:
						o.write(str(taxa_per_og[og][clade]) + ',')
				elif(not relative):
					o.write(str(taxa_per_og[og][clade]) + ',')
				
			o.write('\n')
	
def main():

	bacbin = False
	relative = False
	
	generate_dataframe(bacbin, relative)
		
	
main()
