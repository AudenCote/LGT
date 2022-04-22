import os
import sys
import itertools
#For processing .fasta files
from Bio import SeqIO


#Name of directory containing ReadyToGo files to be processed
rtg_files_handle = '../../ReadyToGo_Files'
og_files_handle = '../../allOG5Files'
#Name of output spreadsheet
output_handle = '../minor_clade_presence_photo_labelled.csv'

def clade_counts(bacbin):
	minor_clade_counts = {}
	
	already_seen = {}
	all_taxa = []
	
	for rtg_file in os.listdir(rtg_files_handle):
		if((rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b') and bacbin):
			all_taxa.append(rtg_file[:10])
		elif(not bacbin):
			all_taxa.append(rtg_file[:10])
						
	all_taxa = list(dict.fromkeys(all_taxa))
		
	for taxon in all_taxa:
		if(taxon[:8].lower() != 'am_tb_he' and taxon[:8].lower() != 'am_tb_hp'):
			try:
				minor_clade_counts[taxon[:5].lower()] += 1
			except KeyError:
				minor_clade_counts.update({taxon[:5].lower() : 1})
		elif(taxon[:8].lower() == 'am_tb_he'):
			try:
				minor_clade_counts['tb_he'] += 1
			except KeyError:
				minor_clade_counts.update({'tb_he' : 1})
		elif(taxon[:8].lower() == 'am_tb_hp'):
			try:
				minor_clade_counts['tb_hp'] += 1
			except KeyError:
				minor_clade_counts.update({'tb_hp' : 1})
		
		if(taxon[:7].lower() == 'sr_ci_s'):
			try:
				minor_clade_counts['sr_ci_s'] += 1
			except KeyError:
				minor_clade_counts.update({'sr_ci_s' : 1})
				
		
	return minor_clade_counts


def generate_dataframe(bacbin, relative, photo):

	photo_taxa = []
	if(photo):
		for line in open('../photo_taxonomy.tsv'):
			line = line.split('\t')
			if('photo' in line[2].lower()):
				photo_taxa.append(line[0].lower().strip())
				
	minor_clade_counts = clade_counts(bacbin)

	clade_codes = {}
	
	if(bacbin):
		for rtg_file in os.listdir(rtg_files_handle):
			if((rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b') or rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'):
	
				if(rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'):
					clade_code = rtg_file[:2].lower()
				elif(rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b'):
					clade_code = rtg_file[:5].lower()
					
				clade_codes.update({ clade_code : 0 })
	elif(not bacbin):
		for rtg_file in os.listdir(rtg_files_handle):
			if(rtg_file[2] == '_' and rtg_file[5] == '_'):
				
				if(photo):
					clade_codes.update({ rtg_file[:5].lower() + '_pho' : 0 })
					clade_codes.update({ rtg_file[:5].lower() + '_het' : 0 })	
				else:
					clade_codes.update({ rtg_file[:5].lower() : 0 })
	
	if(bacbin):
		taxa_per_og = dict.fromkeys(list(itertools.chain.from_iterable([record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta")) if 'OG5_' in record.id[-10:]] for rtg_file in os.listdir(rtg_files_handle) if((rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b') or rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'))))
	elif(not bacbin):
		taxa_per_og = dict.fromkeys(list(itertools.chain.from_iterable([record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta")) if 'OG5_' in record.id[-10:]] for rtg_file in os.listdir(rtg_files_handle) if('ds_store' not in rtg_file.lower()))))
	
	for gene_family in taxa_per_og:
		if(bacbin):
			codes = clade_codes.copy()
			codes.update({ 'tb_he' : 0 })
			codes.update({ 'tb_hp' : 0 })
			taxa_per_og[gene_family] = codes
		else:
			codes = clade_codes.copy()
			codes.update({ 'sr_ci_s' : 0 })
			taxa_per_og[gene_family] = codes
			
	if(bacbin):
		for rtg_file in os.listdir(rtg_files_handle):
			if((rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b') or rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'):
		
				if(rtg_file[0:2] == 'Ba' or rtg_file[0:2] == 'Za'):
					clade_code = rtg_file[:2].lower()
				elif(rtg_file[0:2] != 'Ba' and rtg_file[0:2] != 'Za' and rtg_file[4] == 'b'):
					clade_code = rtg_file[:5].lower()
					
				if(rtg_file[:8].lower() == 'am_tb_he'):
					clade_code = 'tb_he'
				elif(rtg_file[:8].lower() == 'am_tb_hp'):
					clade_code = 'tb_hp'
					
				all_ogs = [record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta"))]
				ogs_in_file = list(dict.fromkeys(all_ogs))

				for og in ogs_in_file:
					taxa_per_og[og][clade_code] += 1
	
	elif(not bacbin):
		for rtg_file in os.listdir(rtg_files_handle):
			if('ds_store' not in rtg_file.lower()):
				if(rtg_file[:7].lower() != 'sr_ci_s'):
					clade_code = rtg_file[:5].lower()
				else:
					clade_code = rtg_file[:7].lower()
				tdc = rtg_file[:10].lower()
							
				all_ogs = [record.id[-10:] for record in list(SeqIO.parse(rtg_files_handle + '/' + rtg_file, "fasta"))]
				ogs_in_file = list(dict.fromkeys(all_ogs))
				
				for og in ogs_in_file:
					if('OG5' in og):
						if(photo):
							if(tdc in photo_taxa):
								taxa_per_og[og][clade_code + '_pho'] += 1
							else:
								taxa_per_og[og][clade_code + '_het'] += 1
						else:
							taxa_per_og[og][clade_code] += 1
																
	with open(output_handle, 'w') as o:	
		o.write(',')
		for clade in taxa_per_og[next(iter(taxa_per_og))]:
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

	bacbin_only = False
	relative = False
	
	#Currently, photosynthesis differentiation only works when running NOT bacterial-bin only
	photo = False
	if(bacbin_only):
		photo = False
	
	generate_dataframe(bacbin_only, relative, photo)
		
	
main()

