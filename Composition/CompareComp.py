import os
import sys
import re
from p4 import *
from Bio import SeqIO
from tqdm import tqdm
from Bio.Blast.NCBIWWW import qblast


var.doRepairDupedTaxonNames = 1


def help():
	
	print('Invalid input')
	exit()


def get_args():
	
	if('--conserved_dir' in sys.argv or '-c' in sys.argv):
		if('--conserved_dir' in sys.argv):
			try:
				conserved_dir = sys.argv[sys.argv.index('--conserved_dir') + 1]
			except IndexError:
				help()
		elif('-c' in sys.argv):
			try:
				conserved_dir = sys.argv[sys.argv.index('-c') + 1]
			except IndexError:
				help()
	else:
		help()
		
	if('--xgt_dir' in sys.argv or '-x' in sys.argv):
		if('--xgt_dir' in sys.argv):
			try:
				xgt_dir = sys.argv[sys.argv.index('--xgt_dir') + 1]
			except IndexError:
				help()
		elif('-x' in sys.argv):
			try:
				xgt_dir = sys.argv[sys.argv.index('-x') + 1]
			except IndexError:
				help()
	else:
		help()
		
	if('--non_recip' in sys.argv):
		non_recip = True
	else:
		non_recip = False
		
	if(conserved_dir.endswith('/')):
		conserved_dir = conserved_dir[:-1]
		
	if(xgt_dir.endswith('/')):
		xgt_dir = xgt_dir[:-1]
		
	if(os.path.isdir(conserved_dir) and os.path.isdir(xgt_dir)):
		return [conserved_dir, xgt_dir, non_recip]
	else:
		help()
		
		
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
		tree.reRoot(biggestBaZa, checkBiRoot=False)
	else:
		if sizes_cladesOp:
			biggestOp = max(sizes_cladesOp, key=sizes_cladesOp.get)
			tree.reRoot(biggestOp, checkBiRoot=False)
		else:
			if sizes_cladesPl:
				biggestPl = max(sizes_cladesPl, key=sizes_cladesPl.get)
				tree.reRoot(biggestPl, checkBiRoot=False)
			else:
				if sizes_cladesAm:
					biggestAm = max(sizes_cladesAm, key=sizes_cladesAm.get)
					tree.reRoot(biggestAm, checkBiRoot=False)
				else:
					if sizes_cladesEx:
						biggestEx = max(sizes_cladesEx, key=sizes_cladesEx.get)
						tree.reRoot(biggestEx, checkBiRoot=False)
					else:
						if sizes_cladesSr:
							biggestSr = max(sizes_cladesSr, key=sizes_cladesSr.get)
							tree.reRoot(biggestSr, checkBiRoot=False)
							
	return tree
	

def get_singleton_seqs(conserved_dir, xgt_dir, non_recip):
	
	all_singleton_leaves = []
	
	if(non_recip):
		og_groups = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('og_groups.csv') }
		
	taxaselection_dict = { line.split('\t')[0].strip()[:4] + line.split('\t')[0].strip()[5:] : line.split('\t')[1].strip().split(' ')[0].split('_')[0].split('-')[0] for line in open('taxaselection.tsv') if(line.split('\t')[0].strip()[:10] in [l.strip() for l in open('xgt_taxa')])}
	
	for tree_f in tqdm(os.listdir(xgt_dir)):
		if(tree_f.endswith('.tre') or tree_f.endswith('.tree')):
			fname = xgt_dir + '/' + tree_f
			og = 'OG5_' + fname.split('OG5_')[1][:6]
			
			nexus_to_newick(fname)
			
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
			
			recips_list = []
			
			if(non_recip):
				try:
					if(og_groups[og] != 'EGT' and og_groups[og] != 'Complex' and og_groups[og] != 'Am_Ex'):
						recips_list = [og_groups[og]]
					elif(og_groups[og] == 'EGT'):
						recips_list = ['Sr', 'Pl', 'EE']
					elif(og_groups[og] == 'Am_Ex'):
						recips_list = ['Am', 'Ex']
					else:
						continue
				except KeyError:
					continue
								
			for node in tree.iterNodesNoRoot(): 	
				if node.getNChildren() == 0:						
					taxon_of_interest = node.name[:10]
					full_leaf_name = node.name
						
					if(taxon_of_interest[:2] not in recips_list and taxon_of_interest[:2] != 'Ba' and taxon_of_interest[:2] != 'Za'):												
						num_sister_taxa = 0; other_rat = 0
						while num_sister_taxa < 6:
							sister_taxa = list(dict.fromkeys([leaf[:10] for leaf in tree.getAllLeafNames(node.parent)]))
							
							num_sister_taxa = len(sister_taxa)
							
							n = 0
 							for sister_taxon in sister_taxa:
 								if(sister_taxon[:2] != taxon_of_interest[:2]):
 									n += 1
 										
 							other_rat = float(n) / float(len(sister_taxa))
 							 
							node = node.parent
							
						n = 0; maj_n = 0; mins = []
						for sis_tax in sister_taxa:
							if(sis_tax[:10] != taxon_of_interest and sis_tax[:5] == taxon_of_interest[:5]):
								if(taxaselection_dict[sis_tax[:4] + sis_tax[5:10]] == taxaselection_dict[taxon_of_interest[:4] + taxon_of_interest[5:10]]):
									n = 2
								n += 1
								
							if(sis_tax[:10] != taxon_of_interest and sis_tax[:2] == taxon_of_interest[:2]):
								maj_n += 1
								
							if(sis_tax[:2] == taxon_of_interest[:2]):
								mins.append(sis_tax[:5])
								
						if(n > 1 or (maj_n > 2 and len(dict.fromkeys(mins)) < 3) or (maj_n > 4)):
							continue
							
						if(other_rat >= .67):
							if(maj_n == 0):
								all_singleton_leaves.append((full_leaf_name, 'S'))
							else:
								all_singleton_leaves.append((full_leaf_name, 'D'))
	
	all_singleton_leaves = list(dict.fromkeys(all_singleton_leaves))
							
	with open('singletons_to_gc3.csv', 'w') as o1:
		with open('singletons_in_genus.csv', 'w') as o2:
			for leaf in all_singleton_leaves:
				if(len([lf for lf in taxaselection_dict if(taxaselection_dict[lf] == taxaselection_dict[leaf[0][:4] + leaf[0][5:10]])]) == 1):
					o1.write(leaf[0] + ',' + leaf[1] + '\n')
				else:
					o2.write(leaf[0] + ',' + leaf[1] + '\n')
			
	with open('conserved_to_gc3.csv', 'w') as o:	
		for tree_f in os.listdir(conserved_dir):
			if(tree_f.endswith('.tre') or tree_f.endswith('.tree')):
				fname = conserved_dir + '/' + tree_f
				og = 'OG5_' + fname.split('OG5_')[1][:6]
				
				var.trees = []			
															
				correct_tree(fname)													
				tree_file = fname.split('.tre')[0] + '_temp.tre'
				read(tree_file) 
				try:
					tree = var.trees[0] 	
				except:
					continue					
													
				os.remove(tree_file)
			
				for node in tree.iterNodesNoRoot(): 	
					if node.getNChildren() == 0:
						if(node.name[:10] in [leaf[0][:10] for leaf in all_singleton_leaves]):				
							o.write(node.name + '\n')
			
			
def get_GC3s():

	iters = ['singletons', 'conserved']
	
	m = open('gc3_master.csv', 'w')
	m.write('Seq,Type,GC3-Degen,ExpWrightENc,ObsWrightENc_6Fold,ObsWrightENc_No6Fold,ObsWeightedENc_6Fold,ObsWeightedENc_No6Fold\n')
	
	sd_dict = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('singletons_to_gc3.csv') }
	
	for type in iters:
		print('\nCompiling GC3 data for data type: ' + type + '\n')
		with open(type + '_seqs.fasta', 'w') as o:
			for seqid in open(type + '_to_gc3.csv'):
				seqid = seqid.split(',')[0].strip()
				
				for rtg_file in os.listdir('../../ReadyToGo_NTD'):
					if(seqid[:10] in rtg_file):
						#print(seqid.replace('LKH', '-LKH'))
						try:
							ntd = [record.seq for record in SeqIO.parse('../../ReadyToGo_NTD/' + rtg_file, 'fasta') if record.id == seqid.replace('LKH', '-LKH')][0]	
							o.write('>' + seqid + '\n' + str(ntd) + '\n\n')
						except IndexError:
							continue
							
					
		os.system('python3 CUB.py ' + type + '_seqs.fasta ' + type + '_comp universal')
		
		for l, line in enumerate(open(type + '_comp/Spreadsheets/' + type + '_comp.CompTrans.ENc.Raw.tsv')):
			line = [cell.strip() for cell in line.split('\t') if cell.strip() != '']
			if(l != 0):
				if(type != 'singletons'):
					m.write(line[0] + ',' + type + ',' + ','.join(line[6:]) + '\n')
				elif line[0] in sd_dict:
					m.write(sd_dict[line[0]] + '_' + line[0] + ',' + type + ',' + ','.join(line[6:]) + '\n')
		
		
	m.close()
		
		
def get_RSCUs():

	iters = ['singletons', 'conserved']
	
	sd_dict = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('singletons_to_gc3.csv') }
		
	if(not os.path.isdir('all_seq_fastas')):
		os.mkdir('all_seq_fastas')
	if(not os.path.isdir('RSCU_tables')):
		os.mkdir('RSCU_tables')
		
	m = open('rscu_master.csv', 'w')
	
	rscu_dict = { }
	
	for type in iters:
		print('\nCompiling RSCU data for data type: ' + type + '\n')
		for seq in SeqIO.parse(type + '_seqs.fasta', 'fasta'):
			with open('all_seq_fastas/' + seq.id + '.fasta', 'w') as o:
				o.write('>' + seq.description + '\n' + str(seq.seq))
				
			if(not os.path.isdir('RSCU_tables/' + seq.id)):
				os.system('python3 CUB.py all_seq_fastas/' + seq.id + '.fasta ' + 'RSCU_tables/' + seq.id + ' universal')
			
			rscu_dict.update({ seq.id + '_' + type : { line.split('\t')[0].strip() : line.split('\t')[2].strip() for l, line in enumerate(open('RSCU_tables/' + seq.id + '/Spreadsheets/' + seq.id + '.RSCU.tsv')) if l != 0 } })
			
	codons_list = []
	for seq in rscu_dict:
		for codon in rscu_dict[seq]:
			if(codon not in codons_list):
				codons_list.append(codon)
			
	m.write('Seq,Type,' + ','.join(codons_list) + '\n')
	for seq in rscu_dict:
		if(sum([float(val) for val in rscu_dict[seq].values()]) > 0):
			if('_'.join(seq.split('_')[:-1]) not in sd_dict):
				m.write('_'.join(seq.split('_')[:-1]) + ',' + seq.split('_')[-1])
			else:
				m.write(sd_dict['_'.join(seq.split('_')[:-1])] + '_' + '_'.join(seq.split('_')[:-1]) + ',' + seq.split('_')[-1])
			for codon in codons_list:
				if(codon in rscu_dict[seq]):
					m.write(',' + str(rscu_dict[seq][codon]))
				else:
					m.write(',0')
			m.write('\n')		
			
	m.close()
	
	
def blast_transcriptomes():

	if(not os.path.isdir('BLASTn_records')):
		os.mkdir('BLASTn_records')
		
	if(not os.path.isdir('BLASTx_records')):
		os.mkdir('BLASTx_records')
# 
# 	all_seqs = list(SeqIO.parse('singletons_seqs.fasta', 'fasta'))
# 	print('\nBLAST-ing Singleton Sequences...\n')
# 	for record in tqdm(all_seqs):
# 		if(not os.path.isfile('BLASTn_records/' + record.id)):
#  			with open('BLASTn_records/' + record.id, 'w') as o:
#  				results = qblast(program = 'blastn', sequence = record.seq, database = 'nr', hitlist_size = 10)	
#  				o.write(results.read())
# 				results.close()
# 		if(not os.path.isfile('BLASTx_records/' + record.id)):
#  			with open('BLASTx_records/' + record.id, 'w') as o:
#  				results = qblast(program = 'blastx', sequence = record.seq, database = 'nr', hitlist_size = 10)	
#  				o.write(results.read())
# 				results.close()
	
	blast_types = ['BLASTn', 'BLASTx']
	
	for blast_type in blast_types:			
		blast_hits = { }
	
		for blast_file in os.listdir(blast_type + '_records'):
			lines = [line for line in open(blast_type + '_records/' + blast_file).readlines()]
			qlen = 0
			blast_hits.update({ blast_file : { } })
			for i, l1 in enumerate(lines):
				if('Iteration_query-len' in l1):
					qlen = float(l1.split('>')[1].split('<')[0])
				elif('Hit_def' in l1):
					print(blast_file, l1)
					org = ' '.join(l1.split('>')[1].split('[')[-1].split(']')[0].split(' ')[:2])
					id = len = query_from = query_to = ''
					for l2 in lines[i:]:
						if('Hsp_identity' in l2):
							id = float(l2.split('>')[1].split('<')[0])
						elif('Hsp_query-to' in l2):
							query_to = float(l2.split('>')[1].split('<')[0])
						elif('Hsp_query-from' in l2):
							query_from = float(l2.split('>')[1].split('<')[0])
						elif('Hsp_align-len' in l2):
							len = float(l2.split('>')[1].split('<')[0])
							break
							
							
					blast_hits[blast_file].update({ org : { 'id' : int(round(id / len, 2) * 100), 'cov' : int(round(float(abs(query_to - query_from)) / float(qlen), 2)* 100) } })
				
		with open(blast_type + '_hits_master_transcriptomes.csv', 'w') as o:
			o.write('Query,Hit,ID,Cov\n')
			for sing_seq in blast_hits:
				for hit in blast_hits[sing_seq]:
					o.write(sing_seq + ',' + hit + ',' + str(blast_hits[sing_seq][hit]['id']) + ',' + str(blast_hits[sing_seq][hit]['cov']) + '\n')
				

def blast_genomes():

	if(not os.path.isdir('BLASTp_records')):
		os.mkdir('BLASTp_records')
	
	genomic_taxa = [line.split('\t')[0].strip()[:10] for line in open('genomic_taxa.tsv') if line.split('\t')[0].strip() != '']
	
	genomic_singletons = [line.split(',')[0].strip() for line in open('singletons_to_gc3.csv') if line.split(',')[0].strip()[:10] in genomic_taxa]
				
	# for singleton in tqdm(genomic_singletons):
# 		#OrthoMCL (allOG5Files)
# 		if(singleton[6:10] == singleton[6:10].lower()):
# 			try:
# 				og_file = '../../allOG5Files/' + [og_file for og_file in os.listdir('../../allOG5Files') if(singleton.split('OG5_')[-1][:6] in og_file)][0]
# 				aa = [rec.seq for rec in SeqIO.parse(og_file, 'fasta') if(rec.description == singleton)][0]
# 			except IndexError:
# 				continue
# 		#ReadyToGo_AA
# 		else:
# 			try:
# 				rtg_file = '../../ReadyToGo_AA/' + [rtg_file for rtg_file in os.listdir('../../ReadyToGo_AA') if(singleton[:10] in rtg_file)][0]
# 				aa = [rec.seq for rec in SeqIO.parse(rtg_file, 'fasta') if(rec.description == singleton)][0]
# 			except IndexError:
# 				continue
# 		
# 		if(not os.path.isfile('BLASTp_records/' + singleton)):
#  			with open('BLASTp_records/' + singleton, 'w') as o:
#  				results = qblast(program = 'blastp', sequence = aa, database = 'nr', hitlist_size = 10)	
#  				o.write(results.read())
# 				results.close()
				
	blast_hits = { }
	
	for blast_file in os.listdir('BLASTp_records'):
		lines = [line for line in open('BLASTp_records/' + blast_file).readlines()]
		qlen = 0
		blast_hits.update({ blast_file : { } })
		for i, l1 in enumerate(lines):
			if('Iteration_query-len' in l1):
				qlen = float(l1.split('>')[1].split('<')[0])
			elif('Hit_def' in l1):
				print(blast_file, l1)
				org = ' '.join(l1.split('>')[1].split('[')[-1].split(']')[0].split(' ')[:2])
				id = len = query_from = query_to = ''
				for l2 in lines[i:]:
					if('Hsp_identity' in l2):
						id = float(l2.split('>')[1].split('<')[0])
					elif('Hsp_query-to' in l2):
						query_to = float(l2.split('>')[1].split('<')[0])
					elif('Hsp_query-from' in l2):
						query_from = float(l2.split('>')[1].split('<')[0])
					elif('Hsp_align-len' in l2):
						len = float(l2.split('>')[1].split('<')[0])
						break
							
							
				blast_hits[blast_file].update({ org : { 'id' : int(round(id / len, 2) * 100), 'cov' : int(round(float(abs(query_to - query_from)) / float(qlen), 2)* 100) } })
				
		with open('BLASTp_hits_master_genomes.csv', 'w') as o:
			o.write('Query,Hit,ID,Cov\n')
			for sing_seq in blast_hits:
				for hit in blast_hits[sing_seq]:
					o.write(sing_seq + ',' + hit.replace(',', ' ') + ',' + str(blast_hits[sing_seq][hit]['id']) + ',' + str(blast_hits[sing_seq][hit]['cov']) + '\n')
				
				
def make_plots():

	taxa = list(dict.fromkeys([line.strip()[:10] for line in open('singletons_to_gc3.csv')]))
	
	for taxon in taxa:
		os.system('/Library/Frameworks/R.framework/Versions/4.0/Resources/Rscript PlotComps_Noise.r ' + taxon)
		
		
def summarise_transcriptomes():

	genomic_taxa = [line.split('\t')[0].strip()[:10] for line in open('genomic_taxa.tsv') if line.split('\t')[0].strip() != '']
	
	o = open('TranscriptomeSingletonSummary.csv', 'w')
	o.write('Singleton Seq,Single-/Doubleton,Alone in Genus,Nucleotide Sequence Length,Coverage,Compositional Outlier,Sisters,BLASTn Hits (ID Cov),BLASTx Hits (ID Cov)\n')
	
	singletons_to_gc3 = [line.split(',')[0].strip() for line in open('singletons_to_gc3.csv')]
	sd_dict = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('singletons_to_gc3.csv') }
	genus_dict = { line.split(',')[0].strip() : line.split(',')[2].strip() for line in open('singletons_to_gc3.csv') }
	
	print('\nSummarizing transcriptome singleton data...')
	for singleton in tqdm(singletons_to_gc3):
		if(singleton[:10] not in genomic_taxa):
		
			if(sd_dict[singleton] == 'S'):
				singledouble = 'Singleton'
			elif(sd_dict[singleton] == 'D'):
				singledouble = 'Doubleton'
				
			if(genus_dict[singleton] == 'alone'):
				alone = 'Yes'
			elif(genus_dict[singleton] == 'friends'):
				alone = ''
			
			try:
				length = str(len(str([rec.seq for rec in SeqIO.parse('singletons_seqs.fasta', 'fasta') if rec.description == singleton][0])))
			except:
				length = 'SNF'
			
			try:
				coverage = str(int(singleton.split('Cov')[-1].split('_')[0]))
			except:
				coverage = 'NA'
				
			if(singleton in [line.strip() for line in open('outliers.csv')]):
				outlier = 'Yes'
			else:
				outlier = ''
				
			sisters = ' '.join(list(dict.fromkeys([sister.strip()[:10] for sister in line for line in open('Contamination_report/XGT/report_walk_contamination_single-all.txt') for sister in line.split('\t')[4].split(',') if(singleton in line.split('\t')[2] and sister.strip() != '')])))
			
			blastn_hits = ' '.join(list(dict.fromkeys(sorted([line.split(',')[1].strip() + ' (' + str(line.split(',')[2].strip()) + ' ' + str(line.split(',')[3].strip()) + ') ' for line in open('BLASTn_hits_master_transcriptomes.csv') if singleton in line], key = lambda x : int(x.split('(')[-1].split(' ')[0].strip()), reverse = True))))
			
			blastx_hits = ' '.join(list(dict.fromkeys(sorted([line.split(',')[1].strip() + ' (' + str(line.split(',')[2].strip()) + ' ' + str(line.split(',')[3].strip()) + ') ' for line in open('BLASTx_hits_master_transcriptomes.csv') if singleton in line], key = lambda x : int(x.split('(')[-1].split(' ')[0].strip()), reverse = True))))
			
			o.write(singleton + ',' + singledouble + ',' + alone + ',' + length + ',' + coverage + ',' + outlier + ',' + sisters + ',' + blastn_hits + ',' + blastx_hits + '\n')
			
			#if(s >= 2):
			#s	exit()
		
	o.close()
	
	
def summarise_genomes():

	if(not os.path.isdir('Scratch')):
		os.mkdir('Scratch')

	genomic_taxa = [line.split('\t')[0].strip()[:10] for line in open('genomic_taxa.tsv') if line.split('\t')[0].strip() != '']
	
	o = open('GenomeSingletonSummary.csv', 'w')
	o.write('Singleton Seq,Data Source,Accession/Link,Annotated,Single-/Doubleton,Alone in Genus,Nucleotide Sequence Length,Sisters,BLASTp Hits (ID Cov)\n')
	
	genomic_singletons = [line.split(',')[0].strip() for line in open('singletons_to_gc3.csv') if(line.split(',')[0].strip()[:10] in genomic_taxa and 'OG5_' in line.split(',')[0].strip())]
	sd_dict = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('singletons_to_gc3.csv') }
	genus_dict = { line.split(',')[0].strip() : line.split(',')[2].strip() for line in open('singletons_to_gc3.csv') }
	
	print('\nSummarising data for all genomic singletons...')
	for singleton in tqdm(genomic_singletons):
	
		if(singleton[6:10].lower() == singleton[6:10]):
			source = 'OrthoMCL'
		else:
			source = 'GenBank'
		
		if(sd_dict[singleton] == 'S'):
			singledouble = 'Singleton'
		elif(sd_dict[singleton] == 'D'):
			singledouble = 'Doubleton'
						
		if(genus_dict[singleton] == 'alone'):
			alone = 'Yes'
		elif(genus_dict[singleton] == 'friends'):
			alone = 'No'
		
		flag = False
		for rtg_file in os.listdir('../../ReadyToGo_AA'):
			if(singleton[:10] in rtg_file):
				flag = True
				if('GCA_' in rtg_file):
					accession = 'GCA_' + rtg_file.split('GCA_')[-1].split('_')[0].replace('_', '')
				elif('GCA' in rtg_file and 'GCA_' not in rtg_file):
					accession = rtg_file[rtg_file.index('GCA'):].split('_')[0]
				elif('GCF_' in rtg_file):
					accession = 'GCF_' + rtg_file.split('GCF_')[-1].split('_')[0].replace('_', '')
				elif('GCF' in rtg_file and 'GCF_' not in rtg_file):
					accession = rtg_file[rtg_file.index('GCA'):].split('_')[0]
				else:
					accession = ''
				
				try:
					length = str(len([rec.seq for rec in SeqIO.parse('../../ReadyToGo_AA/' + rtg_file, 'fasta') if rec.description == singleton][0]) * 3)
				except IndexError:
					length = 'SNF'
					
		if(not flag):
			try:
				length = str(len([rec.seq for rec in SeqIO.parse('../../allOG5Files/OG5_' + singleton.split('OG5_')[-1][:6], 'fasta') if rec.description == singleton][0]) * 3)	
			except IndexError:
				length = 'SNF'
			
			accession = [line.split('\t')[4] for line in open('data_sources_OrthoMCL-5.txt') if line.split('\t')[0] == singleton[6:10]][0].strip()
					
		if('GCA' not in accession and 'GCF' not in accession):
			annotated = 'No'
		else:
			os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + accession[:3] + '/' + accession[4:7] + '/' + accession[7:10] + '/' + accession[10:13] + ' Scratch > /dev/null')
			
			annotated = 'No'
			for dir in os.listdir('Scratch'):
				if(os.path.isdir('Scratch/' + dir)):
					for subdir in os.listdir('Scratch/' + dir):
						if(os.path.isdir('Scratch/' + dir + '/' + subdir)):
							for file in os.listdir('Scratch/' + dir + '/' + subdir):
								if('cds_from_genomic' in file):
									annotated = 'Yes'
				
			os.system('rm -rf Scratch')
			os.mkdir('Scratch')
						
		sisters = ' '.join(list(dict.fromkeys([sister.strip()[:10] for sister in line for line in open('Contamination_report/XGT/report_walk_contamination_single-all.txt') for sister in line.split('\t')[4].split(',') if singleton in line.split('\t')[2] if sister.strip() != ''])))
		
		blastp_hits = ' '.join(list(dict.fromkeys(sorted([line.split(',')[1].strip() + ' (' + str(line.split(',')[2].strip()) + ' ' + str(line.split(',')[3].strip()) + ') ' for line in open('BLASTp_hits_master_genomes.csv') if singleton in line], key = lambda x : int(x.split('(')[-1].split(' ')[0].strip()), reverse = True))))
		
		if(source == 'OrthoMCL'):
			annotated = 'Yes'
		
		o.write(singleton + ',' + source + ',' + accession + ',' + annotated + ',' + singledouble + ',' + alone + ',' + length + ',' + sisters + ',' + blastp_hits + '\n')	
		
		#if(i >= 2):
			#exit()	
		
	o.close()


def main():
	
	conserved_dir, xgt_dir, non_recip = get_args()
	
	#get_singleton_seqs(conserved_dir, xgt_dir, non_recip)
	
	#get_GC3s()
	
	#get_RSCUs()
	
	#blast_transcriptomes()
	
	#blast_genomes()
	
	#make_plots()
	
	#summarise_genomes()
	
	summarise_transcriptomes()
	
	
main()
	