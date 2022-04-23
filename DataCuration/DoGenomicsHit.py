import os
import sys

from subprocess import check_output
import gzip
from shutil import copyfileobj

from time import sleep
import itertools
from tqdm import tqdm

from Bio import SeqIO
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML
import math

from VEuPath_REST import VEuPath_REST
import requests
import selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options

chrome_options = Options()
chrome_options.add_argument("--headless")


def help():

	print('Aaaaah!')
	exit()

def get_accession_numbers():

	print('\nRetrieving Accession Numbers...\n\n######################\n')
	
	acc_by_tax = { }
	
	#Getting accession numbers from NCBI using Entrez search
	c = 0
	for l, line in enumerate(open('genomic_taxa.tsv')):
		if(l != 0 and ' '.join(line.split('\t')[1].strip().split(' ')[:2]) not in acc_by_tax):
			c += 1
			print(str(c) + '. ' + ' '.join(line.split('\t')[1].strip().split(' ')[:2]))
			output = str(check_output('esearch -db assembly -query "' + ' '.join(line.split('\t')[1].strip().split(' ')[:2]) + '" | efetch -format docsum | tee >(grep "Accession") > /dev/null', shell = True, executable='/bin/bash'))
			try:
				acc_by_tax.update({ ' '.join(line.split('\t')[1].strip().split(' ')[:2]) : output[output.index('GCF'):output.index('GCF') + 15] })
			except:
				acc_by_tax.update({ ' '.join(line.split('\t')[1].strip().split(' ')[:2]) : 'NA' })
				
	print('\nDone retrieving Accession Numbers.\n\n######################\n')
	
	print(acc_by_tax)
	
	return acc_by_tax
			
			
def query_NCBI(seqid):

	#acc_by_tax = { line.split(',')[0].strip() : line.split(',')[1].strip() for line in open('../Spreadsheets/AccByTax.csv') if(line.split(',')[0].strip() != '' and line.split(',')[1].strip() != '') }
	
	acc_by_tax = {'Entamoeba dispar': 'GCF_000209125.1', 'Entamoeba histolytica': 'GCF_000208925.1', 'Entamoeba invadens': 'GCF_000330505.1', 'Acanthamoeba castellanii': 'GCF_000313135.1', 'Acytostelium subglobosum': 'GCF_000787575.1', 'Dictyostelium discoideum': 'GCF_000004695.1', 'Thecamonas': 'GCF_000142905.1', 'Guillardia theta': 'GCF_000002975.1', 'Bodo saltans': 'NA', 'Leishmania braziliensis': 'GCF_000002845.2', 'Leishmania infantum': 'GCF_000002875.2', 'Leishmania major': 'GCF_000002725.2', 'Leishmania mexicana': 'GCF_000234665.1', 'Leptomonas pyrrhocoris': 'GCF_001293395.1', 'Phytomonas sp': 'NA', 'Strigomonas culicis': 'NA', 'Trypanosoma brucei': 'GCF_000002445.1', 'Trypanosoma congolense': 'NA', 'Trypanosoma cruzi': 'GCF_000209065.1', 'Trypanosoma vivax': 'NA', 'Giardia instestinalis': 'GCF_000002435.2', 'Giardia lamblia': 'GCF_000002435.2', 'Spironucleus salmonicida': 'NA', 'Naegleria gruberi': 'GCF_000004985.1', 'Tritrichomonas foetus': 'NA', 'Trichomonas vaginalis': 'GCF_000002825.2', 'Monosiga brevicollis': 'GCF_000002865.3', 'Aspergillus fumigatus': 'GCF_000002655.1', 'Allomyces macrogynus': 'NA', 'Emericella nidulans': 'GCF_000149205.1', 'Aspergillus oryzae': 'GCF_000184455.2', 'Batrachochytrium dendrobatidis': 'GCF_000203795.1', 'Candida albicans': 'GCF_000182965.3', 'Conidiobolus coronatus': 'NA', 'Candida glabrata': 'GCF_000002545.3', 'Coccidioides immitis': 'GCF_000149335.2', 'Cryptococcus neoformans': 'GCF_000091045.1', 'Coccidioides posadasii': 'GCF_000151335.2', 'Debaryomyces hansenii': 'GCF_000006445.2', 'Dacryopinax primogenitus': 'NA', 'Enterocytozoon bieneusi': 'GCF_000209485.1', 'Encephalitozoon cuniculi': 'GCF_000091225.1', 'Eremothecium gossypii': 'GCF_000091025.4', 'Encephalitozoon intestinalis': 'GCF_000146465.1', 'Gonapodya prolifera': 'NA', 'Gibberella zeae': 'GCF_000240135.2', 'Kluyveromyces lactis': 'GCF_000002515.2', 'Laccaria bicolor': 'GCF_000143565.1', 'Lichtheimia corymbifera': 'NA', 'Mitosporidium daphniae': 'GCF_000760515.2', 'Magnaporthe grisea': 'GCF_004355905.1', 'Melampsora laricipopulina': 'GCF_000204055.1', 'Nosema Apis': 'NA', 'Nosema ceranae': 'GCF_000988165.1', 'Neurospora crassa': 'GCF_000182925.1', 'Nematocida parisii': 'GCF_000250985.1', 'Phanerochaete chrysosporium': 'NA', 'Pseudoloma neurophilia': 'NA', 'Pichia stipitis': 'GCF_000209165.1', 'Rozella allomycis': 'NA', 'Rhizopus delemar': 'NA', 'Rhizophagus irregularis': 'GCF_000439145.1', 'Rhizopus microsporus': 'GCF_002708625.1', 'Saccharomyces cerevisiae': 'GCF_000146045.2', 'Saitoella complicata': 'GCF_001661265.1', 'Schizosaccharomyces pombe': 'GCF_000002945.1', 'Spizellomyces punctatus': 'GCF_000182565.1', 'Trichophyton rubrum': 'GCF_000151425.1', 'Wallemia sebi': 'NA', 'Yarrowia lipolytica': 'GCF_000002525.2', 'Capsaspora owczarzaki': 'GCF_000151315.2', 'Sphaeroforma arctica': 'GCF_001186125.1', 'Aedes aegypti': 'GCF_002204515.2', 'Anopheles gambiae': 'GCF_000005575.2', 'Apis mellifera': 'GCF_003254395.2', 'Acyrthosiphon pisum': 'GCF_005508785.1', 'Brugia malayi': 'GCF_000002995.2', 'Bombyx mori': 'GCF_014905235.1', 'Caenorhabditis briggsae': 'GCF_000004555.1', 'Caenorhabditis elegans': 'GCF_000002985.6', 'Crassostrea gigas': 'GCF_902806645.1', 'Ciona intestinalis': 'GCF_000224145.3', 'Canis lupus': 'GCF_000002285.5', 'Culex pipiens': 'GCF_016801865.1', 'Capitella teleta': 'NA', 'Drosophila melanogaster': 'GCF_000001215.2', 'Daphnia pulex': 'NA', 'Danio rerio': 'GCF_000002035.6', 'Equus caballus': 'GCF_002863925.1', 'Gallus gallus': 'GCF_016699485.2', 'Helobdella robusta': 'GCF_000326865.1', 'Homo sapiens': 'GCF_000001405.3', 'Ixodes scapularis': 'GCF_016920785.1', 'Monodelphis domestica': 'GCF_000002295.2', 'Macaca mulatta': 'GCF_003339765.1', 'Mus musculus': 'GCF_000001635.2', 'Nematostella vectensis': 'GCF_000209225.1', 'Ornithorhynchus anatinus': 'GCF_004115215.2', 'Oikopleura dioica': 'NA', 'Pediculus humanus': 'GCF_000006295.1', 'Pan troglodytes': 'GCF_002880755.1', 'Rattus norvegicus': 'GCF_015227675.2', 'Saccoglossus kowalevskii': 'GCF_000003605.2', 'Schistosoma mansoni': 'GCF_000237925.1', 'Trichoplax adhaerens': 'GCF_000150275.1', 'Thelohanellus kitauei': 'NA', 'Tetraodon nigroviridis': 'NA', 'Takifugu rubripes': 'GCF_901000725.2', 'Trichinella spiralis': 'GCF_000181795.1', 'Arabidopsis thaliana': 'GCF_000001735.4', 'Amborella trichopoda': 'GCF_000471905.2', 'Bathycoccus prasinos': 'GCF_002220235.1', 'Chlamydomonas reinhardtii': 'GCF_000002595.1', 'Elaeis guineensis': 'GCF_000442705.1', 'Helianthus annus': 'GCF_002127325.2', 'Helicosporidium sp': 'NA', 'Klebsormidium nitens': 'NA', 'Micromonas sp': 'NA', 'Marchantia polymorpha': 'NA', 'Nicotiana tabacum': 'GCF_000715135.1', 'Ostreococcus lucimarinus': 'GCF_000092065.1', 'Oryza sativa': 'GCF_001433935.1', 'Ostreococcus tauri': 'GCF_000214015.2', 'Physcomitrella patens': 'GCF_000002425.4', 'Pinus taeda': 'NA', 'Ricinus communis': 'GCF_000151685.1', 'Volvox carteri': 'GCF_000143455.1', 'Cyanidioschyzon merolae': 'GCF_000091205.1', 'Galdieria sulphuraria': 'GCF_000341285.1', 'Babesia bovis': 'GCF_000165395.1', 'Cryptosporidium hominis': 'GCF_000006425.1', 'Cryptosporidium muris': 'GCF_000006515.1', 'Eimeria acervulina': 'GCF_000499425.1', 'Eimeria maxima': 'GCF_000499605.1', 'Eimeria tenella': 'GCF_000499545.2', 'Neospora caninum': 'GCF_000208865.1', 'Plasmodium berghei': 'GCF_900002375.2', 'Plasmodium chabaudi': 'GCF_900002335.2', 'Plasmodium falciparum': 'GCF_000002765.5', 'Plasmodium knowlesi': 'GCF_000006355.2', 'Plasmodium vivax': 'GCF_000002415.2', 'Plasmodium yoelii': 'GCF_900002385.2', 'Theileria annulata': 'GCF_000003225.4', 'Toxoplasma gondii': 'GCF_000006565.2', 'Theileria parva': 'GCF_000165365.1', 'Vitrella brassicaformis': 'NA', 'Ichthyophthirius multifiliis': 'GCF_000220395.1', 'Oxytricha trifallax': 'NA', 'Pseudocohnilembus persalinus': 'NA', 'Paramecium tetraurelia': 'GCF_000141845.1', 'Stentor coeruleus': 'NA', 'Stylonychia lemnae': 'NA', 'Tetrahymena thermophila': 'GCF_000189635.1', 'Symbiodinium microadriaticum': 'NA', 'Perkinsus marinus': 'GCF_000006405.1', 'Bigelowiella natans': 'GCF_000002455.1', 'Paulinella chromatophora': 'NA', 'Reticulomyxa filosa': 'NA', 'Aureococcus anophagefferens': 'GCF_000186865.1', 'Aphanomyces astaci': 'GCF_000520075.1', 'Aphanomyces invadans': 'GCF_000520115.1', 'Blastocystis hominis': 'GCF_000151665.1', 'Ectocarpus siliculosus': 'NA', 'Nannochloropsis gaditana': 'GCF_000240725.1', 'Phytophthora infestans': 'GCF_000142945.1', 'Phytophthora ramorum': 'NA', 'Phaeodactylum tricornutum': 'GCF_000150955.2', 'Saprolegnia diclina': 'GCF_000281045.1', 'Saprolegnia parasitica': 'GCF_000151545.1', 'Thalassiosira pseudonana': 'GCF_000149405.2'}
	
	already_reffed_by_tax = { }
	for gendir in os.listdir('GenomicCDS'):
		if(gendir[:10] == seqid[:10]):
			if(os.path.isfile('GenomicCDS/' + gendir + '/BLAST/annotated_cds.fasta')):
				already_reffed_by_tax.update({ gendir[:10] : 'GenomicCDS/' + gendir + '/BLAST/annotated_cds.fasta' })
				break
				
	org_name = [' '.join(line.split('\t')[1].strip().split(' ')[:2]) for line in open('genomic_taxa.tsv') if seqid[:10] in line][0]
	
	if(acc_by_tax[org_name] != 'NA'):
		if(seqid[:10] not in already_reffed_by_tax):
			print('\nAttempting to retreive annotated CDS for taxon ' + seqid[:10] + ' using NCBI accession ' + acc_by_tax[org_name] + '\n')
			os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + acc_by_tax[org_name][:3] + '/' + acc_by_tax[org_name][4:7] + '/' + acc_by_tax[org_name][7:10] + '/' + acc_by_tax[org_name][10:13] + ' GenomicCDS/' + seqid + ' > /dev/null')
		
			for num_dir in os.listdir('GenomicCDS/' + seqid):
				if(os.path.isdir('GenomicCDS/' + seqid + '/' + num_dir)):
 					#Check if multiple versions of the reference database are returned
					n_versions = 0; acc_dirs = []
					for acc_dir in os.listdir('GenomicCDS/' + seqid + '/' + num_dir):
						if(os.path.isdir('GenomicCDS/' + seqid + '/' + num_dir + '/' + acc_dir)):
							n_versions += 1; acc_dirs.append('GenomicCDS/' + seqid + '/' + num_dir + '/' + acc_dir)
							
					#Get latest version
					if(n_versions > 1):
						max = 0; best = ''
						for adir in acc_dirs:
							if(int(adir.split('/')[-1].split('.')[1][0]) > max):
								max = int(adir.split('/')[-1].split('.')[1][0]); best = adir
						acc_dirs = [best]
											
					for file in os.listdir(acc_dirs[0]):
						if('cds_from_genomic' in file):
							os.system('cp ' + acc_dirs[0] + '/' + file + ' GenomicCDS/' + seqid + '/annotated_cds.fasta.gz')
							break
								
			for subdir in os.listdir('GenomicCDS/' + seqid):
				if(os.path.isdir('GenomicCDS/' + seqid + '/' + subdir)):
					os.system('rm -rf GenomicCDS/' + seqid + '/' + subdir)
					
			if(not os.path.isdir('GenomicCDS/' + seqid + '/BLAST')):
				os.mkdir('GenomicCDS/' + seqid + '/BLAST')
								
			if(os.path.isfile('GenomicCDS/' + seqid + '/annotated_cds.fasta.gz')):
				with gzip.open('GenomicCDS/' + seqid + '/annotated_cds.fasta.gz', 'rb') as f_in:
					with open('GenomicCDS/' + seqid + '/BLAST/annotated_cds.fasta', 'wb') as f_out:
						copyfileobj(f_in, f_out)
						os.system('rm -f ' + str(f_in.name))
		else:
			if(not os.path.isdir('GenomicCDS/' + seqid + '/BLAST')):
				os.mkdir('GenomicCDS/' + seqid + '/BLAST')
				
			os.system('cp ' + already_reffed_by_tax[seqid[:10]] + 'GenomicCDS/' + seqid + '/BLAST/annotated_cds.fasta')
	
						
def query_VEuPath(seqid, prot_seq):

	org_name = [' '.join(line.split('\t')[1].strip().split(' ')[:2]) for line in open('genomic_taxa.tsv') if seqid[:10] in line][0]
			
	print('\nBLAST-ing putative XGT sequence in the VEuPATH database for taxon ' + org_name + '\n')
	
	org_name = [name for name in VEuPath_REST.avail_orgs if(org_name.lower() in name.lower())][0]
	
	r = requests.post(url = VEuPath_REST.BLAST_URL, json = VEuPath_REST.get_post_rqst(str(prot_seq.seq), org_name))
	to_write = r.text

	with open('GenomicCDS/' + seqid + '/veupath_blast_out.txt', 'w') as o:
		o.write(to_write)
			
	blast_seqs = []
	for line in open('GenomicCDS/' + seqid + '/veupath_blast_out.txt'):
		line = line.split('\t')
		if(org_name in line[1]):
			blast_seqs.append(line[0])
		
	taxonomy = [line.split(',')[1].strip() for line in open('taxaselection.csv') if seqid[:10] in line]
	
	if(len(blast_seqs) == 0 or len(taxonomy) == 0):
		return None
		
	taxonomy = taxonomy[0]
		
	if('amoeb' in taxonomy.lower()):
		db_url = 'https://amoebadb.org/amoeba/app/record/genomic-sequence/'
	elif('trypanosom' in taxonomy.lower()):
		db_url = 'https://tritrypdb.org/tritrypdb/app/record/genomic-sequence/'
	elif('giard' in taxonomy.lower()):
		db_url = 'https://giardiadb.org/giardiadb/app/record/genomic-sequence/'
	elif('trichomon' in taxonomy.lower()):
		db_url = 'https://trichdb.org/trichdb/app/record/genomic-sequence/'
	elif('cryptospor' in taxonomy.lower()):
		db_url = 'https://cryptodb.org/cryptodb/app/record/genomic-sequence/'
	elif('fung' in taxonomy.lower()):
		db_url = 'https://fungidb.org/fungidb/app/record/genomic-sequence/'
	elif('microsporid' in taxonomy.lower()):
		db_url = 'https://microsporidiadb.org/micro/app/record/genomic-sequence/'
	elif('piroplasm' in taxonomy.lower()):
		db_url = 'https://piroplasmadb.org/piro/app/record/genomic-sequence/'
	elif('plasmod' in taxonomy.lower()):
		db_url = 'https://plasmodb.org/plasmo/app/record/genomic-sequence/'
	elif('toxoplasm' in taxonomy.lower()):
		db_url = 'https://toxodb.org/toxo/app/record/genomic-sequence/'
	else:
		db_url = 'https://vectorbase.org/vectorbase/app/record/genomic-sequence/'
				
	downloaded_genomes = []
	for blast_seq in blast_seqs:
			
		driver = webdriver.Chrome(executable_path='/Users/katzlab/Desktop/programs/chromedriver', chrome_options=chrome_options)
		driver.implicitly_wait(30)
			
		contig_url = db_url + blast_seq
		driver.get(contig_url)
		try:
			link = driver.find_element(By.XPATH, '//*[text()="Data files"]')
		except selenium.common.exceptions.NoSuchElementException:
			continue
			
		files_url = link.get_attribute('href') + 'fasta/data/'	
		
		driver.get(files_url)
		a_list = driver.find_elements_by_tag_name('iframe')
		src = ''
		for link in a_list:
			if('data' in link.get_attribute('src').lower()):
				src = link.get_attribute('src')
		
		driver.get(src)
		a_list = driver.find_elements_by_tag_name('a')
		file_name = ''
		for link in a_list:
			if('annotated' in link.get_attribute('href').lower() and 'cds' in link.get_attribute('href').lower()):
				file_name = link.get_attribute('href')
					
		driver.quit()
						
		if(file_name in downloaded_genomes):
			continue
		else:
			downloaded_genomes.append(file_name)
		try:
			r = requests.get(file_name, allow_redirects=True)
		except requests.exceptions.MissingSchema:
			continue
				
		os.system('mkdir GenomicCDS/' + seqid + '/' + blast_seq)
		
		open('GenomicCDS/' + seqid + '/' + blast_seq + '/annotated_cds.fasta', 'wb').write(r.content)	
			
			
def make_DBs():

	all_genomic_taxa = { line.split('\t')[0].strip() : ' '.join(line.split('\t')[1].strip().split(' ')[:2]) for line in open('genomic_taxa.tsv') if(line.split('\t')[0].strip() != '' and line.split('\t')[1].strip() != '') }

	problem_ogs = { }
	
	sings = open('genomic_singleton_seqs.csv', 'w')
	
	print('\nMaking databases...\n')
	for l, line in tqdm(enumerate(open('genomic_singletons_to_BLAST.csv'))):
		if(l != 0):
			seq_id = [cell.strip() for cell in line.split(',') if cell.strip() != ''][0]
			if(seq_id[:10] in all_genomic_taxa and ('_len' in seq_id.lower() or ('contig' not in seq_id.lower() and '_len' not in seq_id.lower()))):
				og = 'OG5_' + seq_id.split('OG5_')[-1][:6]; taxon = seq_id[:10]
				print('\nProcessing ' + og + ' in taxon ' + taxon + '\n')
				
				if(not os.path.isdir('GenomicCDS/' + seq_id.replace('PAm', '-PAm'))):
					os.mkdir('GenomicCDS/' + seq_id.replace('PAm', '-PAm'))
					
				veupath = False
				for spec in VEuPath_REST.avail_orgs:
					if(all_genomic_taxa[taxon].lower() in spec.lower()):
						veupath = True
						
				prot_seq = ''
				if(seq_id[6:10] == 'aaeg'):
					seq_id = seq_id.replace('PAm', '-PAm')
				print(seq_id)
				if(seq_id[6:10].lower() == seq_id[6:10]):
					for og_file in os.listdir('../../allOG5Files'):
						if(seq_id.split('OG5_')[-1][:6] in og_file):
							prot = [record for record in SeqIO.parse('../../allOG5Files/' + og_file, 'fasta') if record.id == seq_id.replace('LKH', '-LKH')][0]
							prot_seq = prot
				else:
					for rtg_file in os.listdir('../../ReadyToGo_AA'):
						if(seq_id[:10] in rtg_file):
							try:
								prot = [record for record in SeqIO.parse('../../ReadyToGo_AA/' + rtg_file, 'fasta') if record.id == seq_id.replace('LKH', '-LKH')][0]	
								prot_seq = prot
							except:
								continue
				
				with open('GenomicCDS/' + seq_id.replace('-', '') + '/singleton_seq.fasta', 'w') as o:
					try:
						o.write('>' + prot_seq.description + '\n' + str(prot_seq.seq))
					except AttributeError:
						print('\nNo annotated CDS were retrieved for OG ' + og + '. Skipping this one for now, it may be worth trying to download this genome manually if possible\n')
						sings.write(seq_id + ',No Singleton Sequence Found\n')
						os.system('rm -r GenomicCDS/' + seq_id)
						continue
																			
				#Getting the annotated CDS files from either NCBI or VEuPathDB
				annot_flag = 0
				for dir in os.listdir('GenomicCDS/' + seq_id):
					if(os.path.isfile('GenomicCDS/' + seq_id + '/' + dir + '/annotated_cds.fasta')):
						annot_flag = 1
				if(annot_flag == 0):
					if(veupath):
						query_VEuPath(seq_id, prot_seq)
					else:
						query_NCBI(seq_id)
				
				annot_flag = 0
				for dir in os.listdir('GenomicCDS/' + seq_id.replace('-', '')):
					if(os.path.isfile('GenomicCDS/' + seq_id.replace('-', '') + '/' + dir + '/annotated_cds.fasta')):
						annot_flag = 1
				if(annot_flag == 0):
					print(seq_id)
					
					if(veupath):
						query_NCBI(seq_id)
						annot_flag = 0
						for dir in os.listdir('GenomicCDS/' + seq_id.replace('-', '')):
							if(os.path.isfile('GenomicCDS/' + seq_id.replace('-', '') + '/' + dir + '/annotated_cds.fasta')):
								annot_flag = 1
						if(annot_flag == 0):
							print('\nNo annotated CDS were retrieved for OG ' + og + '. Skipping this one for now, it may be worth trying to download this genome manually if possible\n')
							sings.write(seq_id + ',No Database Sequence Found\n')
							os.system('rm -r GenomicCDS/' + seq_id.replace('-', ''))
					else:
						print('\nNo annotated CDS were retrieved for OG ' + og + '. Skipping this one for now, it may be worth trying to download this genome manually if possible\n')
						sings.write(seq_id + ',No Database Sequence Found\n')
						os.system('rm -r GenomicCDS/' + seq_id.replace('-', ''))
				else:
					sings.write(seq_id + ',Success\n')
					print('\nAnnotated CDS successfully retrieved for OG ' + og + '\n')
			elif(seq_id[:10] in all_genomic_taxa):
				sings.write(seq_id + ',Formatting\n')

def blast_singleton_seqs():
		
	for dir in os.listdir('GenomicCDS'):
		if('ds_store' not in dir.lower()):
			for seqid in os.listdir('GenomicCDS/' + dir):
				if('ds_store' not in seqid.lower() and os.path.isdir('GenomicCDS/' + dir + '/' + seqid)):
								
					os.system('makeblastdb -in GenomicCDS/' + dir + '/' + seqid + '/annotated_cds.fasta -parse_seqids -dbtype nucl -title ' + dir + ' -out GenomicCDS/' + dir + '/' + seqid + '/' + seqid + '_db')
					
					with open('GenomicCDS/' + dir + '/singleton_seq.fasta') as o:
						if(seqid[6:10].lower() == seqid[6:10]):
							for og_file in os.listdir('../../allOG5Files'):
								if(seqid[:10] in og_file):
									prot = [record.seq for record in SeqIO.parse('../../allOG5Files/' + og_file, 'fasta') if record.id == seqid.replace('LKH', '-LKH')][0]	
									o.write('>' + seqid + '\n' + str(prot) + '\n\n')
						else:
							for rtg_file in os.listdir('../../ReadyToGo_AA'):
								if(seqid[:10] in rtg_file):
									prot = [record.seq for record in SeqIO.parse('../../ReadyToGo_AA/' + rtg_file, 'fasta') if record.id == seqid.replace('LKH', '-LKH')][0]	
									o.write('>' + seqid + '\n' + str(prot) + '\n\n')
			
					os.system('tblastn -query GenomicCDS/' + dir + '/singleton_seq.fasta -evalue 1e-10 -max_target_seqs 10 -num_threads 6 -outfmt 6 -db GenomicCDS/' + dir + '/' + seqid + '/' + seqid + '_db -out GenomicCDS/' + dir + '/' + seqid + '/' + seqid + '_blasted.tsv')
					
					
def summarise_hits():

	blasted_info = { }

	for dir in os.listdir('GenomicCDS'):
		if(os.path.idir('GenomicCDS/' + dir)):
			pass
								
									
def main():

	if(not os.path.isdir('GenomicCDS')):
		os.mkdir('GenomicCDS')	
	
	#get_accession_numbers()
		
	#make_DBs()
			
	#blast_singleton_seqs()
	
	summarise_hits()
		
	
	
main()

