import os
import re
import sys
from tqdm import tqdm
from p4 import *
from Bio import SeqIO
from collections import Counter
import xml.etree.ElementTree
from time import sleep


var.doRepairDupedTaxonNames = 1

URL="https://cipresrest.sdsc.edu/cipresrest/v1"
CRA_USER=""
PASSWORD=""
KEY=""


def help():

	print('Invalid input')
	exit()


def get_args():

	if('--ingroup' in sys.argv or '-i' in sys.argv):
		if('--ingroup' in sys.argv):
			try:
				ingroup = sys.argv[sys.argv.index('--ingroup') + 1]
				ingroup_thresh = int(sys.argv[sys.argv.index('--ingroup') + 2])
			except:
				help()
		elif('-i' in sys.argv):
			try:
				ingroup = sys.argv[sys.argv.index('-i') + 1]
				ingroup_thresh = int(sys.argv[sys.argv.index('-i') + 2])
			except IndexError:
				help()
	else:
		help()
		
	if('--outgroup' in sys.argv or '-o' in sys.argv):
		if('--outgroup' in sys.argv):
			try:
				outgroup = sys.argv[sys.argv.index('--outgroup') + 1]
				outgroup_thresh = int(sys.argv[sys.argv.index('--outgroup') + 2])
			except:
				help()
		elif('-o' in sys.argv):
			try:
				outgroup = sys.argv[sys.argv.index('-o') + 1]
				outgroup_thresh = int(sys.argv[sys.argv.index('-o') + 2])
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
		
	diverse_vtg_outgroups = False
	if('--match_ltg_outgroups' in sys.argv):
		pass
	elif('--diverse_vtg_outgroups' in sys.argv):
		diverse_vtg_outgroups = True
		
	if(conserved_dir.endswith('/')):
		conserved_dir = conserved_dir[:-1]
		
	if(xgt_dir.endswith('/')):
		xgt_dir = xgt_dir[:-1]
		
	if(os.path.isdir(conserved_dir) and os.path.isdir(xgt_dir)):
		return [ingroup, ingroup_thresh, outgroup, outgroup_thresh, xgt_dir, conserved_dir, diverse_vtg_outgroups]
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


def get_best_clade(tree, ingroup, ingroup_thresh, constraint_taxa, contam_threshold = 1, nodes_to_exclude = []):

	best_node = None
	best_size = 0
	#contam_threshold = 1
	#Setting this to zero because filtering occurs after the data has been collected; if you would like to filter for best node size raise this value (e.g. absolute_n_threshold = 3)
	absolute_n_threshold = 0
	most_sisters = -1
	
	forbidden_nodes = [node for node in nodes_to_exclude]
	for node in nodes_to_exclude:
		for num in tree.getNodeNumsAbove(node):
			forbidden_nodes.append(tree.node(num))

	for node in tree.iterNodesNoRoot():
		if(len(tree.getAllLeafNames(node)) > absolute_n_threshold and node not in nodes_to_exclude):
			leaves = tree.getAllLeafNames(node)
							
			num = 0.0; dem = 0.0
			
			seen_tax = []
			for leaf in leaves:
				if(ingroup in leaf[:10] and leaf[:10] not in seen_tax):
					if(constraint_taxa == []):
						dem += 1.0
					elif(leaf[:10] in constraint_taxa):
						dem += 1.0
				elif(leaf[:10] not in seen_tax):
					dem += 1.0; num += 1.0;
				seen_tax.append(leaf[:10])
					
			if(constraint_taxa == []):
				if(num <= contam_threshold and dem - num > best_size):
					best_node = node
					best_size = dem - num
			else:
				if(num == 0 and (int(dem - num) == len(constraint_taxa) or int(dem - num) == len(constraint_taxa) - 1) and dem - num >= best_size):
					p_node = node.parent; last_p_node = node
					sisters = 0; out = 0
					while p_node != None:
						p_sisters = 0
						for sis_leaf in tree.getAllLeafNames(p_node):
							if(ingroup in sis_leaf):
								p_sisters += 1
							else:
								out += 1
						if(out > 0):
							break
						else:
							sisters += p_sisters
							last_p_node = p_node
							p_node = p_node.parent
					
					if(dem - num > best_size):
						best_node = node
						best_size = dem - num
					elif(dem - num == best_size and sisters > most_sisters):
						best_node = node
						most_sisters = sisters

	return best_node
			
			
def process_LTG(og, xgt_dir, fname, ingroup, ingroup_thresh, outgroup, outgroup_thresh):

	pre_flag = 0
	for pre_f in os.listdir(xgt_dir):
		if(og in pre_f and 'preguidance' in pre_f):
			pre_flag = 1
			break
						
	if(pre_flag == 0):
		print('\nNo LTG pre-Guidance file found for ' + og + '\n')
		return None

	nexus_to_newick(fname)
	var.trees = []													
	correct_tree(fname)
	tree_file = fname.split('.tre')[0] + '_temp.tre'
	read(tree_file) 
	tree = var.trees[0]			
	os.remove(tree_file)
	#root_by_major_clade(tree)
	
	ingroup_node = get_best_clade(tree, ingroup, ingroup_thresh, [])
	if(ingroup_node != None):
		ingroup_leaves = [leaf[:10] for leaf in tree.getAllLeafNames(ingroup_node) if ingroup in leaf[:10]]
		outgroup_taxa = [leaf[:10] for leaf in tree.getAllLeafNames(0) if outgroup in leaf[:10]]
		outgroup_leaves = [leaf for leaf in tree.getAllLeafNames(0) if outgroup in leaf[:10]]
				
		if(len(dict.fromkeys(outgroup_taxa)) >= outgroup_thresh):
			ingroup_tree = tree.dupeSubTree(ingroup_node, up = True)
	
			paralog_sbl = { }
		
			for node in ingroup_tree.iterNodesNoRoot():
				if(node.name != None):
					if(node.name[:10] in ingroup_leaves):
						if(ingroup_leaves.count(node.name[:10]) > 1):
							if(node.name[:10] not in paralog_sbl):
								paralog_sbl.update({ node.name[:10] : node })
							elif(paralog_sbl[node.name[:10]].br.len > node.br.len):
								paralog_sbl[node.name[:10]] = node
						else:
							paralog_sbl.update({ node.name[:10] : node })
						
			ingroup_nodes_to_keep = paralog_sbl.values()
			ingroup_nodes_to_remove = []
			for node in ingroup_tree.iterNodesNoRoot():
				if(node not in ingroup_nodes_to_keep and node.getNChildren() == 0):
					ingroup_nodes_to_remove.append(node.name)
				
			for nodename in ingroup_nodes_to_remove:
				ingroup_tree.removeNode(nodename, alsoRemoveSingleChildParentNode = False)
			
			taxa_to_use = [node[:10] for node in ingroup_tree.getAllLeafNames(0)] + list(dict.fromkeys(outgroup_taxa))
			leaves_to_use = [node for node in ingroup_tree.getAllLeafNames(0)] + list(dict.fromkeys(outgroup_leaves))
			
			with open('FilteredUnaligned/' + og + '_FilteredLTG.fasta', 'w') as o:
				for rec in SeqIO.parse(xgt_dir + '/' + pre_f, 'fasta'):
					if(rec.description in [leaf.replace('LKH', '-LKH').replace('--LKH', '-LKH') for leaf in leaves_to_use]):
						o.write('>' + rec.description + '\n' + str(rec.seq) + '\n\n')
				
			return taxa_to_use
		else:
			return None
	else:
		return None
			
			
def process_VTGs(ltg_og, ingroup, outgroup, constraint_taxa, conserved_dir, diverse_vtg_outgroups):
	
	vtg_ogs = []
	for tree_f in os.listdir(conserved_dir):
		if((tree_f.endswith('.tre') or tree_f.endswith('.tree')) and 'OG5_' in tree_f):
			fname = conserved_dir + '/' + tree_f
			og = 'OG5_' + fname.split('OG5_')[1][:6]
			
			pre_flag = 0
			for pre_f in os.listdir(conserved_dir):
				if('preguidance' in pre_f and og in pre_f):
					pre_flag = 1
					break
			
			if(pre_flag == 0):
				print('\nNo VTG pre-Guidance file found for ' + og + '\n')
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
			
			constraint_ingroups = [tax for tax in constraint_taxa if ingroup in tax]
			best_monophyletic_ingroup_node = get_best_clade(tree, ingroup, 0, constraint_ingroups)
			
			if(best_monophyletic_ingroup_node != None):
				print('\n\tFound a monophyletic ingroup clade in VTG ' + og + '\n')
				vtg_ogs.append(og)
						
				constraint_outgroups = [tax for tax in constraint_taxa if outgroup in tax]
				selected_outgroups = []
				if(diverse_vtg_outgroups):	
					print('\n\t\tFinging the most diverse VTG outgroup...\n')
					print('\n\t\tGetting patristic distance matrix for ' + og + '\n')
					if(not os.path.isfile('VTGPatristicMatrices/' + og + '.nex')):
						tree.taxNames = tree.getAllLeafNames(0)
						dm = tree.patristicDistanceMatrix()
						dm.writeNexus('VTGPatristicMatrices/' + og + '.nex')
					else:
						dm = DistanceMatrix()
						dm.matrix = []; dm.names = []
					
						file_lines = [line for line in open('VTGPatristicMatrices/' + og + '.nex')]
						for i, l1 in enumerate(file_lines):
							if(l1.strip().startswith('[')):
								l2 = file_lines[i + 1]
								l = i + 1
								while ';' not in l2:
									dm.matrix.append([val.strip() for val in l2.strip().split(' ') if val.strip() != ''][1:])
									dm.names.append(l2.strip().split(' ')[0])
									l += 1
									l2 = file_lines[l]
				
					for node in tree.iterNodesNoRoot():
						if(node.name in list(tree.getAllLeafNames(best_monophyletic_ingroup_node)) and node.name[:10] in constraint_ingroups):
							benchmark_ingroup_leaf = node.name; ingroup_benchmark_dist = 0
							while(node.nodeNum != best_monophyletic_ingroup_node.nodeNum):
								ingroup_benchmark_dist += node.br.len
								node = node.parent
						
							break
				
					benchmark_ingroup_matrix_row = dm.names.index(benchmark_ingroup_leaf)
				
					mc_counts = Counter([leaf[:5] for leaf in tree.getAllLeafName(0) if outgroup in leaf[:10]])
					mcs_ranked = sorted(mc_counts, key = mc_counts.get, reverse = True)
					
					while len(selected_outgroups) < len(constraint_outgroups) and len(selected_outgroups) < len([leaf[:5] for leaf in tree.getAllLeafName(0) if outgroup in leaf[:10]]):
						for mc in mcs_ranked:
							min_dist = float('inf'); best_leaf = None
							for leaf in tree.getAllLeafName(0):
								outgroup_matrix_col = dm.names.index(leaf)
								if(leaf not in selected_outgroup and dm.matrix[benchmark_ingroup_matrix_row][outgroup_matrix_col] - ingroup_benchmark_dist < min_dist):
									min_dist = dm.matrix[benchmark_ingroup_matrix_row][outgroup_matrix_col] - ingroup_benchmark_dist
									best_leaf = leaf
							if(best_leaf != None):
								selected_outgroups.append(best_leaf)
				else:
					print('\n\t\tMatching the LTG outgroup distribution as closely as possible...\n')
					for tax in constraint_outgroups:
						for leaf in tree.getAllLeafNames(0):
							if(tax in leaf):
								selected_outgroups.append(leaf)
					
				with open('FilteredUnaligned/' + og + '_ConstrainedLeavesVTG_' + ltg_og + '_LTG.fasta', 'w') as o:	
					for rec in SeqIO.parse(conserved_dir + '/' + pre_f, 'fasta'):
						if(rec.description in [leaf.replace('LKH', '-LKH').replace('--LKH', '-LKH') for leaf in [l for l in list(tree.getAllLeafNames(best_monophyletic_ingroup_node)) if l[:10] in constraint_ingroups] + selected_outgroups]):
							o.write('>' + rec.description + '\n' + str(rec.seq) + '\n\n')
						
	return vtg_ogs
	
	
def get_completed_Cipres_runs():

	os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL + '/job/' + CRA_USER + ' > status_check_all_temp.xml')

	jobids = []
	for line in open('status_check_all_temp.xml'):
		if('<title>NGBW-JOB-' in line):
			jobids.append(line.split('>')[1].split('<')[0].strip())
	
	completed = []
	for jobid in jobids:
		os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL + '/job/' + CRA_USER + '/' + jobid + ' > ' + jobid + '_status_temp.xml')
		for line in open(jobid + '_status_temp.xml'):
			if('<stage>COMPLETED</stage>' in line):
				completed.append(jobid)
				break
				
		os.remove(jobid + '_status_temp.xml')
		sleep(90)
	
	running = len([jobid for jobid in jobids if jobid not in completed])
				
	return completed, running
	
	
def delete_Cipres_runs(exceptions = [], mode = 'completed'):

	os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL + '/job/' + CRA_USER + ' > deletion_temp.xml')
	
	jobids = []
	for line in open('deletion_temp.xml'):
		if('<title>NGBW-JOB-' in line):
			jobids.append(line.split('>')[1].split('<')[0].strip())
	
	for jobid in jobids:
		if(jobid not in exceptions):
			if(mode == 'completed' and jobid in get_completed_Cipres_runs()):
				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" -X DELETE "+URL + '/job/' + CRA_USER + '/' + jobid)
			elif(mode == 'all'):
				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" -X DELETE "+URL + '/job/' + CRA_USER + '/' + jobid)
				
			sleep(90)
			
	os.remove('deletion_temp.xml')
	
	
def align_and_tree():
	print('\nNow aligning & building all selected trees...\n')

	if(not os.path.isdir('ReducedAligned')):
		os.mkdir('ReducedAligned')
		
	if(not os.path.isdir('ReducedTrees')):
		os.mkdir('ReducedTrees')
		
	jobid_backup_file = open('Cipres_JobID_Backup.csv', 'w')
	jobid_backup_file.write('File,Job.ID\n')
	
	file_by_jobid = { }
	for file in os.listdir('FilteredUnaligned'):
		if(file != '.DS_Store'):
			print('\tSubmitting ' + file + ' to MAFFT on Cipres\n')
			os.system("curl -k -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F metadata.statusEmail=false -F tool=MAFFT_XSEDE -F input.infile_=@./FilteredUnaligned/" +file+" -F vparam.datatype_='protein' -F vparam.outputOrder_='--reorder' -F vparam.which_mafft_='7471' -F vparam.runtime_='168' -F vparam.more_memory_='1' > " + file + "_MAFFT_runinit_temp.xml")
			
			job_id = ''	
			lines = []	
			for line in open(file + "_MAFFT_runinit_temp.xml"):
				lines.append(line)
				if('<title>NGBW-JOB-' in line):
					job_id = line.split('>')[1].split('<')[0].strip()
					break
			
			if(job_id != ''):
				jobid_backup_file.write(file + ',' + job_id + '\n')
				file_by_jobid.update({ job_id : file })
				print('\t\tSuccessfully retrieved Job ID: ' + job_id + '\n')
			else:
				jobid_backup_file.write(file + ',Failed\n')
				print('\t\tMAFFT failed on Cipres for file ' + file + '\n')
			
			os.remove(file + "_MAFFT_runinit_temp.xml")
			#break ####TEMP####
			
	print('Now starting status check & IQ-Tree loop\n')
	running = True; r = 0; downloaded = []; flag = 0
	while running == True or flag == 1:
		r += 1
		sleep(600)
		
		completed, n_running = get_completed_Cipres_runs()
		print('Loop ' + str(r) + ': ' + str(len(completed)) + ' jobs have finished and ' + str(n_running) + ' are still running on Cipres\n')
		if(n_running == 0): 
			running = False
		else:
			running = True
		
		completed = [jobid for jobid in completed if jobid in file_by_jobid and jobid not in downloaded]
		failed_to_keep = []
		
		for jobid in completed:
			print(jobid)
			os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" -O -J "+URL + '/job/' + CRA_USER + '/' + jobid + '/output')
			output_lines = [line for line in open('output')]
			if('MAFFT' in jobid):
				print('MAFFT job ' + file_by_jobid[jobid] + ' (' + jobid + ') completed. Checking output now\n')
				outfile_name = None
				for l, line in enumerate(output_lines):
					if('<filename>output.mafft</filename>' in line):
						outfile_num = output_lines[l-5].split('<')[1].split('/')[-1]
						if('VTG' in file_by_jobid[jobid]):
							outfile_name = 'ReducedAligned/' + file[:10] + '_ReducedAlignedVTG_' + file[-20:]
						else:
							outfile_name = 'ReducedAligned/' + file[:10] + '_ReducedAlignedLTG.fasta'
						
						os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" -O -J "+URL + '/job/' + CRA_USER + '/' + jobid + '/output/' + outfile_num)
						os.system('mv output.mafft ' + outfile_name)
						break
						
				if(outfile_name != None):
					downloaded.append(jobid)
					print('MAFFT output successfully retrieved for ' + file_by_jobid[jobid] + ' (' + jobid + '). Now running IQ-Tree on this alignment\n')
					sleep(90)
					os.system("curl -k -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F metadata.statusEmail=false -F tool=IQTREE_XSEDE_EXPANSE -F vparam.bootstrap_type_=bb -F vparam.specify_runtype_=2 -F vparam.specify_numpatterns_=1 -F input.infile_=@./" + outfile_name + " -F vparam.sequence_type_=AA -F vparam.specify_protmodel_=LG -F vparam.specify_modelrate_=+G -F vparam.num_bootreps_=1000 -F vparam.which_iqtree_=212 > " + outfile_name.split('/')[-1] + "_IQTREE_runinit_temp.xml")
					
					job_id = ''	
					lines = []	
					for line in open(outfile_name.split('/')[-1] + "_IQTREE_runinit_temp.xml"):
						lines.append(line)
						if('<title>NGBW-JOB-' in line):
							job_id = line.split('>')[1].split('<')[0].strip()
							break
			
					if(job_id != ''):
						jobid_backup_file.write(outfile_name + ',' + job_id + '\n')
						file_by_jobid.update({ job_id : outfile_name })
						print('IQTREE successfully run on Cipres for file ' + outfile_name + ', Job ID ' + job_id + '\n')
					else:
						jobid_backup_file.write(outfile_name + ',Failed\n')
						failed_to_keep.append(job_id)
						print('IQTREE failed on Cipres for file ' + outfile_name + '\n')
			
					os.remove(outfile_name.split('/')[-1] + "_IQTREE_runinit_temp.xml")
					
				else:
					print('No MAFFT output file found on Cipres for ' + file_by_jobid[jobid] + ', Job ID ' + jobid + '. Adding to list of failed jobs\n')
					failed_to_keep.append(jobid)
			elif('IQTREE' in jobid):
				print('iqtree in jobid ', jobid)
				print('IQ-Tree job ' + file_by_jobid[jobid] + ' (' + jobid + ') completed. Checking output now\n')
				outfile_name = None
				for l, line in enumerate(output_lines):
					if('<filename>output.treefile</filename>' in line):
						outfile_num = output_lines[l-5].split('<')[1].split('/')[-1]
						if('VTG' in file_by_jobid[jobid]):
							outfile_name = 'ReducedTrees/' + file[:10] + '_ReducedTreeVTG_' + file.split('.fasta')[0][-14:] + '.tre'
						else:
							outfile_name = 'ReducedTrees/' + file[:10] + '_ReducedTreeLTG.tre'
													
						print('Retrieving IQ-Tree output from Cipres for ' + file_by_jobid[jobid] + ', Job ID ' + jobid + '\n')
						os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" -O -J "+URL + '/job/' + CRA_USER + '/' + jobid + '/output/' + outfile_num)
						os.system('mv output.treefile ' + outfile_name)
						break
						
				if(outfile_name != None):
					print('IQ-Tree output successfully retrieved (' + jobid + ').\n')
					downloaded.append(jobid)
				else:
					print('IQ-Tree output file not found (' + jobid + ').\n')
					failed_to_keep.append(jobid)
				
				exit()
			
			os.remove('output')
		
		with open('Cipres_failed_Jobs.txt', 'w') as f:
			for jobid in failed_to_keep:
				f.write(jobid + '\n')
				
		print('End of loop ' + str(r)+ '\n')
		
		completed, n_running = get_completed_Cipres_runs()
		if(n_running == 0): 
			running = False
		else:
			running = True
			
		if(flag == 1):
			break
		elif(flag == 0 and running == False):
			flag = 1
			
		print('loop end', completed, n_running)
		print(downloaded)
		
	print('\nStatus check and IQ-Tree loop finished -- all output files can be found in the ReducedAligned and ReducedTrees folders. Deleting all completed runs.\n')
	delete_Cipres_runs(failed_to_keep, mode = 'completed')
	
	jobid_backup_file.close()
		
		
def calc_ratio(fname, ingroup, outgroup):

	og = fname.split('/')[-1][:10]

	nexus_to_newick(fname)
	var.trees = []													
	correct_tree(fname)
	tree_file = fname.split('.tre')[0] + '_temp.tre'
	read(tree_file) 
	tree = var.trees[0]			
	os.remove(tree_file)
	#root_by_major_clade(tree)
	
	
	ingroup_node = get_best_clade(tree, ingroup, 0, [], 0)
	ingroup_tree = tree.dupeSubTree(ingroup_node, up = True)
	
	avg_ingroup_bl = float(ingroup_tree.getLen() - ingroup_node.br.len)/float(len(list(ingroup_tree.iterNodesNoRoot())) - 1)
		
	tree.reRoot(ingroup_node)
	ingroup_node = get_best_clade(tree, ingroup, 0, [])
	
	print(ingroup_tree.getAllLeafNames(0))
	print(float(ingroup_tree.getLen()) - ingroup_node.br.len)
	print(float(len(list(ingroup_tree.iterNodesNoRoot()))))
	print(avg_ingroup_bl)
			
	print('\n\tGetting patristic distance matrix for ' + og + '\n')
	tree.taxNames = [node.name for node in tree.iterNodesNoRoot() if node.getNChildren() == 0]
	dm = tree.patristicDistanceMatrix()
	
	for node in tree.iterNodesNoRoot():
		if(node.name in list(tree.getAllLeafNames(ingroup_node))):
			benchmark_ingroup_leaf = node.name; ingroup_benchmark_dist = 0
			while node.nodeNum != ingroup_node.nodeNum:
				ingroup_benchmark_dist += node.br.len
				node = node.parent
			ingroup_benchmark_dist += node.br.len
			
			break
	
	benchmark_ingroup_matrix_row = dm.names.index(benchmark_ingroup_leaf)
	
	outgroup_node = get_best_clade(tree, outgroup, 0, [])
	for node in tree.iterNodesNoRoot():
		if(node.name in list(tree.getAllLeafNames(outgroup_node))):
			benchmark_outgroup_leaf = node.name; outgroup_benchmark_dist = 0
			while node.nodeNum != outgroup_node.nodeNum:
				outgroup_benchmark_dist += node.br.len
				node = node.parent
		
			break
	
	benchmark_outgroup_matrix_col = dm.names.index(benchmark_outgroup_leaf)
	dist_to_ingroup = dm.matrix[benchmark_ingroup_matrix_row][benchmark_outgroup_matrix_col] - ingroup_benchmark_dist - outgroup_benchmark_dist
	
	print(dist_to_ingroup)
	
	print('\tDone calculating branch length ratios for ' + og + '\n')
	return len(list(ingroup_tree.iterNodesNoRoot())) - 1, avg_ingroup_bl, dist_to_ingroup

							
def main():

	print('')

	ingroup, ingroup_thresh, outgroup, outgroup_thresh, xgt_dir, conserved_dir, diverse_vtg_outgroups = get_args()
	
	lgt_ogs = ['OG5_' + file.split('OG5_')[-1][:6] for file in os.listdir(xgt_dir) if file.endswith('.tre')]
	conserved_ogs = ['OG5_' + file.split('OG5_')[-1][:6] for file in os.listdir(conserved_dir) if file.endswith('.tre')]
	
	if(not os.path.isdir('FilteredUnaligned')):
		os.mkdir('FilteredUnaligned')
		
	if(not os.path.isdir('VTGPatristicMatrices')):
		os.mkdir('VTGPatristicMatrices')
	
	print('\nGetting constrained VTG trees for every LTG\n')
	with open('VTGs_per_LTG.csv', 'w') as o:
		o.write('LTG,VTG\n')
		lt = 0
		for ltg_tree in os.listdir(xgt_dir):
			if(ltg_tree.endswith('.tre') and 'OG5_' in ltg_tree):
				lt += 1
				og = 'OG5_' + ltg_tree.split('OG5_')[-1][:6]
				print('LTG ' + str(lt) + ': ' + og)
				
				constraint_taxa = process_LTG(og, xgt_dir, xgt_dir + '/' + ltg_tree, ingroup, ingroup_thresh, outgroup, outgroup_thresh)
				if(constraint_taxa != None):
					vtgs = process_VTGs(og, ingroup, outgroup, constraint_taxa, conserved_dir, diverse_vtg_outgroups)
					
					for vtg in vtgs:
						o.write(og + ',' + vtg + '\n')
		
	print('\nAligning & building reduced trees...\n')															
	align_and_tree()
	
	exit()
	
	print('\nLast step: calculating branch length ratios for ALL FILTERED TREES\n')
	with open('branch_length_ratios.csv', 'w') as o:
		o.write('OG,Type,Ingroup.Clade,N.Ingroups,Avg.Ingroup.Bl,Dist.To.Outgroup\n')
		for tre_f in os.listdir('ReducedTrees'):
			if(tre_f.endswith('.treefile')):
				og = tre_f[:10]
				if(og in lgt_ogs):
					type = 'LGT'
				else:
					type = 'Conserved'
				o.write(og + ',' + type + ',' + ingroup + ',' + ','.join([str(val) for val in calc_ratio('ReducedTrees/' + tre_f, ingroup, outgroup)]) + '\n')
	

main()
#delete_Cipres_runs(mode = 'all')
							
						
