#Dependencies
from p4 import *
import os, re
import sys
import csv
from Bio import SeqIO
import subprocess


def bad_script_call():

	print('\nTry inputting alternative tree building methods --iqtree, --mrbayes, or both. If neither are used, the script will simply reorganize the "Subtrees" directory\n')

def directory_error():
	
	print('\nThe "Subtrees" directory could not be found or is improperly organized. This should be a sister directory to the "Scripts" folder, and should contain the folders "Trees", "Preguidance", and "Postguidance" folders. The "Postguidance" folder should contain all new RaxML subtrees and their pre- and postguidance files.\n')
	exit()

def get_args():

	use_iqtree = False
	use_mrbayes = False
	
	try:
		if(sys.argv[1] == '--iqtree'):
			use_iqtree = True
		elif(sys.argv[1] == '--mrbayes'):
			use_mrbayes = True
		else:
			print('\nInvalid tree building method.\n')
			bad_script_call()
			return [ False, False ]
	except IndexError:
		print('\nNo input tree building methods.\n')
		bad_script_call()
		return [ False, False ]
		
	try:
		if(sys.argv[2] == '--iqtree'):
			use_iqtree = True
		elif(sys.argv[2] == '--mrbayes'):
			use_mrbayes = True
		else:
			print('\nInvalid tree building method.\n')
			bad_script_call()
			return [ False, False ]
	except IndexError:
		bad_script_call()
	
	return [ use_iqtree, use_mrbayes ]
	
	
def realign():
	
	for og in os.listdir('../Subtrees'):
		if('ds_store' not in og.lower()):
			for file in os.listdir('../Subtrees/' + og + '/Originals'):
				if('preguidance' in file):
					os.system('mafft --auto --inputorder ../Subtrees/' + og + '/Originals/' + file + ' > ../Subtrees/' + og + '/Realigned/' + file.split('preguidance')[0] + 'Realigned.fas')


def iqtree():

	for og in os.listdir('../Subtrees'):
		if('ds_store' not in og.lower()):
			os.system('mkdir ../Subtrees/' + og + '/IQTree')
			for file in os.listdir('../Subtrees/' + og + '/Realigned'):
				if('Realigned' in file):
					os.system('iqtree -s ../Subtrees/' + og + '/Realigned/' + file + ' -nt 4 -bb 1000 -pre ../Subtrees/' + og + '/IQTree/' + file.split('Realigned')[0] + '_iqtree')
		
	
def mrbayes():

	for og in os.listdir('../Subtrees'):
		if('ds_store' not in og.lower()):
			#os.system('mkdir ../Subtrees/' + og + '/MrBayes')
			for file in os.listdir('../Subtrees/' + og + '/Realigned'):
				if('Realigned' in file):
					#Getting the model to use for MrBayes using the IQTree ModelFinder
					#os.system('iqtree -s ../Subtrees/' + og + '/Realigned/' + file + ' -m MF -mset mrbayes -nt 4 -pre ../Subtrees/' + og + '/MrBayes/' + file.split('Realigned')[0] + 'mrbayes_model')
					
					model = ''
					for line in open('../Subtrees/' + og + '/MrBayes/' + file.split('Realigned')[0] + 'mrbayes_model.iqtree'):
						if('best-fit model according to' in line.lower()):
							model = line.split(':')[1].strip()
							break
							
					if(model == ''):
						print('\nNo model returned for OG ' + og + '\n')
					else:
						subtree_nex_path = '../Subtrees/' + og + '/MrBayes/' + file.split('Realigned')[0] + 'mrbayes_subtree_nex.tre'
						mrbayes_run_path = '../Subtrees/' + og + '/MrBayes/' + file.split('Realigned')[0] + 'mrbayes_run_file.txt'
					
						records = [(record.id, record.seq) for record in list(SeqIO.parse('../Subtrees/' + og + '/Realigned/' + file, 'fasta'))]
					
						with open(subtree_nex_path, 'w') as o:
							o.write('#NEXUS\nbegin taxa;\n\tdimensions ntax=' + str(len(records)) + ';\n\ttaxlabels\n')
							for seq in records:
								o.write('\t\t' + seq[0] + '\n')
							o.write(';\nend;\n\nbegin data;\n\tdimensions ntax=' + str(len(records)) + ' nchar=' + str(len(records[0][1])) + ';\n')
							o.write('\tformat datatype=protein gap=-;\nmatrix\n')
							for seq in records:
								o.write('\n' + seq[0] + '\n' + str(seq[1]) + '\n')
							o.write(';\nend;')
					
						with open(mrbayes_run_path, 'w') as o:
							o.write('begin mrbayes;\n')
							o.write('\tset autoclose=yes nowarn=yes;\n')
							o.write('\texecute ' + subtree_nex_path + ';\n')
							o.write('\tlset applyto=(all) nst=6 rates=gamma;\n')
							o.write('\tprset applyto=(all) aamodelpr=fixed(' + model + ') ratepr=variable;\n')
							o.write('\tunlink shape=(all) pinvar=(all) statefreq=(all) revmat=(all);\n')
							o.write('\tmcmcp ngen=5500000 printfreq=1000 samplefreq=100;\n')
							#o.write('\tnchains=4 savebrlens=yes;\n')
							o.write('\tmcmc file=' + '../Subtrees/' + og + '/MrBayes/' + file.split('Realigned')[0] + '_mrbayes;\n')
							o.write('\tsumt filename=' + '../Subtrees/' + og + '/MrBayes/' + file.split('Realigned')[0] + '_mrbayes' + ' burnin=5000 contype=halfcompat;\nend;')
						
						os.system('mb ' + mrbayes_run_path + ' > ' + mrbayes_run_path.split('run_file')[0] + 'log_file.txt')
					
					exit()
					break
					
	
def organize_files():
	
	if(os.path.isdir('../Subtrees/Subtrees') and os.path.isdir('../Subtrees/PreGuidance') and os.path.isdir('../Subtrees/Postguidance')):
		
		failed = []
		
		for tree in os.listdir('../Subtrees/Subtrees'):
			og = 'OG5_' + tree.split('OG5_')[-1][:6]
			preguidance = ''
			postguidance = ''
		
			for pg_file in os.listdir('../Subtrees/Preguidance'):
				if(og in pg_file):
					preguidance = '../Subtrees/Preguidance/' + pg_file
					
			for new_file in os.listdir('../Subtrees/Postguidance'):
				if(og in new_file and 'guidance' in new_file and '95gapTrimmed' not in new_file and '.tre' not in new_file):
					postguidance = '../Subtrees/Postguidance/' + new_file
												
			if(preguidance != '' and postguidance != ''):
				os.system('mkdir ../Subtrees/' + og)
				os.system('mkdir ../Subtrees/' + og + '/Originals')
				os.system('mkdir ../Subtrees/' + og + '/Realigned')
				
				os.system('mv ../Subtrees/Trees/' + tree + ' ../Subtrees/' + og + '/Originals/' + tree)
				os.system('mv ' + preguidance + ' ../Subtrees/' + og + '/Originals/' + preguidance.split('/')[-1])
				os.system('mv ' + postguidance + ' ../Subtrees/' + og + '/Originals/' + postguidance.split('/')[-1])
			else:
				print(og + ' FAILED -- Critical files could not be found\n')
				failed.append(og)
				
		if(len(failed) == 0):
			os.system('rm -r ../Subtrees/Subtrees')
			os.system('rm -r ../Subtrees/Preguidance')
			os.system('rm -r ../Subtrees/Postguidance')
	else:
		directory_error()


def main():

	use_iqtree, use_mrbayes = get_args()

	#organize_files()
	
	#realign()
	
	if(use_iqtree):
		iqtree()
		
	if(use_mrbayes):
		mrbayes()
	

main()
