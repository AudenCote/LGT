#!/usr/bin/env python3
# coding=utf-8

#Author: Auden Cote-L'Heureux
#Contact at acotelheureux@smith.edu, audenemil@gmail.com, or on Slack if a Katzlabber.
#Retrieves and wrangles GO term data from PFam domains & sequence identifers using the InterPro2GO online database: http://current.geneontology.org/ontology/external2go/interpro2go (I downloaded this and saved it as "pfams_to_go_terms_database.txt")
#Last updated 05/20/21
#Notes: should technically be able to run with python 2 or 3


import os
import sys
from goatools.obo_parser import GODag
from goatools.godag.go_tasks import get_go2parents
from tqdm import tqdm

#You can read about GODag and other GOATools here: https://link.springer.com/protocol/10.1007/978-1-4939-3743-1_16 (or at their GitHub, which is pretty good)
godag = GODag('go-basic.obo',
              optional_attrs={'relationship'})
optional_relationships = set()
go2parents_isa = get_go2parents(godag, optional_relationships)
                        

#Given a PFam domain id, return any corresponding GO terms
def get_go(pfam):
	
	raw_gos = []
	for line in open('pfams_to_go_terms_database.txt'):
		if(pfam in line):
			raw_gos.append(line.split('; ')[-1].strip())
	
	return raw_gos
	
	
#Given a list of GO terms, get all of the GO terms that are above each GO term in the hierarchy
def up_hierarchy(raw_gos):

	go_list = raw_gos
	
	for go in raw_gos:
		while 1 == 1:
			try:
				parents = go2parents_isa[go]
			except KeyError:
				break
			for parent in parents: go_list.append(parent)
			go = parent
		
	return list(dict.fromkeys(go_list))
	

#Wrapper function to wrangle all of the data, calls above functions
def convert_pfams():
	
	for og_file in tqdm(os.listdir('../Unique_Pfams')):
		if(os.path.isfile('../Unique_Pfams/' + og_file) and 'ds_store' not in og_file.lower()):
			og = 'OG5_' + og_file.split('OG5_')[-1][:6]
			with open('../GoTermsFromPFams/' + og + '.tsv', 'w') as o:
				o.write('Seq\tGO Terms\n')
				seen_seqs = ['ORF', 'Protein']
				for l1 in open('../Unique_Pfams/' + og_file):
					seq = l1.split('\t')[0]
					if(seq not in seen_seqs):
						seen_seqs.append(seq)
						pfams = []
						for l2 in open('../Unique_Pfams/' + og_file):
							if(seq in l2):
								pfams.append(l2.split('\t')[1])
													
						raw_gos = []
						for pfam in pfams:
							raw_gos += get_go(pfam)
							
						go_list = up_hierarchy(raw_gos)
						
						o.write(seq + '\t' + ','.join(go_list) + '\n')
						
						
def main():
	
	if(not os.path.isdir('../GOTermsFromPFams')):
		os.mkdir('../GOTermsFromPFams')
	
	convert_pfams()
		
	
main()