# LGT
Public repository for code used in Cote-L'Heureux et al. 2022, "Old Genes in New Places: A Taxon-Rich Analysis of Interdomain Lateral Gene Transfer Events" (https://journals.plos.org/plosgenetics/article/figures?id=10.1371/journal.pgen.1010239) These scripts are not intended as standalone ready-to-use tools, and are tailored to our database naming system in which each taxon gets a ten-digit code, e.g. Am_di_Acas indicates Acanthamoeba CAStellanii (DIscosea, AMoebozoa).

## Folder structure:
### BranchLengthComparisons
There are two scripts used to conduct the comparison of relative branch-lengths between LTG and VTG trees. LossVsLGT_TreeComparison.py subsamples conserved gene trees (must be provided) to match the taxonomic distribution of each LTG, re-aligns the remaining sequences using MAFFT and constructs gene trees using IQ-Tree through the CIPRES API (user information required). The arguments to the script are as follows:

--ingroup: the major clade of eukaryotes (e.g. "Ex" for Excavata)
--outgroup: the major clade of prokaryotes (e.g. "Ba Za" for bacteria and archaea)
--xgt_dir: a folder of LTG trees and alignments
--conserved_dir: a folder of VTG trees and alignments to be subsampled
--match_ltg_outgroups vs. --diverse_vtg_outgroups: whether the prokaryotes should match the LTG distribution exactly or whether the broadest range of prokaryotes should be kept that matches the number of sequences in the given LTG tree

The second script, BLRatioDistribution.r, is an example of how the dataframe output by LossVsLGT_TreeComparison.py could be processed into plots. This script will require hard-coded customizations depending on the recipients being investigated (and dataframe name, etc.)

### DataCuration
This folder contains the most general scripts used in the curation of genomic and transcriptomic data. Some of this work was done manually in Excel spreadsheets, and these scripts are tools to create dataframes that require manual curation. DoGenomicsHit.py queries 1) NCBI and 2) VEuPathDB for each genomic taxon, downloads the assembly, and BLASTs the LTG sequence against the genome to allow analysis of the genomic context of the putatively transferred gene. We then used the annotations where available to determine the other genes that lay on the hit scaffold, or else BLASTed several thousand base pairs on each end of the scaffold for further curation. A spreadsheet with each genomic sequence to BLAST and the accessions for all genomic taxa must be provided; the spreadsheet names can be found in the script. 

CompareComps.py is largely a wrapper for CUB.py, which contains several functions for compositonal analysis of ORFs. It has the following arguments, and compares the composition of conserved genes to that of putative LTGs to assess the probability of a putative LTG sequence being a contaminant:

--conserved_dir: path to folder of alignments of conserved genes
--xgt_dir: path to folder of alignments of putative LTGs
--non_recip: whether to include analysis of non-eukaryotic sequences

PlotComps.r uses Mahalonobis distance to iteratively determine whether putative LTGs are significantly from the conserved genes (i.e. with distance in the top 1%) and then plots the points with coloration of 1) conserved vs. putatively transferred genes and 2) significantly outlying sequences.

### GeneFamilySelection
The scripts genus_presence.py and minor_clade_presence.py are very similar; both take in a folder of sequence data, with all sequences for each taxon initially included in the study and each sequence assigned an OG (paths hardcoded at the top of the script). They then output a spreadsheet with the number of times each genus (or minor clade, e.g. Discosea) occurs in each OG. These spreadsheets can then be used to rank OGs by exclusive presence -- i.e. the concentration of the OG in a specific genus (or minor clade) of eukaryotes. gene_families_by_major_clade_species.py is very similar, but counts by major clade (e.g. Amoebozoa) and only counts each species once (we have multiple samples from several species in our data). OG_presence_plots.r take as input the output spreadsheets from any of the above scripts (path hard coded in script) and plots the proportion of eukaryotes with the OG that are in a given major clade by the number of prokaryotes with the OG (putative LTGs lie in the top righthand corner of this plot). Several example dataframes are also provided in the Dataframes folder.

### FunctionalCategorization
These scripts all were used in analyzing the functional distributions of groups of LTGs. This was done by first querying EggNOG using FunctionalCategorization.py, which has two arguments:

--input_dir: the folder of input trees & sequences with which to query EggNOG
--slim: whether or not to slim the returned Gene Ontology (GO) terms using the goslim_generic.obo database

The script also includes some functionality to subsample sequences on the tree in a taxonomically informed manner to search to minimize computational time, but this was not applied in the end. The script outputs a database with the presence of Gene Ontology terms across all sequences and edits the input Newick strings to colorize sequences based on the group of GO terms returned. GOFromPFam.py uses GOATools to return and then slim a list of GO terms given a PFAM identifier, which were widely returned by EggNOG. SummaryAndMetadata.py returns some summary statistics of the spreadsheets output by FunctionalCategorization.py, and also aids in the mapping of functional terms onto phylogenetic trees. It does this by acting as a wrapper for the R script MapFunc2Tree.r, which creates a heatmap figure of functional-term presence across each putative LTG tree. SummarizeGOs.py creates summary dataframes using the dataframes output by FunctionalCategorization.py and was only used to aid in development of the functional analysis process and later manual curation stages. query-revigo.r was used as an alternative tool for slimming GO term lists.

### Topology
AUTest/MultiTransfers.py detects non-monophyletic eukaryotic clades and tests for their monophyly by approximately-unbiased (AU) testing using IQ-Tree. This script was used to resolve cases of OGs that have undergone multiple putative lateral transfer events.

CladeAnalysisPipelines contains a series of scripts that were used to analyse the taxonomic composition of gene trees; these scripts are numbered, and generally scripts with later numbers depend on outputs from scripts before. These scripts were used to characterize 

0_GetSubtrees.py: given a eukaryotic group (e.g. Amoebozoa) and a phylogenetic tree, locate and save the sequences in the largest clades that contain only eukaryotes of that clade and meet certain other criteria

1_CladeSizes.py: returns a spreadsheet with a list of clade sizes for each eukaryotic groups for each input gene tree

2_MissingSisters.py: for each OG, clusters all eukaryotic sequences at a given identity and BLASTs one representative sequence from each cluster against the NCBI nr database to detect any eukaryotes not in our database that have the OG (and potentially indicate that there is not actually exclusive presence)

3_TaxonomicAnalysis.py: requires an input spreadsheet with a semi-colon separated taxonomic string for each species (e.g. Eukarya; Amoebozoa; Discosea; Acanthamoeba; Acanthamoeba castellanii) and counts the number of times each taxonomic level appear across all trees (this was used, for example, to add the number of transfers at each taxonomic levels in Opisthokonts to Figure 3).

tree_singleton_identification.py: identifies eukaryotic sequences with no immediate eukaryotic sisters (i.e. "singletons" nested in prokaryotes). These sequences are either very recent LTGs, or, more likely, contamination.
