This directory contains files related to MitoCOGs constructed at Eugene Koonin's group at the National Center for Biotechnology Information (NCBI), National Library of Medicine (NLM), National Institutes of Health (NIH).

#-----------------------------------------------------------------------------
Citation



##############################################################################
1.	Files
##############################################################################

MitoCOGs.txt                                      1107595   Jul  30  2014
mitochondrial_proteomes.txt                        843229   Jul  30  2014
MitoCOGs_sequences.fa                            13147634   Jul  30  2014
MitoCOGs_annotation.txt                              5030   Jul  30  2014
nuclear_encoded_orthologs.txt                       41628   Jul  30  2014
nuclear_proteomes.txt                               12947   Jul  30  2014
nuclear_encoded_orthologs_sequences.fa             457660   Jul  30  2014
intron_gain_loss_estimation_tree.tre                  980   Jul  30  2014
concatenated_intron_presence_absence_matrix.txt     32015   Jul  30  2014



#-----------------------------------------------------------------------------
1.1.	MitoCOGs.txt

Contains the list of orthologous mitochondrial-encoded proteins (MitoCOGs) in a comma-separated format.

<seqID>,<ncbi_taxID>,<protein_length>,<MitoCOG_ID>

* Example:
5819100,9770,347,MitoCOG0001

* Comments:
If the data was downloaded from NCBI, seqID is NCBI protein GI. Otherwise, seqID is a internally assigned ID.

* Example:
intID9000000053,392300,546,MitoCOG0122



#-----------------------------------------------------------------------------
1.2.	mitochondrial_proteomes.txt

Contains list of mitochondrial proteomes (2,486) in a comma-separated format. 

<species_name>,<ncbi_taID>,<ncbi_taxonomy>

* Example:
Jakoba libera,143017,cellular organisms; Eukaryota; Jakobida; Jakobidae; Jakoba



#-----------------------------------------------------------------------------
1.3.	MitoCOGs_sequences.fa

Contains protein sequence of orthologous mitochondrial-encoded proteins (34,740) in FASTA format.

* Example
>gi|251831117|ref|YP_003024036.1| NADH dehydrogenase subunit 5 [Homo sapiens]
MTMHTTMTTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTMFMCLDQEVIISNWHWATTQTTQLSLSFKLDYFSMMFIPVALFVTWSIMEFSLWYMNSDPNINQFFKYLLIFLITMLILVTANNLFQLFIGWEGVGIMSFLLISWWYA

* Comments:
If the data was downloaded from NCBI, defline of the FASTA includes protein GI, accession and a description of the sequence. 

If the data was downloaded from non-NCBI sources, the defline contains only the internal ID. 

* Example
>lcl|intID10000000002 unnamed protein product



#-----------------------------------------------------------------------------
1.4.	MitoCOGs_annotation.txt

Contains the list of MitoCOGs and their function in a comma-separated format.

<MitoCOG_ID>,<annotation>

* Example:
MitoCOG0001,NADH dehydrogenase subunit 2



#-----------------------------------------------------------------------------
1.5.	nuclear_encoded_orthologs.txt

Contains the list of orthologous nuclear-encoded proteins (1,317) of MitoCOGs in a comma-separated format as described in 1.1.



#-----------------------------------------------------------------------------
1.6.	nuclear_proteomes.txt

Contains list of nuclear proteomes (57) in a comma-separated format as described in 1.2.



#-----------------------------------------------------------------------------
1.7.	nuclear_encoded_orthologs_sequences.fa

Contains protein sequence of orthologous nuclear-encoded proteins (1,317) of MitoCOGs (1,317) in FASTA format as described in 1.3.



#-----------------------------------------------------------------------------
1.8.	intron_gain_loss_estimation_tree.tre

This is the phylogenetic tree file in Newick format used for estimating intron gain and loss in the nuclear-encoded proteins. 



#-----------------------------------------------------------------------------
1.9.	concatenated_intron_presence_absence_matrix.txt

This is the concatenated intron presence absence matrix used for estimating intron gain and loss in the nuclear-encoded proteins. 
















