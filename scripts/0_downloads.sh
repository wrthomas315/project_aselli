#!/bin/bash

#References
#Get S. araneus transcriptome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_rna.fna.gz ~/project_aselli/data/0_refs
#index it
gunzip ~/project_aselli/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.fna.gz
kallisto index -i ~/project_aselli/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.idx ~/project_aselli/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.fna

#Get G. gallus transcriptome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.5_GRCg6a/GCF_000002315.5_GRCg6a_rna.fna.gz ~/project_aselli/data/0_refs
#index it
gunzip ~/project_aselli/data/0_refs/GCF_000002315.5_GRCg6a_rna.fna.gz
kallisto index -i ~/project_aselli/data/0_refs/GCF_000002315.5_GRCg6a_rna.idx ~/project_aselli/data/0_refs/GCF_000002315.5_GRCg6a_rna.fna

#Chicken sequencing data
#get useable SRA ids
awk '{print $1}' ~/project_aselli/data/1_ids/ids_bof.txt | grep "SRR" > ~/project_aselli/data/1_ids/bof4SRA.txt
# get sequencing
while read j; do
	echo $j
	prefetch -v $j
	fastq-dump -I --split-files --outdir "~/project_aselli/data/3_RNAseq/" $j".sra"
done <~/project_aselli/data/1_ids/bof4SRA.txt
