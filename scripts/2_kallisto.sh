#!/bin/bash

#Map chicken
while read j; do
	echo $j
	kallisto quant -i ~/project_aselli/data/0_refs/GCF_000002315.6_GRCg6a_rna.idx -o "~/project_aselli/data/5_TranscriptAbundances/"$j "~/project_aselli/data/4_trimmedRNA/$j"_1.trimmed.fastq" "~/project_aselli/data/4_trimmedRNA/"$j"_2.trimmed.fastq"
done <~/project_aselli/data/1_ids/bof4SRA.txt

#Map Sorex
while read j; do
	echo $j
	kallisto quant -i ~/project_aselli/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.idx -o "~/project_aselli/data/5_TranscriptAbundances/"$j "~/project_aselli/data/4_trimmedRNA/$j"_1.trimmed.fastq" "~/project_aselli/data/4_trimmedRNA/"$j"_2.trimmed.fastq"
done <~/project_aselli/data/1_ids/poa4fastp.txt
