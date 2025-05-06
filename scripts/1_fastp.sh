#!/bin/bash

#Paired for chicken
while read j; do
        echo $j
        fastp --in1 "~/project_aselli/data/3_RNAseq/"$j"_1.fastq" --in2 "~/project_aselli/data/3_RNAseq/"$j"_2.fastq" --out1 "~/project_aselli/data/4_trimmedRNA/"$j"_1.trimmed.fastq" --out2 "~/project_aselli/data/4_trimmedRNA/"$j"_2.trimmed.fastq" --cut_tail --html $j".html" --json $j".json" 2> $j".log"
done <~/project_aselli/data/1_ids/bof4SRA.txt

#Paired for Sorex
awk '{print $1}' ~/ShrewProjects/github/project_aselli/data/1_ids/ids.txt | grep -v "Sample" >poa4fastp.txt
while read j; do
        echo $j
        fastp --in1 "~/project_aselli/data/3_RNAseq/"$j"_1.fastq" --in2 "~/project_aselli/data/3_RNAseq/"$j"_2.fastq" --out1 "~/project_aselli/data/4_trimmedRNA/"$j"_1.trimmed.fastq" --out2 "~/project_aselli/data/4_trimmedRNA/"$j"_2.trimmed.fastq" --cut_tail --html $j".html" --json $j".json" 2> $j".log"
done <~/project_aselli/data/1_ids/poa4fastp.txt
#Note, would need to switch to SRA names if downloading from SRA
