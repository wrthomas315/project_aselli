## *Sorex araneus* Pancreas of Aselli Paper
Scripts and code to reproduce transcriptomic comparison of the pancreas of Aselli and the spleen in *Sorex araneus* (common shrew), the spleen and bursa of Fabricius in *Gallus gallus* (chicken), overlap between the two comparisons, and genes that are also under positive selection in *S. araneus*.

### Goals and strategy

The first objective of this project is to determine if the pancreas of Aselli is a primary immune organ by comparing it to a secondary immune organ, the spleen, using RNA-sequencing data. The second objective is to test whether the pancreas of Aselli is an evolutionary analog of the cloacal bursa of Fabricius in birds, by conducting the same test between the bursa of Fabricius and spleen using chicken RNA-sequencing data found on the Sequencing Read Archives (SRA) and testing if overlap in differential expression is greater than expected by chance. Finally, we see if upregulated genes in the pancreas of Aselli are under positive selection and propose potential functional changes.


### Data, Code, and Analysis


#### Download transcriptomes for S. araneus and G. gallus. Download sequencing data from SRA for G. gallus bursa and spleen.

Note: Can skip this step if you want to begin with provided transcript abundances. Also, *S. araneus* data was local for us, so also would need to be downloaded from the SRA.

```
0_downloads.sh
```


#### Transcriptomics
Need to clean sequencing with fastp, pseudo-align to reference transcriptome with Kalisto
```
bash 1_fastp.sh
bash 2_kallisto.sh
```
Then run differential expression analyses, gene set enrichment, create PCAs, volcano plots, test overlap of two data sets, and overlap with positively selected genes.
Note: Done sequentially in Rstudio
```
R 3_poA_diffexp.R          
```

#### Immunoglobulin receptor diversity

Install ImReP according to source https://github.com/Mangul-Lab-USC/imrep

Then run on pancreas of Aselli reads
```
python2 imrep.py --fastq ML6_1.fq ML6_output > ML6output_std
python2 imrep.py --fastq ML7_1.fq ML7_output > ML7output_std
python2 imrep.py --fastq ML8_1.fq ML8_output > ML8output_std
python2 imrep.py --fastq ML9_1.fq ML9_output > ML9output_std
```

Sort outputs to look at diversity across all tissues
```
bash 4_IG_diversity.sh
```
