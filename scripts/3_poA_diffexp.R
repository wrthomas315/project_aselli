### Goals
#1. Translating pancreas of Aselli kalisto transcript counts to gene counts (STEP1-3)
#2. Aselli input into DESeq2, normalize counts, diff exp, fsgea (STEP4)\
#3. Positive selection overlap (STEP5)
#4. Bursa input into DESeq2, normalize counts, diff exp (STEP6)
#5. Make plots for figure 3B (STEP7)
#6. Alignment and domain visualizations (STEP8)
#7. Overlap between two datasets, test for significance with Chi-squared tests (STEP9)




###Step1: First set up all your libraries
library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(tibble)
library( "genefilter" )
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(TCseq)
library(cluster)
library(EnhancedVolcano)
library(fgsea)
library(biomaRt)
library(ggtree)
library(Biostrings)
library(ggmsa)





### STEP2: Create a mechanism for getting transcript to gene
POAtxdb <- makeTxDbFromGFF("~/project_aselli/data/0_refs/GCF_027595985.1_mSorAra2.pri_genomic.gtf",format = "auto")
POA_k <- keys(POAtxdb, keytype="TXNAME")       
POAtx2gene <- AnnotationDbi::select(POAtxdb, POA_k, "GENEID", "TXNAME")
POAtx2gene <- POAtx2gene[!duplicated(POAtx2gene[,1]),]
POAtx2gene <- na.omit(POAtx2gene)





### STEP3: Import kallisto quantifications and write out
POA_h <- read.table("~/project_aselli/data/1_ids/ids.txt", header = T)
POA_files_h <-file.path("~/project_aselli/data/5_TranscriptAbundances", POA_h$Sample_name, "abundance.tsv")
names(POA_files_h) <- paste0("sample_", POA_h$Sample_name)
all(file.exists(POA_files_h))
POA_h.count.tsv <- tximport(POA_files_h, type = "kallisto", tx2gene = POAtx2gene, ignoreAfterBar=TRUE)
POA_h.tpm.tsv <- tximport(POA_files_h, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = POAtx2gene, ignoreAfterBar=TRUE)
#Write abundance outputs
write.table(POA_h.tpm.tsv$abundance, "~/project_aselli/analysis/DESeq/Aselli/poa_spl.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(POA_h.count.tsv$abundance, "~/project_aselli/analysis/DESeq/Aselli/poa_spl.count.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


#Normalize and set up DESeq
poa_id <- factor(c(POA_h$Sample_name))
poa_organs <- factor(c(POA_h$Organ))
poa_sex <- factor(c(POA_h$Sex))
poa_frame <-cbind(as.data.frame(poa_id),as.data.frame(poa_organs),as.data.frame(poa_sex))

#Make DESeq onject
colnames(POA_h.count.tsv$counts) <- c("poa_1","poa_2","poa_3","poa_4","spl_1","spl_2","spl_3","spl_4")
dds_poa <- DESeqDataSetFromMatrix(round(POA_h.count.tsv$counts), DataFrame(poa_frame), ~poa_sex + poa_organs)
mcols(dds_poa) <- cbind(mcols(dds_poa), row.names(POA_h.count.tsv$counts))
rownames(dds_poa) <- row.names(POA_h.count.tsv$counts)
dds_poa <- DESeq(dds_poa)
vst_dds_poa <- vst(dds_poa)
#Generate PCA
pca_dds_poa<- plotPCA(vst_dds_poa,intgroup=c("poa_organs"), ntop=5000, returnData=TRUE)
ggplot(pca_dds_poa, aes(x = PC1, y = PC2, color = factor(poa_organs))) +
  scale_color_manual(values = c("#009E73", "#E69F00"))+
  xlim(-60, 60) +
  ylim(-30, 25) +
  geom_point(size=10)+
  theme_bw()
#Look at distances if desired
poasampleDists <- dist(t(assay(vst_dds_poa)))
poasampleDistMatrix <- as.matrix(poasampleDists)
colnames(poasampleDistMatrix) <- NULL
##make the heatmap
pheatmap(poasampleDistMatrix, clustering_distance_rows=poasampleDists,
         clustering_distance_cols = poasampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
hyptopVarGenes <- head( order( rowVars( assay(vst_dds_poa) ), decreasing=TRUE ), 40 )
heatmap.2( assay(vst_dds_poa)[ hyptopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)






# STETP 4: DESEQ - No log fold changes as we are using a ranked gene set enrichment
lfc <- 0
poAres <- results(dds_poa, contrast = c("poa_organs","poA","Spleen"))
poAresSig <- subset(poAres,poAres$padj<.05)
poAresSigLog <- subset(poAresSig,abs(poAresSig$log2FoldChange)>=lfc)
poAup <- subset(poAresSigLog,(poAresSigLog$log2FoldChange)>=0)
poAdown <- subset(poAresSigLog,(poAresSigLog$log2FoldChange)<=0)
length(poAup$log2FoldChange)
length(poAdown$log2FoldChange)
#Write outputs for
#Total results
write.table(poAres, "~/project_aselli/analysis/DESeq/Aselli/poa_DESeq.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#Significant result
write.table(poAresSig, "~/project_aselli/analysis/DESeq/Aselli/poa_Sig.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#Sig upregulated
write.table(poAup, "~/project_aselli/analysis/DESeq/Aselli/poaUpsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#Sig downregulated
write.table(poAdown, "~/project_aselli/analysis/DESeq/Aselli/poaDownsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

####CREATE VOLCANO PLOT FOR FIGURE
poAresX <-  poAres
for (i in 1:length(poAresX$padj)) {
  if  (poAresX$padj[i]<1e-40 & !is.na (poAresX$padj[i])) {
    poAresX$padj[i] <- 1e-40
  }
  if (poAresX$log2FoldChange[i]>10 & !is.na (poAresX$log2FoldChange[i])) {
    poAresX$log2FoldChange[i] <- 10
  }
  if (poAresX$log2FoldChange[i]< -10 & !is.na (poAresX$log2FoldChange[i])) {
    poAresX$log2FoldChange[i] <- -10
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(poAresX$log2FoldChange) == 10, 17,
  ifelse(poAresX$padj==1e-40, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- '<log1e-15'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  poAresX$padj > 0.05, 'grey',
  ifelse(poAresX$log2FoldChange <= -1.58, 'red',
         ifelse(poAresX$log2FoldChange >= 1.58, 'blue',
                ifelse(poAresX$log2FoldChange >= 0, 'blue',
                       'red'))))
keyvals
keyvals[is.na(keyvals)] <- 'grey'
str(keyvals)
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
length(keyvals)
length(keyvals.shape)
#
EnhancedVolcano(poAresX,
                lab = rownames(poAresX),
                xlim=c(-11 ,11),
                ylim=c(0,40),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = c("PRDM1", "XBP1", "IRF4", "PAX5", "BACH2", "BCL6"),
                pCutoff = .05,
                FCcutoff = 15,
                colCustom = keyvals,
                colAlpha = .5,
                boxedLabels = TRUE,
                pointSize = 4.0,
                legendPosition = 'none',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.5)

#fsgea
poA_DESeq <- read_delim("~/project_aselli/analysis/DESeq/poa_DESeq.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)
poA_hres <- poA_DESeq %>% 
  dplyr::select(Gene, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(stat))
library(tidyverse)
poa_hranks <- deframe(poA_hres)
fgsea_poA_DESeq <- fgsea(pathways=gmtPathways("~/project_aselli/data/0_refs/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), poa_hranks) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_poA_DESeqTidy <- fgsea_poA_DESeq %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_poA_DESeqTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
fgsea_poA_DESeqTidy
poA_fsgea <-ggplot(subset(fgsea_poA_DESeqTidy,padj<0.05 & (NES< -1.8 | NES>0)), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
poA_fsgea
fgsea_poA_DESeqTidy2 <- apply(fgsea_poA_DESeqTidy,2,as.character)
write.table(fgsea_poA_DESeqTidy2, file='~/project_aselli/analysis/DESeq/fsgea.txt', quote=FALSE, sep='\t')

#Quick plots for...
#Maintain B cells The genetic network controlling plasma cell differentiation
plotCounts(dds_poa, gene="PAX5", intgroup="poa_organs") #nope!
plotCounts(dds_poa, gene="BACH2", intgroup="poa_organs") #yup
plotCounts(dds_poa, gene="BCL6", intgroup="poa_organs") #meh!
#facilitate plasma cell differentiation
plotCounts(dds_poa, gene="IRF4", intgroup="poa_organs")
plotCounts(dds_poa, gene="XBP1", intgroup="poa_organs")
plotCounts(dds_poa, gene="PRDM1", intgroup="poa_organs")
#differentiation is key for pancreas of Aselli






### STEP 5: Overlap with positively selected genes
shrewtotalPSG <- read_table("~/project_aselli/data/0_refs/shrewtotal.txt")
shrewspecificPSG <-read_table("~/project_aselli/data/0_refs/shrewspecific.txt")

#Not shrew specific
nss_tot <- Reduce(intersect, list(shrewtotalPSG$p_adjusted,rownames(poAup))) #64 total
write.table(nss_tot, "~/project_aselli/analysis/DESeq/overlap/nss_tot.txt", na = "NA", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
sas_tot <-Reduce(intersect, list(shrewspecificPSG$p_adjusted,rownames(poAup))) #21 total
write.table(sas_tot, "~/project_aselli/analysis/DESeq/overlap/sas_tot.txt", na = "NA", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)




#STEP 6: bursa of Fabricius
#gallus gallus
galtxdb <- makeTxDbFromGFF("~/project_aselli/data/refs/GCF_000002315.5_GRCg6a_genomic.gtf",format = "auto")
gal_k <- keys(galtxdb, keytype="TXNAME")       
galtx2gene <- AnnotationDbi::select(galtxdb, gal_k, "GENEID", "TXNAME")
galtx2gene <- galtx2gene[!duplicated(galtx2gene[,1]),]
galtx2gene <- na.omit(galtx2gene)

### Import kallisto quantifications and write out
gal_h <- read.table("~/project_aselli/data/ids/ids_bof.txt", header = T)
gal_files_h <-file.path("~/project_aselli/data/5_TranscriptAbundances", gal_h$ID, "abundance.tsv")
names(gal_files_h) <- paste0("sample_", gal_h$ID)
all(file.exists(gal_files_h))
gal_h.count.tsv <- tximport(gal_files_h, type = "kallisto", tx2gene = galtx2gene, ignoreAfterBar=TRUE)
gal_h.tpm.tsv <- tximport(gal_files_h, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = galtx2gene, ignoreAfterBar=TRUE)
write.table(gal_h.tpm.tsv$abundance, "~/project_aselli/analysis/DESeq/Fabricius/gal_spl.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(gal_h.count.tsv$abundance, "~/project_aselli/analysis/DESeq/Fabricius/gal_spl.count.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

#Normalize and set up DESeq
gal_id <- factor(c(gal_h$ID))
gal_organs <- factor(c(gal_h$Organ))
gal_sex <- factor(c(gal_h$Sex))
gal_breed <- factor(c(gal_h$Breed))
gal_frame <-cbind(as.data.frame(gal_id),as.data.frame(gal_organs),as.data.frame(gal_sex),as.data.frame(gal_breed))

#Make DESeq onject
colnames(gal_h.count.tsv$counts) <- c("bof_1","bof_2","bof_3","bof_4","bof_5","spl_1","spl_2","spl_3")
dds_gal <- DESeqDataSetFromMatrix(round(gal_h.count.tsv$counts), DataFrame(gal_frame), ~gal_breed + gal_organs)
mcols(dds_gal) <- cbind(mcols(dds_gal), row.names(gal_h.count.tsv$counts))
rownames(dds_gal) <- row.names(gal_h.count.tsv$counts)
dds_gal <- DESeq(dds_gal)
vst_dds_gal <- vst(dds_gal)
#Generate PCA
pca_dds_gal<- plotPCA(vst_dds_gal,intgroup=c("gal_organs"), ntop=1000, returnData=TRUE)
ggplot(pca_dds_gal, aes(x = PC1, y = PC2, color = factor(gal_organs))) +
  geom_point(size=2)+
  theme_bw()
#Look at distances if desired
galsampleDists <- dist(t(assay(vst_dds_gal)))
galsampleDistMatrix <- as.matrix(galsampleDists)
colnames(galsampleDistMatrix) <- NULL
##make the heatmap
pheatmap(galsampleDistMatrix, clustering_distance_rows=galsampleDists,
         clustering_distance_cols = galsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
hyptopVarGenes <- head( order( rowVars( assay(vst_dds_gal) ), decreasing=TRUE ), 40 )
heatmap.2( assay(vst_dds_gal)[ hyptopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

subset(galres, rownames(galres)=="BACH2")
subset(poAres, rownames(poAres)=="BACH2")
#DESEQ - No LFC to keep the same as used in pancreas of Aselli
lfc <- 0
galres <- results(dds_gal, contrast = c("gal_organs","BoF","Spleen"))
galresSig <- subset(galres,galres$padj<.05)
galresSigLog <- subset(galresSig,abs(galresSig$log2FoldChange)>=lfc)
galup <- subset(galresSigLog,(galresSigLog$log2FoldChange)>=0)
galdown <- subset(galresSigLog,(galresSigLog$log2FoldChange)<=0)
#determine how many significant resultes
length(galup$log2FoldChange)
length(galdown$log2FoldChange)
#look at overlap between data sets
Reduce(intersect, list(rownames(galup),rownames(poAup)))
Reduce(intersect, list(rownames(galdown),rownames(poAup)))
Reduce(intersect, list(rownames(galdown),rownames(poAdown)))
Reduce(intersect, list(rownames(galup),rownames(poAdown)))
#visualize overlap with an upset plot (For figure 3A)
poAbursagenelist <- list(poAup = rownames(poAup),poAdown = rownames(poAdown),galup = rownames(galup),galdown = rownames(galdown))
upset_poAbursagenelist <- fromList(poAbursagenelist)
upset(upset_poAbursagenelist,nsets = 4,order.by = "freq",)

#Visualize key genes and compare to pancreas of Aselli
#differentiation into plastma cells is poA
plotCounts(dds_gal, gene="IRF4", intgroup="gal_organs")
plotCounts(dds_poa, gene="IRF4", intgroup="poa_organs")
plotCounts(dds_gal, gene="XBP1", intgroup="gal_organs")
plotCounts(dds_poa, gene="XBP1", intgroup="poa_organs")
plotCounts(dds_gal, gene="PRDM1", intgroup="gal_organs")
plotCounts(dds_poa, gene="PRDM1", intgroup="poa_organs")
#maintaining b cell is bursa like bone marrow
plotCounts(dds_gal, gene="PAX5", intgroup="gal_organs")
plotCounts(dds_poa, gene="PAX5", intgroup="poa_organs")
plotCounts(dds_gal, gene="BACH2", intgroup="gal_organs")
plotCounts(dds_poa, gene="BACH2", intgroup="poa_organs")
plotCounts(dds_gal, gene="BCL6", intgroup="gal_organs")
plotCounts(dds_poa, gene="BCL6", intgroup="poa_organs")
#Write ouputs
write.table(galres, "~/project_aselli/analysis/DESeq/Fabricius/gal_DESeq.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(galresSig, "~/project_aselli/analysis/DESeq/Fabricius/gal_Sig.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(galup, "~/project_aselli/analysis/DESeq/Fabricius/galUpsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(galdown, "~/project_aselli/analysis/DESeq/Fabricius/galDownsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)




#STEP 7 - Make figures for Figure 3B, run for each gene here shown PAX5
# Extracting PAX5 rows from both dataframes
POA_PAX5 <- subset(POA_h.tpm.tsv$abundance, rownames(POA_h.tpm.tsv$abundance) == "PAX5")
gal_PAX5 <- subset(gal_h.tpm.tsv$abundance, rownames(gal_h.tpm.tsv$abundance) == "PAX5")

# Transpose data to make samples into rows
POA_PAX5_t <- t(POA_PAX5)
gal_PAX5_t <- t(gal_PAX5)

# Convert to data frames and add a sample identifier
POA_df <- data.frame(Sample = rownames(POA_PAX5_t), Count = as.numeric(POA_PAX5_t[,1]), stringsAsFactors = FALSE)
gal_df <- data.frame(Sample = rownames(gal_PAX5_t), Count = as.numeric(gal_PAX5_t[,1]), stringsAsFactors = FALSE)

# Combine the two datasets
combined_df <- rbind(POA_df, gal_df)
combined_df$Type <- c("Aselli", "Aselli", "Aselli", "Aselli",
                      "Spleen_s", "Spleen_s", "Spleen_s", "Spleen_s",
                      "BoF", "BoF", "BoF", "BoF", "BoF",
                      "Spleen_g", "Spleen_g", "Spleen_g")
summary_df <- combined_df %>%
  group_by(Type) %>%
  summarize(
    Mean = mean(Count),
    SD = sd(Count),
    .groups = 'drop'
  )
# Add a 'Group' variable to control grouping in the plot
summary_df <- summary_df %>%
  mutate(Group = ifelse(Type %in% c("Aselli", "Spleen_s"), "Group1", "Group2"))

# Reorder 'Type' factor for correct order within groups
summary_df$Type <- factor(summary_df$Type, levels = c("Aselli", "Spleen_s", "BoF", "Spleen_g"))

# Plot with grouped bars and gap between the groups
ggplot(summary_df, aes(x = Group, y = Mean, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge2(.9), color = "black") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                position = position_dodge(width = .9), width = 0.2, color = "black") +
  scale_fill_manual(values = c("Aselli" = "#009E73", "Spleen_s" = "#E69F00", 
                               "BoF" = "#D81B60", "Spleen_g" = "#E69F00")) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0, limits =c(0,20))) + 
  theme(legend.position = "none",           
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),      
        axis.ticks = element_blank()
  )
####




#STEP 8 Alignment and domain visualizations
#PRDM1
#amino acid boces
boxes_4R <- read_delim("~/ShrewProjects/project_aselli/analysis/PRDM1/boxes_4R.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
#Now this is where things get even more janky
#the bar plots will NOT recognize your species and put them aligned with it
#but will do so based on numerical order in which they appear
#set seq from 1 to how many species you have
#then order of species from bottom branch to top branch
#can be found by continuous printing phylogeny with labeled branch (line 8), ir printing your phylogeny with names but ggtree hates printing words
spec_key <- data.frame(
  A = seq(1, 40, by = 1),
  B = c("microcebus_murinus","macaca_mulatta","papio_anubis","hylobates_moloch","homo_sapiens","pan_troglodytes","tupaia_chinensis","oryctolagus_cuniculus","sciurus_vulgaris","peromyscus_leucopus","acomys_russatus","mus_musculus","rattus_norvegicus","solenodon_paradoxus","erinaceus_europaeus","suncus_etruscus","sorex_araneus","condylura_cristata","talpa_occidentalis","scalopus_aquaticus","rhinolophus_ferrumequinum","phyllostomus_discolor","pipistrellus_pipistrellus","myotis_myotis","manis_javanica","ceratotherium_simum_cottoni","camelus_ferus","sus_scrofa","phocoena_sinus","cervus_elaphus","ovis_orientalis","bos_gaurus","lynx_canadensis","vulpes_lagopus","ailuropoda_melanoleuca","martes_zibellina","lutra_lutra","mustela_erminea","mustela_nigripes","mustela_putorius_furo"))
merged_df <- merge(boxes_4R, spec_key, by.x = "species", by.y = "B", all.x = TRUE)
merged_df
order_ENST1 <- as.data.frame(cbind(merged_df$species,merged_df$start_position,merged_df$domain,merged_df$end_position,merged_df$species,merged_df$A))
order_ENST1
order_ENST1$V2<-as.numeric(order_ENST1$V2)
order_ENST1$V3<-as.character(order_ENST1$V3)
order_ENST1$V4<-as.numeric(order_ENST1$V4)
order_ENST1$V5<-as.character(order_ENST1$V5)
order_ENST1$V6<-as.numeric(order_ENST1$V6)
# Plot using ggplot2 with just ticks, see how it looks without bars
PRDM1snp_data<-read_table("~/project_aselli/analysis/PRDM1/PRDM1_rfigure_input.txt", 
                          col_names = FALSE)
PRDM1_sor <- subset(order_ENST1, V1 == "sorex_araneus")
PRDM1_sor$size <- ifelse(PRDM1_sor$V3 == "NO", 15, 20)  # "NO" segments thinner (10), others thicker (20)
PRDM1_sor$V4 <- PRDM1_sor$V4 +1
####
ggplot(PRDM1_sor, aes(x = V2, xend = V4, y = 1, yend = 1, color = V3,, size = size)) +
  # Segment bar
  geom_segment() +
  # Mutation lollipops with adjusted position
  geom_segment(
    data = PRDM1snp_data,
    aes(x = X2, xend = X2, y = 1 + 10 / 2, yend = 1 + 10 / 2 + 0.5),  # Stick starts at top of the bars
    color = "black",  # Stick color
    inherit.aes = FALSE
  ) +
  geom_point(
    data = PRDM1snp_data,
    aes(x = X2, y = 1 + 10 / 2 + 0.5),  # Head just above the stick
    color = "red",  # Head color
    size = 3,
    inherit.aes = FALSE
  ) +
  # Custom colors for segments
  scale_color_manual(values = c("NO" = "grey", "SET" = "#FFC20A", "ZF" = "#0C7BDC")) +
  scale_size_identity() +  # Use the size column directly
  theme_minimal() +
  labs(
    title = "Segment Bar with Mutations for Sorex araneus",
    x = "Position",
    y = NULL,  # Remove y-axis label
    color = "Segment Label"
  ) +
  theme(
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove grid lines
    legend.position = "none"
  )
#PRDM1
###Graphics for genes
ENST00000648754_tree <- read.tree("~/project_aselli/analysis/PRDM1/noFSENST00000648754.macse.aligned.fa.taper.fa.nh")
#phylopic depreciated but still run to create needed frame for fore
ENST00000648754_tree_images <- ggimage::phylopic_uid(ENST00000648754_tree$tip.label)
ENST00000648754_tree_images$labeling <-  ENST00000648754_tree_images$name
ENST00000648754_tree_images$labeling
ENST00000648754_tree_images$fore <- c(rep("Background",32),rep("Foreground",1),rep("Background",5))
#check to make sure foreground and background are w/ correct species
ENST00000648754_tree_images
head(ENST00000648754_tree_images)
ENST00000648754_phy <-ggtree(ENST00000648754_tree,branch.length="none") %<+% ENST00000648754_tree_images + 
  #geom_tiplab(aes(label=labeling,color=fore), offset=-6, size =4)+
  scale_color_manual(values=c("black","lightslateblue"))
ENST00000648754_phy
ENST00000648754strings_aa <- readAAStringSet("~/project_aselli/analysis/PRDM1/seqs/AA/noFSENST00000648754.macse.aligned.taper.aa")
#Choose which sites you want to look at
ENST00000648754strings_aa_dat_1 <- tidy_msa(ENST00000648754strings_aa,290,300)
ENST00000648754strings_aa_dat_2 <- tidy_msa(ENST00000648754strings_aa,475,485)
ENST00000648754strings_aa_dat_3 <- tidy_msa(ENST00000648754strings_aa,566,576)
facet_widths(ENST00000648754_phy + geom_facet(geom = geom_msa, data = ENST00000648754strings_aa_dat_1,  panel = 'msa1')+geom_facet(geom = geom_msa, data = ENST00000648754strings_aa_dat_2,  panel = 'msa2')+geom_facet(geom = geom_msa, data = ENST00000648754strings_aa_dat_3,  panel = 'msa3')+theme_tree2(legend.position = 'none'),c(.10,.2,.2,.2))


###PTPRCAP
####PTPRCAP Boxes
#amino acid boces
PTPRCAPboxes_4R <- read_delim("~/project_aselli/analysis/PTPRCAP/PTPRCAP_boxes_4R.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

PTPRCAPsnp_data<-read_table("~/project_aselli/analysis/PTPRCAP/PTPRCAP_4R_ticks.txt", 
                          col_names = TRUE)

PTPRCAPboxes_4R$size <- ifelse(PTPRCAPboxes_4R$domain %in% c("NO", "I", "D"), 15, 20)####
PTPRCAPboxes_4R$end_position <- PTPRCAPboxes_4R$end_position +1
# Plot with mutations as lollipops (corrected y positioning)
PTPRCAPsnp_data
ggplot(PTPRCAPboxes_4R, aes(x = start_position, xend = end_position, y = 1, yend = 1, color = domain, size = size)) +
  # Segment bar
  geom_segment() +
  # Mutation lollipops with adjusted position
  geom_segment(
    data = PTPRCAPsnp_data,
    aes(x = boxes, xend = boxes, y = 1 + 10 / 2, yend = 1 + 10 / 2 + 0.5),  # Stick starts at top of the bars
    color = "black",  # Stick color
    inherit.aes = FALSE
  ) +
  geom_point(
    data = PTPRCAPsnp_data,
    aes(x = boxes, y = 1 + 10 / 2 + 0.5, color = `5indel`),  # Color by 5indel category
    size = 3,
    inherit.aes = FALSE
  ) +
  # Define colors for both segment domains and mutations
  scale_color_manual(values = c("NOT" = "red", "NO"="grey", "YES" = "#56B4E9", "PTPRCAP" = "#CC79A7", "EXTRA" = "#F0E442", "TRANS"="#009E73")) +
  scale_size_identity() +  # Use the size column directly
  theme_minimal() +
  labs(
    title = "Segment Bar with Mutations for Sorex araneus",
    x = "Position",
    y = NULL,  # Remove y-axis label
    color = "Segment Label"
  ) +
  theme(
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove grid lines
    legend.position = "none"
  )
#PTPRCAP
###Graphics for genes
ENST00000326294_tree <- read.tree("~/project_aselli/analysis/PTPRCAP/tree/noFSENST00000326294.macse.aligned.fa.taper.fa.nh")
#phylopic depreciated but still run to create needed frame for fore
ENST00000326294_tree_images <- ggimage::phylopic_uid(ENST00000326294_tree$tip.label)
ENST00000326294_tree_images$labeling <-  ENST00000326294_tree_images$name
ENST00000326294_tree_images$labeling
ENST00000326294_tree_images$fore <- c(rep("Background",8),rep("Foreground",1),rep("Background",24),rep("Foreground",1),rep("Background",5))
#check to make sure foreground and background are w/ correct species
ENST00000326294_tree_images
ENST00000326294_phy <-ggtree(ENST00000326294_tree,branch.length="none") %<+% ENST00000326294_tree_images + 
  #geom_tiplab(aes(label=labeling,color=fore), offset=0.04, size =3)+
  scale_color_manual(values=c("black","lightslateblue"))
ENST00000326294_phy
#read in MSA as strings
ENST00000326294strings_nt <- readAAStringSet("~/project_aselli/analysis/PTPRCAP/NT/noFSENST00000326294.macse.aligned.fa.taper.fa")
ENST00000326294strings_aa <- readAAStringSet("~/project_aselli/analysis/PTPRCAP/AA/noFSENST00000326294.macse.aligned.fa.taper.aa")
#set values for amino acids or nucleotides you would want to see
ENST00000326294strings_aa_dat_1 <- tidy_msa(ENST00000326294strings_aa,22,58)
#ENST00000326294strings_aa_dat_2 <- tidy_msa(ENST00000326294strings_aa,188,208), could look at second
facet_widths(ENST00000326294_phy + geom_facet(geom = geom_msa, data = ENST00000326294strings_aa_dat_1,  panel = 'msa')+theme_tree2(legend.position = 'none'),c(.1,.6))




#STEP 9 : Testing if overlap between the differentially expressed genes are significant
#### First look at correlation
# Convert DESeqResults to data frames
poAres_df <- as.data.frame(poAresSig)
galres_df <- as.data.frame(galresSig)

# Extract gene names
poAres_df$gene <- rownames(poAres_df)
galres_df$gene <- rownames(galres_df)

# Merge data frames by gene
merged_df <- merge(poAres_df[, c("gene", "log2FoldChange")], 
                   galres_df[, c("gene", "log2FoldChange")], 
                   by = "gene", suffixes = c("_poA", "_gal"))

# Remove NA and infinite values
merged_df <- na.omit(merged_df)
merged_df <- merged_df[is.finite(merged_df$log2FoldChange_poA) & is.finite(merged_df$log2FoldChange_gal), ]

# Check the number of overlapping genes
print(paste("Number of overlapping genes:", nrow(merged_df)))

# Compute Pearson correlation
if (nrow(merged_df) > 1) {  # Ensure enough data points
  correlation <- cor(merged_df$log2FoldChange_poA, merged_df$log2FoldChange_gal, method = "pearson")
} else {
  correlation <- NA
}

# Print correlation
print(paste("Pearson correlation:", round(correlation, 3)))

# Plot only if correlation is defined
if (!is.na(correlation)) {
  library(ggplot2)
  ggplot(merged_df, aes(x = log2FoldChange_poA, y = log2FoldChange_gal)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(x = "log2FoldChange (Pancreas vs Spleen)", 
         y = "log2FoldChange (Bursa vs Spleen)", 
         title = paste("Correlation: ", round(correlation, 3))) +
    theme_minimal()
}

#Then use a hypogeometric test for total overlap+same direction greater than expected by chance
total_genes <- 12367
n_sim <- 100000

# Make a vector of gene IDs
all_genes <- paste0("gene", 1:total_genes)

# Number of up/down regulated genes
pancreas_up <- 1542
pancreas_down <- 2666
bursa_up <- 4254
bursa_down <- 3670

# Observed same-direction overlap
observed_overlap <- 543 + 964

# Simulation: sample independently for each tissue
simulated_overlaps <- replicate(n_sim, {
  # Randomly sample genes for each group
  pancreas_up_genes <- sample(all_genes, pancreas_up)
  pancreas_down_genes <- sample(setdiff(all_genes, pancreas_up_genes), pancreas_down)
  
  bursa_up_genes <- sample(all_genes, bursa_up)
  bursa_down_genes <- sample(setdiff(all_genes, bursa_up_genes), bursa_down)
  
  # Count overlap in same directions
  overlap_up_up <- length(intersect(pancreas_up_genes, bursa_up_genes))
  overlap_down_down <- length(intersect(pancreas_down_genes, bursa_down_genes))
  
  overlap_up_up + overlap_down_down
})

# P-value is how often simulated overlap ≥ observed
p_value <- mean(simulated_overlaps >= observed_overlap)


###Chi squared for overlap with PSG and upregulated
#shrew poa upregulated
rownames(poAup)
#shrew specific psgs
shrewspecificPSG$p_adjusted
#total diff exp
rownames(poAres)
#total PSG gene list
abs_for <- read_table("~/project_aselli/data/0_refs/output.txt", 
                      col_names = FALSE)
#import transcript to gene and add to column
LT_output_4hp <- read_table("~/project_aselli/data/0_refs/LT_output_4hp.txt", 
                            col_names = FALSE)
abs_for$X5 <- NA
for (i in seq_len(nrow(abs_for))) {
  # Find the matching row
  match_row <- match(abs_for[i, "X1"], LT_output_4hp$X3)
  
  # If a match found, assign the corresponding value to X5
  if (!is.na(match_row)) {
    abs_for[i, "X5"] <- LT_output_4hp[match_row, "X1"]
  }
}
#use gene ID to get ensemble gene name
library(dbplyr)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
abs_for_vec <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = abs_for$X5, 
                     mart = ensembl)
abs_for$X6 <- NA
# Iterate over rows in absrel results
for (i in seq_len(nrow(abs_for))) {
  # Find the matching row in ensembl data frame
  match_row <- match(abs_for[i, "X5"], abs_for_vec$ensembl_gene_id)
  
  # If a match found, assign the corresponding value to X6
  if (!is.na(match_row)) {
    abs_for[i, "X6"] <- abs_for_vec[match_row, "external_gene_name"]
  }
}
ara.sub <-subset(abs_for, !is.na(X3) & X2 == "sorex_araneus")
#Get overlaps
significant_overlap_PSG <- length(Reduce(intersect, list(shrewspecificPSG$p_adjusted,rownames(poAup))))
shrew_poAup <- length(intersect(rownames(poAup), ara.sub$X6))
shrewPSG <- length(intersect(shrewspecificPSG$p_adjusted, rownames(poAres)))
overlap_tot <- length(intersect(ara.sub$X6, rownames(poAres)))

# Calculate expected overlap under independence assumption
expected_overlap <- (shrew_poAup * shrewPSG) / overlap_tot

contingency_PSG <- matrix(c(
  significant_overlap_PSG, shrew_poAup-significant_overlap_PSG,
  shrewPSG-significant_overlap_PSG, overlap_tot - (shrew_poAup + shrewPSG - significant_overlap_PSG)),
  nrow = 2, byrow = TRUE)

# Perform Chi-squared test
chi_PSG <- chisq.test(contingency_PSG, correct = FALSE)
#X-squared = 1.1003, df = 1, p-value = 0.2942