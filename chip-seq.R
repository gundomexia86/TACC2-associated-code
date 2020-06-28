library(ChIPseeker)
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

HDAC1_hepg2_2 <- readPeakFile('HDAC1-hepg2-2-ENCFF384XDG.bed.gz')
HDAC1_hepg2 <- readPeakFile('HDAC1-hepg2-ENCFF768XZW.bed.gz')
HDAC1_k562 <- readPeakFile('HDAC1-k562-ENCFF357NRL.bed.gz')
HDAC2_hepg2 <- readPeakFile('HDAC2-hepg2-ENCFF503AZK.bed.gz')
HDAC2_k562_2 <- readPeakFile('HDAC2-k562-2-ENCFF490VVG.bed.gz')
HDAC2_k562 <- readPeakFile('HDAC2-k562-ENCFF713HRG.bed.gz')
HMG20b_hepg2 <- readPeakFile('HMG20b-hepg2-ENCFF010AHT.bed.gz')
MTA1_k562 <- readPeakFile('MTA1-k562-ENCFF769GYG.bed.gz')
MTA2_k562_2 <- readPeakFile('MTA2-k562-2-ENCFF662MVX.bed.gz')
MTA2_k562 <- readPeakFile('MTA2-k562-ENCFF410EDG.bed.gz')
RCOR1_hepg2 <- readPeakFile('RCOR1-hepg2-ENCFF706RIA.bed.gz')

# 
peaks <- list(HDAC1_hepg2_2=HDAC1_hepg2_2,
              HDAC1_hepg2=HDAC1_hepg2,
              HDAC1_k562=HDAC1_k562,
              HDAC2_hepg2=HDAC2_hepg2,
              HDAC2_k562_2=HDAC2_k562_2,
              HDAC2_k562=HDAC2_k562,
              HMG20b_hepg2=HMG20b_hepg2,
              MTA1_k562=MTA1_k562,
              MTA2_k562_2=MTA2_k562_2,
              MTA2_k562=MTA2_k562,
              RCOR1_hepg2=RCOR1_hepg2)

#setting promotor region
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#gene ID annotation using annotatePeak
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE,annoDb="org.Hs.eg.db")

#visualization
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci relative to TSS")

tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color="red")


plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), 
            #conf=0.95,resample=500, 
            facet="row")


######gene_function_annotation
library(tidyr)
library(dplyr)
library(ggplot2)

df_HDAC1_hepg2_2 <- data.frame(peakAnnoList$HDAC1_hepg2_2)
df_HDAC1_hepg2 <- data.frame(peakAnnoList$HDAC1_hepg2)
df_HDAC1_k562 <- data.frame(peakAnnoList$HDAC1_k562)
df_HDAC2_hepg2 <- data.frame(peakAnnoList$HDAC2_hepg2)
df_HDAC2_k562_2 <- data.frame(peakAnnoList$HDAC2_k562_2)
df_HDAC2_k562 <- data.frame(peakAnnoList$HDAC2_k562)
df_HMG20b_hepg2 <- data.frame(peakAnnoList$HMG20b_hepg2)
df_MTA1_k562 <- data.frame(peakAnnoList$MTA1_k562)
df_MTA2_k562_2 <- data.frame(peakAnnoList$MTA2_k562_2)
df_MTA2_k562 <- data.frame(peakAnnoList$MTA2_k562)
df_RCOR1_hepg2 <- data.frame(peakAnnoList$RCOR1_hepg2)


write.csv(df_HDAC1_hepg2_2,file = 'peaks_annotation_HDAC1_hepg2_2.csv', row.names = F)
write.csv(df_HDAC1_hepg2,file = 'peaks_annotation_HDAC1_hepg2.csv', row.names = F)
write.csv(df_HDAC1_k562,file = 'peaks_annotation_HDAC1_k562.csv', row.names = F)
write.csv(df_HDAC2_hepg2,file = 'peaks_annotation_HDAC2_hepg2.csv', row.names = F)
write.csv(df_HDAC2_k562_2,file = 'peaks_annotation_HDAC2_k562_2.csv', row.names = F)
write.csv(df_HDAC2_k562,file = 'peaks_annotation_HDAC2_k562.csv', row.names = F)
write.csv(df_HMG20b_hepg2,file = 'peaks_annotation_HMG20b_hepg2.csv', row.names = F)
write.csv(df_MTA1_k562,file = 'peaks_annotation_MTA1_k562.csv', row.names = F)
write.csv(df_MTA2_k562_2,file = 'peaks_annotation_MTA2_k562_2.csv', row.names = F)
write.csv(df_MTA2_k562,file = 'peaks_annotation_MTA2_k562.csv', row.names = F)
write.csv(df_RCOR1_hepg2,file = 'peaks_annotation_RCOR1_hepg2.csv', row.names = F)

gene_HDAC1_hepg2_2 <- peakAnnoList$HDAC1_hepg2_2@anno$geneId
gene_HDAC1_hepg2 <- peakAnnoList$HDAC1_hepg2@anno$geneId
gene_HDAC1_k562 <- peakAnnoList$HDAC1_k562@anno$geneId
gene_HDAC2_hepg2 <- peakAnnoList$HDAC2_hepg2@anno$geneId
gene_HDAC2_k562_2 <- peakAnnoList$HDAC2_k562_2@anno$geneId
gene_HDAC2_k562 <- peakAnnoList$HDAC2_k562@anno$geneId
gene_HMG20b_hepg2 <- peakAnnoList$HMG20b_hepg2@anno$geneId
gene_MTA1_k562 <- peakAnnoList$MTA1_k562@anno$geneId
gene_MTA2_k562_2 <- peakAnnoList$MTA2_k562_2@anno$geneId
gene_MTA2_k562 <- peakAnnoList$MTA2_k562@anno$geneId
gene_RCOR1_hepg2 <- peakAnnoList$RCOR1_hepg2@anno$geneId


# Run GO enrichment analysis 
go <- enrichGO(gene = gene_RCOR1_hepg2, 
               OrgDb = org.Hs.eg.db, 
               ont = "all" )
barplot(go, split = 'ONTOLOGY') + facet_grid(ONTOLOGY~., scale="free")
go_df <- data.frame(go)
write.csv(go_df, file = 'go_RCOR1_hepg2.csv')


#kegg
kegg <- enrichKEGG(gene= gene_RCOR1_hepg2,
                   organism = 'hsa', pvalueCutoff = 0.05)
barplot(kegg)
#dotplot(kegg)
kegg_df <- as.data.frame(kegg)
write.csv(kegg_df, file = 'kegg_RCOR1_hepg2.csv')













