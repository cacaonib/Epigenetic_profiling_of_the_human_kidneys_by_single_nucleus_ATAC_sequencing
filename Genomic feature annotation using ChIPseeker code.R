library(ChIPseeker)
library(ensembldb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(UpSetR)
library(ggupset)
# Load DAR bed files  
files <- c("./Endo_DAR.bed.gz","./Podo_DAR.bed.gz","./PT_DAR.bed.gz","./LOH_DAR.bed.gz","./DCT_DAR.bed.gz","./CNT_DAR.bed.gz",
           "./CD-PC_DAR.bed.gz","./CD-a-IC_DAR.bed.gz","./CD-b-IC_DAR.bed.gz","./Fibro_DAR.bed.gz","./Macro_DAR.bed.gz", "./Lympho_DAR.bed.gz")
names(files) = levels(x.sp@cluster)
print(files)
peak <- readPeakFile(files[[4]])
peakAnno <- annotatePeak(files[[3]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
upsetplot(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
# KEGG pathway 
library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)
gene_df <- as.data.frame(gene)
# Visualization all cell types 
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
peakAnnoList
