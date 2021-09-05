library(Seurat)
library(rlang)
library(magrittr)
library(dplyr)
library(MASS)
library(GenomicRanges)
library(SnapATAC)
library(ggplot2)
library(viridis)
library(ggpubr)

# Load RNA matrix 
nat_com <- read.table("/home/donggun/hkid_sc_RNA_publish/nature_com_GSE121862/GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv")

# Subset for each cell type in scRAC-seq matrix & creation of pseudo-bulk matrix
nat_com_C3 <- grep('C2',names(nat_com))
Podo_matrix <- nat_com[,nat_com_C3]
exprMat <- rowSums(Podo_matrix)
idx <- RNA_ATAC_gene$gene %in% c("NPHS2","CLIC5","PLA2R1", "NPHS2") # marker genes 
exprMat <- log10(exprMat+1)
exprMat_data <- as.data.frame(exprMat)
exprMat_data$gene <- row.names(exprMat_data)

# Creat overlaped gene matrix in snATAC-seq snap object
genes = read.table("./gencode.v33.annotation.bed")
genes.gr = GRanges(genes[, 1], 
                   IRanges(genes[, 2], genes[, 3]), name = genes[, 5])
marker.genes <- row.names(nat_com)
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)]
x.sp = createGmatFromMat(
  obj = x.sp, 
  input.mat = "bmat",
  genes = genes.sel.gr,
  do.par = TRUE,
  num.cores = 10
)
x.sp = scaleCountMatrix(
  obj = x.sp, 
  cov = x.sp@metaData$passed_filters + 1,
  mat = "gmat",
  method = "RPM"
)
x.sp = runMagic(
  obj = x.sp,
  input.mat = "gmat",
  step.size = 3
)

# Subset for each cell type in snap object & creation of pseudo-bulk matrix 
ATAC_celltype <- subset(x.sp, x.sp@cluster %in% c("Podo"))
ATAC_gmat<- as.matrix(ATAC_celltype@gmat)
ATAC_gmat <- t(ATAC_gmat)
ATAC_gmat[1:5,1:5]
gene_name <- row.names(ATAC_gmat)
ATAC_gmat <- rowSums(ATAC_gmat)
ATAC_gmat <- log10(ATAC_gmat+1)
ATAC_matrix <- ATAC_gmat
ATAC_matrix[1:5]
ATAC_matrix <- as.data.frame(ATAC_matrix)
ATAC_matrix$gene <- gene_name

# Caculation of Pearson Correlation Coefficient using scRNA-seq, snATAC-seq pseudo-bulk matrix
RNA_ATAC_gene <- merge(x = exprMat_data, y = ATAC_matrix, by= "gene")
colnames(RNA_ATAC_gene)<- c("gene","RNA", "ATAC")
ggscatter(RNA_ATAC_gene, x = "ATAC", y = "RNA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ATAC", ylab = "RNA")
res <- cor.test(RNA_ATAC_gene$RNA, RNA_ATAC_gene$ATAC, method = "pearson")
res
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
RNA_ATAC_gene$Density <- get_density(RNA_ATAC_gene$ATAC, RNA_ATAC_gene$RNA, n = 100)
RNA_ATAC_gene$gene_1 <- RNA_ATAC_gene$gene
RNA_ATAC_gene$gene_1[!idx] <- NA
ggplot(RNA_ATAC_gene, aes(x = ATAC, y = RNA)) + geom_point(aes(ATAC, RNA, color = Density)) + scale_color_viridis()+  xlab("log10(scATAC-seq + 1)") + ylab("log10(scRNA-seq + 1)") + geom_label_repel(aes(label = gene_1), box.padding   = 0.35, point.padding = 0.5, segment.color = 'Black') 
ggsave("./nature_com_GSE121862/correlation/correlation_plot/Podo_correlation.pdf", width = 7, height = 6, units = "in")
