library(SnapATAC)
library(GenomicRanges)
library(pheatmap)
library(SnapATAC);
library(viridisLite);
library(viridis)
library(ggpubr)
library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(MASS)
library(leiden)
library(ggplot2)
library(rhdf5)

#Load the aggregated kidney snap file & barcodes file
x.sp = createSnap(
  file="/node01data/dg/batch_effect/aggr_output/aggr.snap",
  sample="hkid_aggr",
  num.cores=1
);
barcodes = read.csv(
  "/node01data/dg/batch_effect/aggr_output/kidney_ATAC_aggr/outs/singlecell.csv",
  head=TRUE
);

# Add the sample number to snap object   
ATAC_sample <- sapply(as.character(x.sp@barcode), function(x) {  strsplit(x, split = "-")[[1]][2]}, simplify = TRUE)
x.sp@metaData$sample_no <- ATAC_sample

# Barcode selection 
barcodes = barcodes[2:nrow(barcodes),];
promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
UMI = log(barcodes$passed_filters+1, 10);
FRIP = barcodes$peak_region_fragments/barcodes$passed_filters
data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
barcodes$promoter_ratio = promoter_ratio;
data_cutoff <- data
data_cutoff <- data[which(UMI >= 3 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6 & FRIP > 0.3),];
data_cutoff$cut_off <- "cut_off"
ggplot() + geom_point(data=data, aes(x=UMI, y=promoter_ratio),color="grey", size=0.1) +
  geom_point(data=data_cutoff, aes(x=UMI, y=promoter_ratio),color = "red", size=0.1)  + theme(text = element_text(size=20))+
  ggtitle("") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio") +
  geom_hline(yintercept = 0.15, colour="black", linetype="dashed", size = 0.1) + 
  geom_hline(yintercept = 0.6, colour="black", linetype="dashed", size = 0.1)+
  geom_vline(xintercept = 3, colour="black", linetype="dashed", size = 0.1)+
  geom_vline(xintercept = 5, colour="black", linetype="dashed", size = 0.1)
barcodes.sel = barcodes[which(UMI >= 3 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6 & FRIP > 0.3),];
rownames(barcodes.sel) = barcodes.sel$barcode;
x.sp = x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),];
x.sp@metaData = barcodes.sel[x.sp@barcode,];

# Add cell-by-bin matrix with filtering unwanted bin (Atfer add the bin size to snap file in linux snaptools)
showBinSizes("/node01data/dg/batch_effect/aggr_output/aggr.snap");
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=10);
x.sp = makeBinary(x.sp, mat="bmat");
black_list = read.table("/node01data/dg/human_Kidney_ATAC_normal_disease_raw/HN00142306_10X_RawData_Outs/hg38.blacklist.bed.gz");
black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3]));
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse = "|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];

# Dimension reduction & clustering with correcting batch effect 
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
);
library(harmony)
x.sp@sample <- x.sp@metaData$sample_no
x.sp = runHarmony(
  obj=x.sp, 
  eigs.dim=1:16, 
  meta_data=x.sp@sample # sample index
);
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:16,
  k=15
);
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  seed.use=10,
  resolution=2.3
);
x.sp@metaData$cluster = x.sp@cluster;
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:16, 
  method="umap",
  seed.use=10
);

# Visulization
Palette <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(4, "Paired"))
df <- data.frame(UMAP1 = x.sp@umap[, 1],
                 UMAP2 = x.sp@umap[, 2],
                 sample = x.sp@metaData$sample_no,
                 cluster = x.sp@cluster, 
                 Read_Depth = log(x.sp@metaData[,"passed_filters"]+1,10),
                 FRiP = x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters)
ggplot(df, aes(x = UMAP1 , y = UMAP2, col = cluster)) +
  geom_point(size = 0.1, alpha = 1) + 
  scale_color_manual(values = Palette[1:20]) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) +
  coord_fixed(ratio = 1)
  
# Add the cell-by-gene matrix
genes = read.table("/node01data/dg/heritage_2/previous_infile/gencode.v33.annotation.bed")
x.sp = addBmatToSnap(x.sp)
genes.gr = GRanges(genes[, 1], 
                   IRanges(genes[, 2], genes[, 3]), name = genes[, 5])
                   marker.genes = c("KDR", "NPHS2", "SLC5A2", "UMOD", "SLC12A3", "TRPV5", "AQP2", "AQP6", "SLC26A4", "COL1A1", "C1QA", "CD3E");
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
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
  obj=x.sp,
  input.mat="gmat",
  step.size=3
);

# correlation among cell type (heretical cluster & correlation heatmap)
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
plot(hc, hang=-1, xlab="Cell Type");
dist_mat = matrix(cor(t(do.call(rbind, ensemble.ls))),
                  nrow = 12,
                  ncol = 12,
                  byrow = TRUE)
rownames(dist_mat) <- levels(x.sp@cluster)
colnames(dist_mat) <- levels(x.sp@cluster)
library(gplots)
pheatmap(mat = dist_mat,
         color = bluered(100),
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         treeheight_col = 0,
         cellwidth = 20,
         cellheight = 20,
         show_colnames = FALSE)

# Peak calling
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 200)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/opt/anaconda2/bin/snaptools",
    path.to.macs="/node01data/dg/anaconda3/bin/macs2",
    gsize="hs",
    buffer.size=500, 
    num.cores=1,
    macs.options="--nomodel --shift -75 --ext 150 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=15);
peaks.names = system("ls | grep narrowPeak", intern=TRUE);

peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
            
# Add cell-by-peak matrix to snap object (Atfer add the peak list to snap file in linux snaptools)
x.sp = addPmatToSnap(x.sp);
x.sp = makeBinary(x.sp, mat="pmat");

# Calling DAR
idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    bcv=0.4,
    test.method="exactTest",
    seed.use=10
  );
  DARs$FDR = p.adjust(DARs$PValue, method="BH");
  idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue;
    PValues[DARs$logFC < 0] = 1;
    idy = order(PValues, decreasing=FALSE)[1:2000];
    rm(PValues); # free memory
  }
  idy
})
idy.ls
names(idy.ls) = levels(x.sp@cluster)
