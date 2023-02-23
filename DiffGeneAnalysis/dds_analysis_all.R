library("DESeq2")
library("IHW")
library("EnhancedVolcano")
library("vsn")
library("RColorBrewer")
library("pheatmap")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("reshape")
library("biomaRt")
defaultW <- getOption("warn") 

fname<- paste0("read.counts/all")
dir.create(file.path(fname)) 

total.counts <- read.table("read.counts/renamed.total.reads.txt", sep = "\t", header = T, row.names = 1)
targets <- read.table(paste0("read.counts/targets/experiment.renamed.txt"), header = T, sep = "\t", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = total.counts, colData = targets,design = ~condition)
dds$condition <- relevel(dds$condition, ref = targets$condition[1])

dds <- DESeq(dds)
res <- results(dds, filterFun=ihw)

# Remove NaNs
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$pvalue),]

png(paste0(fname,"/ma_plot_all.png"))
plotMA(res, ylim=c(-2,2))
#ggsave(paste0(fname,"/ma_plot_",set1,"_vs_",set2,".pdf"), plot=maplot, width=10, height=10, units="in")
dev.off()

# Summarize
summary(res)

# Plot volcanos

# volc<-EnhancedVolcano(res, lab = res@rownames, x = 'log2FoldChange',y = 'pvalue', xlim = c(-12,12),
                      # title = paste0("DGE all ", " (reference ",targets$condition[1],")"),
                      # subtitle=" ",
                      # pCutoff = 0.05,
                      # FCcutoff = 1,
                      # pointSize = 3.0,
                      # labSize = 3.0,    
                      # legendLabels=c('Not sig.',
                                     # 'Log (base 2) FC','p-value',
                                     # 'p-value & Log (base 2) FC'),
                      # legendPosition = 'right',
                      # legendLabSize = 12,
                      # legendIconSize = 4.0,
                      # labCol = 'black',
                      # labFace = 'bold',
                      #boxedLabels = TRUE,
                      # colAlpha = 4/5,
                      # xlab = bquote(~Log[2]~ 'fold change'),
                      # selectLab = c('ENSG00000169896',
                                    # 'ENSG00000129654',
                                    # 'ENSG00000162692',
                                    # 'ENSG00000119508',
                                    # 'ENSG00000142319',
                                    # 'ENSG00000179520'),
                      # drawConnectors = TRUE,
                      # widthConnectors = 0.25,
                      # colConnectors = 'black'
# )
# ggsave(paste0(fname,"/","DGE_all.png"), plot=volc, width=15, height=9, units="in")

rld <- rlog(dds,blind = FALSE)

# meanSdPlot(assay(rld))
# dev.off()
dist.rld <- dist(t(assay(rld)))
matrix.dist.rld <- as.matrix(dist.rld)


rownames(matrix.dist.rld) <- paste(rld$condition)
colnames(matrix.dist.rld) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(paste0(fname,"/","pheatmap_rld_all.png"))
pheatmap(matrix.dist.rld,
         clustering_distance_rows=dist.rld,
         clustering_distance_cols=dist.rld,
         col=colors)
dev.off()

# PCA
#pdf(paste0(fname,"/","pca_rld_",set1,"_vs_",set2,".pdf"))
plotpca<-plotPCA(rld, intgroup=c("condition"))
ggsave(file=paste0(fname,"/","allPCA.png"), plot=plotpca, width=7, height=7, units="in")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
mat <- assay(rld[select,]) 
mat.res <- as.matrix(res)

# Genes appeared in up- and down-regulated GO pathways in C2 vs Pac0.5
# pdf(paste0(fname,"/","ENSG00000169896.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000169896") 
# dev.off()
# pdf(paste0(fname,"/","ENSG00000129654.pdf"),  width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000129654") 
# dev.off()
# pdf(paste0(fname,"/","ENSG00000162692.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000162692") 
# dev.off()
# pdf(paste0(fname,"/","ENSG00000119508.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000119508") 
# dev.off()
# pdf(paste0(fname,"/","ENSG00000142319.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000142319") 
# dev.off()
# pdf(paste0(fname,"/","ENSG00000179520.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000179520") 
# dev.off()


## Plot gene counts for type of migration
# pdf(paste0(fname,"/","E-cadherin.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000039068") 
# dev.off()
# pdf(paste0(fname,"/","N-cadherin.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000170558") 
# dev.off()
# pdf(paste0(fname,"/","Vimentin.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000187650") 
# dev.off()
# pdf(paste0(fname,"/","Twist1.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000122691") 
# dev.off()
# pdf(paste0(fname,"/","Snail.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000124216") 
# dev.off()
# pdf(paste0(fname,"/","Slug.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000019549") 
# dev.off()
# pdf(paste0(fname,"/","Zeb1.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000148516") 
# dev.off()
# pdf(paste0(fname,"/","KRAS.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000133703") 
# dev.off()
# pdf(paste0(fname,"/","RB1.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000139687") 
# dev.off()
# pdf(paste0(fname,"/","ARP3.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000249577") 
# dev.off()
# pdf(paste0(fname,"/","RAP1.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000229509") 
# dev.off()
# pdf(paste0(fname,"/","KIF14.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000118193") 
# dev.off()
# pdf(paste0(fname,"/","EZR.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000092820") 
# dev.off()
# pdf(paste0(fname,"/","TRPM7.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000092439") 
# dev.off()
# pdf(paste0(fname,"/","AKT3.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000117020") 
# dev.off()
# pdf(paste0(fname,"/","AKT1.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000142208") 
# dev.off()
# pdf(paste0(fname,"/","AKT2.pdf"), width=7, height=7)
# d<-plotCounts(dds, gene="ENSG00000105221") 
# dev.off()


# t_mig_genes <- c("ENSG00000039068","ENSG00000170558","ENSG00000187650",
#                  "ENSG00000122691","ENSG00000124216","ENSG00000019549",
#                  "ENSG00000148516","ENSG00000133703","ENSG00000139687",
#                  "ENSG00000249577","ENSG00000229509","ENSG00000118193","ENSG00000092820")


# sp.genes <- c("PDLIM7", "MYH9") #these genes are commonly appeared in:
# #"cell-substrate junction" "focal adhesion" "cell leading edge" "adherens junction" "microtubule"  
# sp.genes <- c(sp.genes,"PALLD","PXN","ACTN1")  
# #these genes are commonly appeared in:
# #"cell-substrate junction" "focal adhesion" "cell leading edge" "adherens junction"
# #together with the previous genes
# hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
# mapping <- getBM(
#   attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
#   filters = 'hgnc_symbol',
#   values = sp.genes,
#   mart = hsmart
# )
# 
# pdf(paste0(fname,"/",mapping$hgnc_symbol[1],"pdf"), width=7, height=7)
# d<-plotCounts(dds, gene=mapping$ensembl_gene_id[1]) 
# dev.off()
# pdf(paste0(fname,"/",mapping$hgnc_symbol[2],"pdf"), width=7, height=7)
# d<-plotCounts(dds, gene=mapping$ensembl_gene_id[2]) 
# dev.off()
# pdf(paste0(fname,"/",mapping$hgnc_symbol[3],"pdf"), width=7, height=7)
# d<-plotCounts(dds, gene=mapping$ensembl_gene_id[3]) 
# dev.off()
# pdf(paste0(fname,"/",mapping$hgnc_symbol[4],"pdf"), width=7, height=7)
# d<-plotCounts(dds, gene=mapping$ensembl_gene_id[4]) 
# dev.off()
# pdf(paste0(fname,"/",mapping$hgnc_symbol[5],"pdf"), width=7, height=7)
# d<-plotCounts(dds, gene=mapping$ensembl_gene_id[5]) 
# dev.off()
# genes <- c('ENSG00000169896',
#            'ENSG00000129654',
#            'ENSG00000162692',
#            'ENSG00000119508',
#            'ENSG00000142319',
#            'ENSG00000179520',
#            mapping$ensembl_gene_id,
#            mapping$ensembl_gene_id)
# ndds<-dds[genes,]
# normalized_counts <- as.data.frame(counts(ndds, normalized=TRUE)) 
# normalized_counts$genes <- row.names(normalized_counts)
# melted_counts <- melt(normalized_counts)
# 
# goi<-ggplot(melted_counts) +
#   geom_point(aes(x = variable, y = value, color = genes), position=position_jitter(w=0.1,h=0)) +
#   geom_line(aes(x = variable, y = value, color = genes)) +
#   scale_y_log10() +
#   xlab("Sample") +
#   ylab("Normalized Counts") +
#   ggtitle("Genes of interest") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(plot.title=element_text(hjust=0.5))
# ggsave(paste0(fname,"/","goi.pdf"), plot=goi, width=7, height=7, units="in")



## MAPK pathway genes
mapk <- c("AKT3", "CSF1", "DUSP1", "EFNA1", "EFNA5", "ELK4", "EPHA2", "EREG", "FGF5", "FGFR4", "HSPA1A", 
          "HSPA1B", "HSPA6", "MAP3K20", "MAP3K7", "MAPK9", "PDGFC", "RPS6KA1", "RPS6KA2", "SRF", "TGFB2", "TGFBR2")

mapk.mapping <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mapk,
  uniqueRows=TRUE,
  mart = hsmart
)
mapk.mapping <- mapk.mapping[-c(1,12,13,14,15,17,18,19,20,27),]
mapk.dds<-dds[mapk.mapping$ensembl_gene_id,]

mapk.counts       <- as.data.frame(counts(mapk.dds, normalized=TRUE)) 
mapk.counts$max   <- apply(mapk.counts[,c(1:12)], 1, max) 
mapk.counts       <- mapk.counts/mapk.counts$max
mapk.counts       <- subset(mapk.counts, select=-max)
mapk.counts$genes <- mapk.mapping$hgnc_symbol[row.names(mapk.counts) %in% mapk.mapping$ensembl_gene_id]
mapk.melted_counts <- melt(mapk.counts)
colnames(mapk.melted_counts) = c("genes","sample","norm.count")

heat.mapk.goi<-ggplot(mapk.melted_counts, aes(sample, genes)) + geom_tile(aes(fill=norm.count)) +
  ggtitle("MAPK signalling pathway")+
  theme_classic(base_size=20)+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  scale_fill_distiller(palette = "RdBu")+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8))
ggsave(paste0(fname,"/","mapk_heat_norm.pdf"), plot=heat.mapk.goi, width=7, height=7, units="in")

mapk.counts_mean <- mapk.counts 
mapk.counts_mean$C1<- rowMeans(subset(mapk.counts, select = c(1,2)))
mapk.counts_mean$C2<- rowMeans(subset(mapk.counts, select = c(3,4)))
mapk.counts_mean$Pac0.0005<- rowMeans(subset(mapk.counts, select = c(5,6)))
mapk.counts_mean$Pac0.005<- rowMeans(subset(mapk.counts, select = c(7,8)))
mapk.counts_mean$Pac0.05<- rowMeans(subset(mapk.counts, select = c(9,10)))
mapk.counts_mean$Pac0.5<- rowMeans(subset(mapk.counts, select = c(11,12)))
mapk.counts_mean  <- subset(mapk.counts_mean ,select=-c(1:12))
mapk.melted_counts_mean <- melt(mapk.counts_mean)
colnames(mapk.melted_counts_mean) = c("genes","sample","norm.count")

heat.mapk.goi.mean<-ggplot(mapk.melted_counts_mean, aes(x=sample, y=genes)) + geom_tile(aes(fill=norm.count)) +
  ggtitle("MAPK signalling pathway")+
  theme_classic(base_size=20)+
  scale_fill_distiller(palette = "RdBu")+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  #theme(axis.text.x = element_text(angle = 10)) +
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8))
ggsave(paste0(fname,"/","mapk_heat_norm_mean.pdf"), plot=heat.mapk.goi.mean, width=7, height=7, units="in")


## TGF-beta pathway
tgf <- c("ATF3", "BRCA1", "CCNE1", "CDK1", "FOSB", "MAP2K6", "RBL1", "SIK1", "TP73")
tgf.mapping <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = tgf,
  uniqueRows=TRUE,
  mart = hsmart
)
#tgf.mapping <- mapk.mapping[-c(1,12,13,14,15,17,18,19,20,27),]
tgf.dds<-dds[tgf.mapping$ensembl_gene_id,]

tgf.counts       <- as.data.frame(counts(tgf.dds, normalized=TRUE)) 
tgf.counts$max   <- apply(tgf.counts[,c(1:12)], 1, max) 
tgf.counts       <- tgf.counts/tgf.counts$max
tgf.counts       <- subset(tgf.counts, select=-max)
tgf.counts$genes <- tgf.mapping$hgnc_symbol[row.names(tgf.counts) %in% tgf.mapping$ensembl_gene_id]
tgf.melted_counts <- melt(tgf.counts)
colnames(tgf.melted_counts) = c("genes","sample","norm.count")

heat.tgf<-ggplot(tgf.melted_counts, aes(sample, genes)) + geom_tile(aes(fill=norm.count)) +
  ggtitle("TGF-beta signalling pathway")+
  theme_classic(base_size=20)+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  scale_fill_distiller(palette = "RdBu")+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8))
ggsave(paste0(fname,"/","tgf_heat_norm.pdf"), plot=heat.tgf, width=7, height=7, units="in")

tgf.counts_mean <- tgf.counts 
tgf.counts_mean$C1<- rowMeans(subset(tgf.counts, select = c(1,2)))
tgf.counts_mean$C2<- rowMeans(subset(tgf.counts, select = c(3,4)))
tgf.counts_mean$Pac0.0005<- rowMeans(subset(tgf.counts, select = c(5,6)))
tgf.counts_mean$Pac0.005<- rowMeans(subset(tgf.counts, select = c(7,8)))
tgf.counts_mean$Pac0.05<- rowMeans(subset(tgf.counts, select = c(9,10)))
tgf.counts_mean$Pac0.5<- rowMeans(subset(tgf.counts, select = c(11,12)))
tgf.counts_mean  <- subset(tgf.counts_mean ,select=-c(1:12))
tgf.melted_counts_mean <- melt(tgf.counts_mean)
colnames(tgf.melted_counts_mean) = c("genes","sample","norm.count")

heat.tgf.mean<-ggplot(tgf.melted_counts_mean, aes(x=sample, y=genes)) + geom_tile(aes(fill=norm.count)) +
  ggtitle("TGF-beta signalling pathway")+
  theme_classic(base_size=20)+
  scale_fill_distiller(palette = "RdBu")+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  #theme(axis.text.x = element_text(angle = 10)) +
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8))
ggsave(paste0(fname,"/","tgf_heat_norm_mean.pdf"), plot=heat.tgf.mean, width=7, height=7, units="in")




## Structural analysis
af_genes <- read.csv(file = "genes_for_actin_filaments_C2vsPac0_5.csv")
colnames(af_genes) <- "actin_filaments"
fas_genes <- read.csv(file = "genes_for_focal_adhesion_sites_C2vsPac0_5.csv")
colnames(fas_genes) <- "focal_adhesion_sites"
pm_genes <- read.csv(file = "genes_for_plasma_membrane_C2vsPac0_5.csv")
colnames(pm_genes) <- "plasma_membrane"
sp_genes <- read.csv(file = "genes_for_secreted_proteins_C2vsPac0_5.csv")
colnames(sp_genes) <- "secreted_proteins"

total_genes <- c(af_genes$actin_filaments,fas_genes$focal_adhesion_sites,pm_genes$plasma_membrane,sp_genes$secreted_proteins)
unique_genes <- unique(total_genes)

map_genes <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = unique_genes,
  mart = hsmart
)

all_genes <- dds@rowRanges@partitioning@NAMES
test      <- map_genes$ensembl_gene_id %in% all_genes
filtered_map_genes <- map_genes[test,]
gdds      <- dds[filtered_map_genes$ensembl_gene_id,]
nottest   <- all_genes %in% map_genes$ensembl_gene_id
notdds    <- dds[!nottest,]

nt <- results(notdds)
nt <- nt[!is.na(nt$padj),]
nt <- nt[!is.na(nt$pvalue),]
nt <- nt[nt$pvalue<0.05,]

notdds <- dds[rownames(nt),]


library("tidyr")
library("dplyr")

normalized_counts_genes     <- as.data.frame(counts(gdds, normalized=TRUE)) 
normalized_counts_genes$max <- apply(normalized_counts_genes[,c(1:12)], 1, max)
normalized_counts_genes     <- normalized_counts_genes/normalized_counts_genes$max
normalized_counts_genes     <- subset(normalized_counts_genes, select=-max)
normalized_counts_genes_mean <- normalized_counts_genes
normalized_counts_genes_mean$C1<- rowMeans(subset(normalized_counts_genes, select = c(1,2)))
normalized_counts_genes_mean$C2<- rowMeans(subset(normalized_counts_genes, select = c(3,4)))
normalized_counts_genes_mean$Pac0.0005<- rowMeans(subset(normalized_counts_genes, select = c(5,6)))
normalized_counts_genes_mean$Pac0.005<- rowMeans(subset(normalized_counts_genes, select = c(7,8)))
normalized_counts_genes_mean$Pac0.05<- rowMeans(subset(normalized_counts_genes, select = c(9,10)))
normalized_counts_genes_mean$Pac0.5<- rowMeans(subset(normalized_counts_genes, select = c(11,12)))
normalized_counts_genes_mean <- subset(normalized_counts_genes_mean,select=-c(1:12))
#normalized_counts_genes_mean <- normalized_counts_genes_mean[order(normalized_counts_genes_mean$C1, decreasing = TRUE), ]

normalized_counts_genes_rest     <- as.data.frame(counts(notdds, normalized=TRUE)) 
normalized_counts_genes_rest$max <- apply(normalized_counts_genes_rest[,c(1:12)], 1, max)
normalized_counts_genes_rest     <- normalized_counts_genes_rest/normalized_counts_genes_rest$max
normalized_counts_genes_rest     <- subset(normalized_counts_genes_rest, select=-max)
normalized_counts_genes_mean_rest <- normalized_counts_genes_rest
normalized_counts_genes_mean_rest$C1<- rowMeans(subset(normalized_counts_genes_rest, select = c(1,2)))
normalized_counts_genes_mean_rest$C2<- rowMeans(subset(normalized_counts_genes_rest, select = c(3,4)))
normalized_counts_genes_mean_rest$Pac0.0005<- rowMeans(subset(normalized_counts_genes_rest, select = c(5,6)))
normalized_counts_genes_mean_rest$Pac0.005<- rowMeans(subset(normalized_counts_genes_rest, select = c(7,8)))
normalized_counts_genes_mean_rest$Pac0.05<- rowMeans(subset(normalized_counts_genes_rest, select = c(9,10)))
normalized_counts_genes_mean_rest$Pac0.5<- rowMeans(subset(normalized_counts_genes_rest, select = c(11,12)))
normalized_counts_genes_mean_rest <- subset(normalized_counts_genes_mean_rest,select=-c(1:12))
#normalized_counts_genes_mean_rest <- normalized_counts_genes_mean_rest[order(normalized_counts_genes_mean_rest$C1, decreasing = TRUE), ]


mat <- as.data.frame(normalized_counts_genes_mean)
mat$genes <- rownames(mat)
mat <- mat %>% gather(key='samples', value='value', -genes)

mat_rest <- as.data.frame(normalized_counts_genes_mean_rest)
mat_rest$genes <- rownames(mat_rest)
mat_rest <- mat_rest %>% gather(key='samples', value='value', -genes)

heat_str<-ggplot(mat, aes(samples, genes)) + geom_tile(aes(fill=value)) + 
          ggtitle("Expression of genes related to actin filaments, focal adhesion sites,\n plasma membrane, protein secretion")+
          theme_classic(base_size=24)+
          #remove extra space
          #scale_alpha( trans = "log" ) +
          theme(axis.text.x = element_text(angle = 90)) +
          theme(axis.text.y=element_text(size=3))

ggsave(file=paste0(fname,"/","heatmap_genes_for_structure.pdf"), plot=heat_str, width=18, height=18, units="in")


heat_rest<-ggplot(mat_rest, aes(samples, genes)) + geom_tile(aes(fill=value)) + 
  ggtitle("Rest of the differentially expressed genes")+
  theme_classic(base_size=24)+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank()) 

ggsave(file=paste0(fname,"/","heatmap_genes_for_rest.pdf"), plot=heat_rest, width=18, height=18, units="in")



## Type of migration: EMT
indiv <- c("TWIST1","SNAI1","SNAI2","CDH1","CDH2","ZEB1","VMAC")

map_mig_genes <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = indiv,
  mart = hsmart
)

emt.dds<-dds[map_mig_genes$ensembl_gene_id,]

emt.counts       <- as.data.frame(counts(emt.dds, normalized=TRUE)) 
emt.counts$max   <- apply(emt.counts[,c(1:12)], 1, max) 
emt.counts       <- emt.counts/emt.counts$max
emt.counts       <- subset(emt.counts, select=-max)
emt.counts$genes <- map_mig_genes$hgnc_symbol[row.names(emt.counts) %in% map_mig_genes$ensembl_gene_id]
emt.melted_counts <- melt(emt.counts)
colnames(emt.melted_counts) = c("genes","sample","norm.count")

heat.emt<-ggplot(emt.melted_counts, aes(sample, genes)) + geom_tile(aes(fill=norm.count)) +
  ggtitle("EMT related genes")+
  theme_classic(base_size=20)+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  scale_fill_distiller(palette = "RdBu")+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8))
ggsave(paste0(fname,"/","emt_heat_norm.pdf"), plot=heat.emt, width=7, height=7, units="in")

emt.counts_mean <- emt.counts 
emt.counts_mean$C1<- rowMeans(subset(emt.counts, select = c(1,2)))
emt.counts_mean$C2<- rowMeans(subset(emt.counts, select = c(3,4)))
emt.counts_mean$Pac0.0005<- rowMeans(subset(emt.counts, select = c(5,6)))
emt.counts_mean$Pac0.005<- rowMeans(subset(emt.counts, select = c(7,8)))
emt.counts_mean$Pac0.05<- rowMeans(subset(emt.counts, select = c(9,10)))
emt.counts_mean$Pac0.5<- rowMeans(subset(emt.counts, select = c(11,12)))
emt.counts_mean  <- subset(emt.counts_mean ,select=-c(1:12))
emt.melted_counts_mean <- melt(emt.counts_mean)
colnames(emt.melted_counts_mean) = c("genes","sample","norm.count")

heat.emt.mean<-ggplot(emt.melted_counts_mean, aes(x=sample, y=genes)) + geom_tile(aes(fill=norm.count)) +
  ggtitle("EMT related genes")+
  theme_classic(base_size=20)+
  scale_fill_distiller(palette = "RdBu")+
  #remove extra space
  #scale_alpha( trans = "log" ) +
  #theme(axis.text.x = element_text(angle = 10)) +
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8))
ggsave(paste0(fname,"/","emt_heat_norm_mean.pdf"), plot=heat.emt.mean, width=7, height=7, units="in")
