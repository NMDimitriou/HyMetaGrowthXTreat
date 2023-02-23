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
library("plyr")
defaultW <- getOption("warn") 
options(warn = -1)


args = commandArgs(trailingOnly=TRUE)
set1 <- args[1]
set2 <- args[2]

set1 <- "C2"
set2 <- "Pac0.5"

fname<- paste0("read.counts/",set1,"_vs_",set2)
#
dir.create(file.path(fname)) 

total.counts <- read.table("read.counts/renamed.total.reads.txt", sep = "\t", header = T, row.names = 1)

c1 <- grep(set1,colnames(total.counts), fixed = TRUE)
c2 <- grep(set2,colnames(total.counts), fixed = TRUE)
iso.counts <- as.matrix(total.counts[,c(c1,c2)])

targets <- read.table(paste0("read.counts/targets/target_",set1,"_vs_",set2,".txt"), header = T, sep = "\t", row.names = 1)
#

###########################################################
# Differential gene expression analysis
##########################################################

# Construct the input data set from compact files (_counts_compact.csv, HeLa_counts_compact.csv)
dds <- DESeqDataSetFromMatrix(countData = iso.counts, colData = targets,design = ~condition)
#
dds$condition <- relevel(dds$condition, ref = targets$condition[1])
# Run the DGE analysis
dds <- DESeq(dds)

# Filter results with IHW (see Ignatiadis et al.)
res <- results(dds, filterFun=ihw)

# Remove NaNs
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$pvalue),]


png(paste0(fname,"/ma_plot_",set1,"_vs_",set2,".png"))
plotMA(res, ylim=c(-2,2))
#ggsave(paste0(fname,"/ma_plot_",set1,"_vs_",set2,".pdf"), plot=maplot, width=10, height=10, units="in")
dev.off()

# Plot Histograms of p-values
hist(res$pvalue, main="res.", xlab="p-values",cex.main=1.5)
dev.off()
# Summarize
summary(res)

# Change annotation from ensembl to gene symbol for better visualization
####   mat.res <- as.matrix(res)
#mat.res <- merge(mat.res,geneIDs,by.x="row.names",by.y="V1")

# Plot volcanos
volc<-EnhancedVolcano(res, lab = res@rownames, x = 'log2FoldChange',y = 'pvalue', xlim = c(-12,12),
                title = paste0("DGE ", set1,"_vs_",set2," (reference ",set1,")"),
                subtitle=" ",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,    
                legendLabels=c('Not sig.',
                               'Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                colAlpha = 4/5,
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = c('ENSG00000169896',
                              'ENSG00000129654',
                              'ENSG00000162692',
                              'ENSG00000119508',
                              'ENSG00000142319',
                              'ENSG00000179520'),
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                colConnectors = 'black'
                )
ggsave(paste0(fname,"/","DGE ", set1,"_vs_",set2,".png"), plot=volc, width=15, height=9, units="in")

# Apply regularized trnasformation to plot heatmaps on the distance between samples
rld <- rlog(dds,blind = FALSE)

meanSdPlot(assay(rld))
dist.rld <- dist(t(assay(rld)))
matrix.dist.rld <- as.matrix(dist.rld)


rownames(matrix.dist.rld) <- paste(rld$condition)
colnames(matrix.dist.rld) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(paste0(fname,"/","pheatmap_rld_",set1,"_vs_",set2,".png"))
#png(paste0(fname,"/","pheatmap_rld_all.png"))
pheatmap(matrix.dist.rld,
         clustering_distance_rows=dist.rld,
         clustering_distance_cols=dist.rld,
         col=colors)
dev.off()

# PCA
#pdf(paste0(fname,"/","pca_rld_",set1,"_vs_",set2,".pdf"))
plotpca<-plotPCA(rld, intgroup=c("condition"))
ggsave(paste0(fname,"/","pca_rld_",set1,"_vs_",set2,".png"), plot=plotpca, width=7, height=3, units="in")
#ggsave(file="allPCA.png", plot=plotpca, width=7, height=7, units="in")
#dev.off()
# Plot heatmaps for the top 20 differential expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
mat <- assay(rld[select,]) 
#mat <- merge(mat,geneIDs,by.x="row.names",by.y="V1")
#mat<-mat[-c(1)]
#rownames(mat) <- mat$V2
#mat<-mat[-c(7)]
df <- as.data.frame(colData(dds)["condition"])

#png(paste0(fname,"/","pheatmap_mat_",set1,"_vs_",set2,".png"))
#pmap<-pheatmap(mat, annotation_col = df)
#ggsave(paste0(fname,"/","pheatmap_mat_",set1,"_vs_",set2,".png"), plot=pmap, width=10, height=10, units="in")
#dev.off()

mat.res <- as.matrix(res)

######################################################
# Enrichment analysis
######################################################
#colnames(mat.res)[1] <- "ENSEMBL"
#colnames(mat.res)[9] <- "SYMBOL"
#gene <- mat.res[,c(1,9)]
#gene$ENSEMBL <- substr(mat.res$ENSEMBL,1,15)

################################################
# GO O.R. TEST
################################################

#Select up and down regulated genes
up   <- res[which(res$log2FoldChange > 0 & res$padj < .05),]
down <- res[which(res$log2FoldChange < 0 & res$padj < .05),]
write.table(up  , file=paste0(fname,"/","upregulated_",set1,"_vs_",set2,".txt")  , sep = "\t", row.names = T)
write.table(down, file=paste0(fname,"/","downregulated_",set1,"_vs_",set2,".txt"), sep = "\t", row.names = T)

up.gene   <- as.matrix(row.names(up))
down.gene <- as.matrix(row.names(down))
#gene <- as.matrix(row.names(mat.res))

# GO over-representation test. Here, you replace "gene" with "up.genes" and "down.genes" for the enrichment
ego.up <- enrichGO(gene           = up.gene,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    #qvalueCutoff  = 0.05,
                    readable      = TRUE)

ego.down <- enrichGO(gene        = down.gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   #qvalueCutoff  = 0.05,
                   readable      = TRUE)

#Visualize
#pdf(paste0(fname,"/","barplot_ego_up_",set1,"_vs_",set2,".pdf"))
#barplot(ego.up,showCategory = 20) + ggtitle(paste0("GO ", set1,"_vs_",set2," up in (",set2,")") )
#dev.off()

#pdf(paste0(fname,"/","barplot_ego_down_",set1,"_vs_",set2,".pdf"))
#barplot(ego.down,showCategory = 20) + ggtitle(paste0("GO ", set1,"_vs_",set2," up in (",set2,")"))
#dev.off()

################################################
# GO O.R Visualization
################################################

#pdf(paste0(fname,"/","dotplot_go_up_",set1,"_vs_",set2,".pdf"))
try(dotup<-dotplot(ego.up, showCategory=30) + ggtitle(paste0("GO ", set1,"_vs_",set2," (up in ",set2,")")), silent = TRUE)
try(ggsave(paste0(fname,"/","dotplot_go_up_",set1,"_vs_",set2,".png"), plot=dotup, width=10, height=10, units="in"), silent = TRUE)
#dev.off()

#pdf(paste0(fname,"/","dotplot_go_down_",set1,"_vs_",set2,".pdf"))
try(dotdown<-dotplot(ego.down, showCategory=30) + ggtitle(paste0("GO ", set1,"_vs_",set2," (down in ",set2,")")), silent = TRUE)
try(ggsave(paste0(fname,"/","dotplot_go_down_",set1,"_vs_",set2,".png"), plot=dotdown, width=10, height=10, units="in"), silent = TRUE)
#dev.off()

goi_terms=c("cell-substrate junction","focal adhesion",
            "cell leading edge","adherens junction","microtubule",
            "actomyosin")


allterms.ego.down <- ego.down@result$Description

allgenes.ego.down <- ego.down@result$geneID[allterms.ego.down %in% goi_terms]
allgenes.ego.down <- strsplit(allgenes.ego.down,"/")
intersect.genes.down<-Reduce(intersect,list(allgenes.ego.down[[1]],allgenes.ego.down[[2]],allgenes.ego.down[[3]],allgenes.ego.down[[4]]))
lg.ego.down <- ldply (allgenes.ego.down, data.frame)

unique.count = function(train.data, all.numeric=FALSE) {                                                                                                                                                                                                 
  # first convert each row in the data.frame to a string                                                                                                                                                                              
  train.data.str = apply(train.data, 1, function(x) paste(x, collapse=','))                                                                                                                                                           
  # use table to index and count the strings                                                                                                                                                                                          
  train.data.str.t = table(train.data.str)                                                                                                                                                                                            
  # get the unique data string from the row.names                                                                                                                                                                                     
  train.data.str.uniq = row.names(train.data.str.t)                                                                                                                                                                                   
  weight = as.numeric(train.data.str.t)                                                                                                                                                                                               
  # convert the unique data string to data.frame
  if (all.numeric) {
    train.data.uniq = as.data.frame(t(apply(cbind(train.data.str.uniq), 1, 
                                            function(x) as.numeric(unlist(strsplit(x, split=","))))))                                                                                                    
  } else {
    train.data.uniq = as.data.frame(t(apply(cbind(train.data.str.uniq), 1, 
                                            function(x) unlist(strsplit(x, split=",")))))                                                                                                    
  }
  names(train.data.uniq) = names(train.data)                                                                                                                                                                                          
  list(data=train.data.uniq, weight=weight)                                                                                                                                                                                           
}  


# GO Gene set enrichment analysis
#mat.res <- as.data.frame(mat.res)
#filter <- mat.res[mat.res$pvalue<0.5,]
#list <- filter[,2] # we add 2 (log2FoldChange)
#add ENSEMBL IDs as a separate column to the list object
#filter$ENSEMBL <- row.names(filter)
#names(list) <- as.character(substr(filter$ENSEMBL,1,15))
#GOlist <- sort(list,decreasing = TRUE)


################################################
# GSEA/GO
################################################

####For up-regulates genes in right sample (treatment)
up <- as.data.frame(up)
filter.up <- up[up$pvalue<0.5,]
list.up <- filter.up[,2]
filter.up$ENSEMBL <- row.names(filter.up)
names(list.up) <- as.character(substr(filter.up$ENSEMBL,1,15))
GOlist.up <- sort(list.up,decreasing = TRUE)

######For down-regulated genes in right sample:
down <- as.data.frame(down)
filter.down <- down[down$pvalue<0.5,]
list.down <- filter.down[,2]
filter.down$ENSEMBL <- row.names(filter.down)
names(list.down) <- as.character(substr(filter.down$ENSEMBL,1,15))
GOlist.down <- sort(list.down,decreasing = TRUE)
##

## GSEA
gsea.up <- gseGO(geneList = GOlist.up,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = 'ENSEMBL', 
                 ont          = "ALL",
                 #nPerm        = 100000,
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)


gsea.down <- gseGO(geneList = GOlist.down,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = 'ENSEMBL', 
                 ont          = "ALL",
                 #nPerm        = 100000,
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 #nPermSimple  = 1000000,
                 verbose      = FALSE)


################################################
# GSEA/GO Visualization
################################################

#pdf(paste0(fname,"/","heatplot_gsea_go_down_",set1,"_vs_",set2,".pdf"))
try(hup<-heatplot(gsea.up, foldChange = gsea.up@geneList, showCategory = 50)+
  ggplot2::coord_flip()+
  ggplot2::ylab('Annotations of Biological Processes')+
  ggplot2::xlab('Gene Symbols')+
  ggplot2::ggtitle(paste0('GO Over-representation for up-regulated pathways in ',set2))+
  ggplot2::theme(axis.text.y = element_text(size=2), axis.text.x = element_text(size = 6)), silent = TRUE)
try(ggsave(filename=paste0(fname,"/","heatplot_gsea_go_up_",set1,"_vs_",set2,".pdf"), plot=hup, width=10, height=10, units="in"), silent = TRUE)
#str(hup)
#dev.off()

try(hdown<-heatplot(gsea.down, foldChange = gsea.down@geneList,showCategory=50)+
  ggplot2::coord_flip()+
  ggplot2::ylab('Annotations of Biological Processes')+
  ggplot2::xlab('Gene Symbols')+
  ggplot2::ggtitle(paste0('GO Over-representation for down-regulated pathways in ',set2))+
  ggplot2::theme(axis.text.y = element_text(size=2), axis.text.x = element_text(size = 6)), silent = TRUE)
try(ggsave(filename=paste0(fname,"/","heatplot_gsea_go_down_",set1,"_vs_",set2,".pdf"), plot=hdown, width=10, height=10, units="in"), silent = TRUE)

#pdf(paste0(fname,"/","rigdeplot_gsea_go_up_",set1,"_vs_",set2,".pdf"))
try(rup<-ridgeplot(gsea.up,label_format = 30, core_enrichment = TRUE, showCategory = 30) + 
  ggtitle(paste0("GSEA/GO ", set1,"_vs_",set2," (up in ",set2,")")) +
  theme(text=element_text(size=20)), silent = TRUE)
try(ggsave(filename=paste0(fname,"/","rigdeplot_gsea_go_up_",set1,"_vs_",set2,".png"), plot=rup, width=18, height=10, units="in"), silent = TRUE)
#dev.off()

#pdf(paste0(fname,"/","rigdeplot_gsea_go_down_",set1,"_vs_",set2,"_copy.pdf"))
#ridgeplot(gsea.down,label_format = 30) + ggtitle(paste0("GSEA/GO ", set1,"_vs_",set2," (down in ",set2,")"))
try(rdown<-ridgeplot(gsea.down,label_format = 30, core_enrichment = TRUE, showCategory = 30) + 
  ggtitle(paste0("GSEA/GO ", set1,"_vs_",set2," (down in ",set2,")")) +
  theme(text=element_text(size=20)), silent = TRUE)
try(ggsave(filename=paste0(fname,"/","rigdeplot_gsea_go_down_",set1,"_vs_",set2,".png"), plot=rdown, width=18, height=10, units="in"), silent = TRUE)
#dev.off()

####################################################
# KEGG enrichment analysis
####################################################
#mat.res$ENSEMBL <- row.names(mat.res)
#for down-reg genes"
up$ENSEMBL   <- row.names(up)
down$ENSEMBL <- row.names(down)

#Replace mat.res$entrez with down$entrez for down-reg genes
up$entrez <- mapIds(org.Hs.eg.db,
                    keys = substr(up$ENSEMBL,1,15),
                    column = "ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals="first")

down$entrez <- mapIds(org.Hs.eg.db,
                             keys = substr(down$ENSEMBL,1,15),
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals="first")


#mat.res <- mat.res[!is.na(mat.res$entrez),]
#for down-reg genes:
up   <- up[!is.na(up$entrez),]
down <- down[!is.na(down$entrez),]

################################################
# KEGG O.R. TEST
################################################

# KEGG over-representation test
#replace with down for down-reg genes
kk.up   <- enrichKEGG(gene         = up$entrez,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)

kk.down <- enrichKEGG(gene        = down$entrez,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)


################################################
# KEGG O.R. Visualization
################################################

#pdf(paste0(fname,"/","dotplot_kegg_up_",set1,"_vs_",set2,".pdf"))
kkupdot<-dotplot(kk.up,showCategory = 20) + ggtitle(paste0("KEGG ", set1,"_vs_",set2," up in (",set2,")"))
ggsave(paste0(fname,"/","dotplot_kegg_up_",set1,"_vs_",set2,".png"), plot=kkupdot, width=10, height=10, units="in")
#dev.off()

#pdf(paste0(fname,"/","dotplot_kegg_down_",set1,"_vs_",set2,".pdf"))
kkdowndot<-dotplot(kk.down, showCategory=30) + ggtitle(paste0("KEGG ", set1,"_vs_",set2," down in (",set2,")"))
ggsave(paste0(fname,"/","dotplot_kegg_down_",set1,"_vs_",set2,".png"), plot=kkdowndot, width=10, height=10, units="in")
#dev.off()


################################################
# Onco enrichR
################################################
library(oncoEnrichR)
library(rmarkdown)
oeDB<-oncoEnrichR::load_db()
down.gene<-data.frame(down.gene)
up.gene<-data.frame(up.gene)
colnames(down.gene) <- "ensembl_gene_id"
colnames(up.gene)   <- "ensembl_gene_id"
down.onc_report <- oncoEnrichR::onco_enrich(query = down.gene$ensembl_gene_id, oeDB = oeDB,query_id_type = "ensembl_gene", show_cell_tissue = T, project_title = paste0(set1,"_vs_",set2," (down in ",set2,")"))
up.onc_report   <- oncoEnrichR::onco_enrich(query = up.gene$ensembl_gene_id  , oeDB = oeDB,query_id_type = "ensembl_gene", show_cell_tissue = T, project_title = paste0(set1,"_vs_",set2," (up in ",set2,")"))
try(oncoEnrichR::write(report = down.onc_report, file = paste0(fname,"/",set1,"_vs_",set2,"_down_report_oncoenrichr.html"), overwrite = T, format = "html"))
try(oncoEnrichR::write(report = up.onc_report  , file = paste0(fname,"/",set1,"_vs_",set2,"_up_report_oncoenrichr.html")  , overwrite = T, format = "html"))
################################################
# GSEA/KEGG 
################################################

# KEGG gene-set enrichment analysis

# filter.up   <- up[up$pvalue<0.5,]
# filter.down <- down[down$pvalue<0.5,]
# 
# filter.up   <- filter.up[!duplicated(filter.up$entrez),]
# filter.down <- filter.down[!duplicated(filter.down$entrez),]
# 
# list.up <- filter.up[,2] 
# names(list.up) <- as.character(filter.up$entrez)
# KEGGlist.up <- sort(list.up,decreasing = TRUE)
# 
# list.down <- filter.down[,2] #for down-reg genes
# names(list.down) <- as.character(filter.down$entrez)
# KEGGlist.down <- sort(list.down,decreasing = TRUE)
# 
# 
# gsek.up<- gseKEGG(geneList= KEGGlist.up,
#                   organism     = 'hsa',
#                   #nPerm        = 100000,
#                   minGSSize    = 120,
#                   pvalueCutoff = 0.05,
#                   verbose      = FALSE)
# 
# gsek.down <- gseKEGG(geneList= KEGGlist.down,
#                organism     = 'hsa',
#                #nPerm        = 100000,
#                minGSSize    = 120,
#                pvalueCutoff = 0.05,
#                verbose      = FALSE)


################################################
# GSEA/KEGG Visualization
################################################

#pdf(paste0(fname,"/","rigdeplot_gsea_kegg_up_",set1,"_vs_",set2,".pdf"))
#kkrup<-ridgeplot(gsek.up) + ggtitle(paste0("GSEA/KEGG ", set1,"_vs_",set2," up in (",set2,")"))
#ggsave(paste0(fname,"/","rigdeplot_gsea_kegg_up_",set1,"_vs_",set2,".png"), plot=kkrup, width=10, height=10, units="in")
#dev.off()

#pdf(paste0(fname,"/","rigdeplot_gsea_kegg_down_",set1,"_vs_",set2,".pdf"))
#kkrdown<-ridgeplot(gsek.down) + ggtitle(paste0("GSEA/KEGG ", set1,"_vs_",set2," down in (",set2,")"))
#ggsave(paste0(fname,"/","rigdeplot_gsea_kegg_down_",set1,"_vs_",set2,".pdf"), plot=kkrdown, width=10, height=10, units="in")
#dev.off()


#############################################
# Enrichment Network of Cancer Genes database
#############################################

# ncg.up <- enrichNCG(
#                       gene = up$entrez,
#                       pvalueCutoff = 0.05,
#                       pAdjustMethod = "BH",
#                       minGSSize = 10,
#                       maxGSSize = 500,
#                       qvalueCutoff = 0.2,
#                       readable = FALSE
#                     )
# 
# ncg.down <- enrichNCG(
#                       gene = down$entrez,
#                       pvalueCutoff = 0.05,
#                       pAdjustMethod = "BH",
#                       minGSSize = 10,
#                       maxGSSize = 500,
#                       qvalueCutoff = 0.2,
#                       readable = FALSE
#                     )
# 
# ##################################
# # VISUALIZATION OF ENRICHED NCGs
# ##################################
# 
# ncgupdot<-dotplot(ncg.up,showCategory = 30) + ggtitle(paste0("NCG ", set1,"_vs_",set2," (up in ",set2,")"))
# ggsave(paste0(fname,"/","dotplot_ncg_up_",set1,"_vs_",set2,".png"), plot=ncgupdot, width=10, height=10, units="in")
# #dev.off()
# 
# #pdf(paste0(fname,"/","dotplot_kegg_down_",set1,"_vs_",set2,".pdf"))
# ncgdowndot<-dotplot(ncg.down, showCategory=30) + ggtitle(paste0("NCG ", set1,"_vs_",set2," (down in ",set2,")"))
# ggsave(paste0(fname,"/","dotplot_ncg_down_",set1,"_vs_",set2,".png"), plot=ncgdowndot, width=10, height=10, units="in")
# 
# 
# ###################################
# # GSE NCG
# ###################################
# 
# gseNCG.up<-gseNCG(
#                   geneList=KEGGlist.up,
#                   exponent = 1,
#                   minGSSize = 10,
#                   maxGSSize = 500,
#                   pvalueCutoff = 0.05,
#                   pAdjustMethod = "BH",
#                   verbose = TRUE,
#                   seed = FALSE,
#                   by = "fgsea"
#                               )
# 
# 
# gseNCG.down<-gseNCG(
#                     geneList=KEGGlist.down,
#                     exponent = 1,
#                     minGSSize = 10,
#                     maxGSSize = 500,
#                     pvalueCutoff = 0.05,
#                     pAdjustMethod = "BH",
#                     verbose = TRUE,
#                     seed = FALSE,
#                     by = "fgsea"
#                                 )
# 
# ################################################
# # GSEA/NCG Visualization
# ################################################
# 
# #pdf(paste0(fname,"/","rigdeplot_gsea_kegg_up_",set1,"_vs_",set2,".pdf"))
# ncgrup<-ridgeplot(gseNCG.up) + ggtitle(paste0("GSEA/NCG ", set1,"_vs_",set2," (up in ",set2,")"))
# ggsave(paste0(fname,"/","rigdeplot_gsea_ncg_up_",set1,"_vs_",set2,".png"), plot=ncgrup, width=10, height=10, units="in")
# #dev.off()
# 
# #pdf(paste0(fname,"/","rigdeplot_gsea_kegg_down_",set1,"_vs_",set2,".pdf"))
# ncgrdown<-ridgeplot(gseNCG.down) + ggtitle(paste0("GSEA/NCG ", set1,"_vs_",set2," (down in ",set2,")"))
# ggsave(paste0(fname,"/","rigdeplot_gsea_ncg_down_",set1,"_vs_",set2,".pdf"), plot=ncgrdown, width=10, height=10, units="in")
options(warn = defaultW)
rm(list = ls()) # Clear workspace
gc(reset=TRUE)
