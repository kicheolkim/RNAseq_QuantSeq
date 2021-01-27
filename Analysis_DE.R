##### Differential expression analysis
invisible(lapply(c("tidyverse", "DESeq2", "pheatmap", "RColorBrewer"), library, character.only = TRUE))

setwd("/data/Analysis")

##### loading dataset #####
gtf_table <- read_tsv("/data/genome/gencode.vM25.pri.annotation.tsv")
gene_featureCounts <- read_csv("gene_featureCounts_output.csv")
sTable <- read_csv("sTable.csv")
CellTypes <- unique(sTable$CellType)

##### convert gene counts matrix
gene_counts <- as.data.frame(gene_featureCounts[,c(7:ncol(gene_featureCounts))])
row.names(gene_counts) <- gene_featureCounts$Geneid
gene_counts <- gene_counts[,sTable$SampleID]





##### DESeq2 #####
setwd("/data/Analysis")
dir.create(paste0("DE"), showWarnings = FALSE)
setwd(paste0("/data/Analysis/DE"))

cds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = sTable, design=~condition)
cds <- cds[ rowSums(counts(cds)) > 0, ]
cds <- DESeq(cds)
cds.res <- results(cds)
summary(cds.res)


de_result <- left_join(as_tibble(as.data.frame(cds.res) %>% rownames_to_column("gene_id")), gtf_table, by="gene_id")
write.csv(de_result, paste0("DESeq2_result.csv"))

tmp <- as.data.frame(cds.res)
result_table <- data.frame(condition = "mutant vs control",
                                 pval05_up = nrow(tmp %>% filter(pvalue<0.05, log2FoldChange>0)),
                                 pval05_down = nrow(tmp %>% filter(pvalue<0.05, log2FoldChange<0)),
                                 fdr_up = nrow(tmp %>% filter(padj<0.1, log2FoldChange>0)),
                                 fdr_down = nrow(tmp %>% filter(padj<0.1, log2FoldChange<0)))
rm(tmp)


cds.rlog <- rlog(cds, blind=FALSE)
plotPCA(cds.rlog, intgroup=c("condition")) + ggtitle("PCA for condition")
ggsave(paste0("PCA-DESeq2.png"), width=7, height=6, units="in")



##### Heatmap for significant miRNAs (pval < 0.01)
### Both (C vs NC)
res.tmp <- de_result
dds.tmp <- cds
sigList <- res.tmp %>% filter(pvalue < 0.01) %>% arrange(pvalue)

### transformation of counts from DESeq2 object
heatmap.expr <- assay(normTransform(dds.tmp))
selected.expr <- heatmap.expr[sigList$gene_id,]
selected.expr <- selected.expr - rowMeans(selected.expr)
row.names(selected.expr) <- sigList$gene_name

### make heatmap
df <- data.frame(Condition=dds.tmp$condition)
rownames(df) <- make.names(dds.tmp$SampleID, unique=TRUE)   #QC$sampleName   #
colnames(selected.expr) <- rownames(df)
anno_colors <- list(Condition=c(brewer.pal(5, "Set2")[1],brewer.pal(5, "Set2")[2]))
names(anno_colors[[1]]) <- unique(df$Condition)

while (!is.null(dev.list()))  dev.off()
pdf(file=paste0("Heatmap_DEsig_pval0.01.pdf"))  #, width=8, height=7)
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 3, fontsize_col = 4, scale = "row", #breaks=seq(-4,4, by=0.05),
         cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE,
         annotation_col=df, annotation_colors=anno_colors, fontsize=7)   #
dev.off()
rm(res.tmp, dds.tmp, heatmap.expr, selected.expr, df, anno_colors)


##### Heatmap for significant miRNAs (FDR < 0.1)
### Both (C vs NC)
res.tmp <- de_result
dds.tmp <- cds
sigList_fdr <- sigList %>% filter(padj < 0.1) %>% arrange(pvalue)

### transformation of counts from DESeq2 object
heatmap.expr <- assay(normTransform(dds.tmp))
selected.expr <- heatmap.expr[sigList_fdr$gene_id,]
selected.expr <- selected.expr - rowMeans(selected.expr)
row.names(selected.expr) <- sigList_fdr$gene_name

### make heatmap
df <- data.frame(Condition=dds.tmp$condition)
rownames(df) <- make.names(dds.tmp$SampleID, unique=TRUE)   #QC$sampleName   #
colnames(selected.expr) <- rownames(df)
anno_colors <- list(Condition=c(brewer.pal(5, "Set2")[1],brewer.pal(5, "Set2")[2]))
names(anno_colors[[1]]) <- unique(df$Condition)

while (!is.null(dev.list()))  dev.off()
pdf(file=paste0("Heatmap_DEsig_fdr0.1.pdf"))  #, width=8, height=7)
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 3, scale = "row", #breaks=seq(-4,4, by=0.05),
         cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE,
         annotation_col=df, annotation_colors=anno_colors, fontsize=7)   #
dev.off()
rm(res.tmp, dds.tmp, heatmap.expr, selected.expr, df, anno_colors)



##### Boxplots
dir.create(paste0("boxplots"), showWarnings = FALSE)
num <- 1
for (g in 1:50){
  geneCounts <- plotCounts(cds, gene = sigList$gene_id[g], intgroup = c("condition"), returnData = TRUE)
  ggplot(geneCounts, aes(x = condition, y = count, color = condition)) +
    geom_boxplot(outlier.shape=1) + scale_y_log10() + geom_point(position=position_jitterdodge(dodge.width=0.9)) +
    ggtitle(paste0(sigList$gene_name[g]," (",sigList$gene_id[g],")")) + xlab("") + ylab("DESeq2 norm. count - log scaled") +
    theme_bw() + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=10))
  ggsave(paste0("boxplots/Boxplot_",num,"-",sigList$gene_name[g],".png"), plot = last_plot(), device="png",
         scale = 1, width = 5, height = 5, units = "in")
  num <- num+1
  rm(geneCounts)
}
rm(num)

write_csv(result_table, "DE_summary.csv")
