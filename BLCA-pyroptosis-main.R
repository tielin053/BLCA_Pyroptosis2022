
# 00 rawdata --------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("00_rawdata")){
  dir.create("00_rawdata")
}
setwd("00_rawdata/")

library(TCGAbiolinks)
library(SummarizedExperiment)
grep("TCGA-BLCA",
     getGDCprojects()$project_id,
     value = T)

query <- GDCquery(project = "TCGA-BLCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
samplesDown <- getResults(query = query,
                          cols = c("cases"))
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
group_list <- data.frame(sample = c(dataSmNT, dataSmTP),
                         group = c(rep("Normal", 19), rep("Tumor", 412)))
write.table(group_list,
            file = "group_list.txt",
            quote = F,
            row.names = F)
group_list <- read.table("group_list.txt", header = T)
head(group_list)
dataSmNT <- group_list$sample[which(group_list$group=="Normal")]
dataSmTP <- group_list$sample[which(group_list$group=="Tumor")]

GDCdownload(query = query)
dataPrep1 <- GDCprepare(query = query,
                        save = TRUE,
                        save.filename = "blca_case1.RData")

dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                      cor.cut = 0.6,
                                      datatype = "unstranded")

dataPrep <- subset(dataPrep, select = c(dataSmNT, dataSmTP))
rownames(dataPrep) <- rowData(dataPrep1)[rownames(dataPrep),]$external_gene_name
dataNorm.blca <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                           geneInfo = geneInfo,
                                           method = "gcContent")
dataFilt.blca.final <- TCGAanalyze_Filtering(tabDF = dataNorm.blca,
                                             method = "quantile",
                                             qnt.cut = 0.25)
write.table(expr_count,
            file= "blca_count.xls",
            quote = F)
expr_count <- read.table("blca_count.xls", check.names = F)
head(expr_count[1:3,1:3])
clinical_data <- GDCquery_clinic(project = "TCGA-BLCA",
                                 type = "clinical")
write.table(clinical_data,file = "clinical_data.xls",sep="\t",row.names = T,col.names = T,quote = F)

library(stringr)
dataSmTP_sub <- unlist(lapply(dataSmTP, function(x) paste(str_split(x, "-")[[1]][0:3],
                                                          collapse ="-")))
tumor_clinical_sample <- intersect(dataSmTP_sub, clinical_data$submitter_id)
expr_tumor <- subset(expr_count, select = dataSmTP)
expr_tumor <- expr_tumor[, dataSmTP]
all(colnames(expr_tumor) == dataSmTP)
colnames(expr_tumor) <- dataSmTP_sub
expr_tumor <- subset(expr_tumor, select = tumor_clinical_sample) {
  library(TCGAbiolinks)
  library(SummarizedExperiment)

  query2 <- GDCquery(project = "TCGA-BLCA",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - FPKM")
  samplesDown <- getResults(query = query,
                            cols = c("cases"))

  GDCdownload(query = query2)
  dataPrep2 <- GDCprepare(query = query2,
                          save = TRUE,
                          save.filename = "blca_case2.RData")
  dataPrep3 <- TCGAanalyze_Preprocessing(object = dataPrep2,
                                        cor.cut = 0.6,
                                        datatype = "HTSeq - FPKM")
  dataPrep3 <- subset(dataPrep3, select = c(dataSmNT, dataSmTP))
  rownames(dataPrep3) <- rowData(dataPrep2)[rownames(dataPrep3),]$external_gene_name
  dataNorm.blca2 <- TCGAanalyze_Normalization(tabDF = dataPrep3,
                                             geneInfo = geneInfo,
                                             method = "geneLength")
  dataFilt.blca.final2 <- TCGAanalyze_Filtering(tabDF = dataPrep3,
                                               method = "quantile",
                                               qnt.cut = 0.25)
  expr_fpkm <- dataFilt.blca.final2
  write.table(expr_fpkm,
              file = "blca_fpkm.xls",
              quote = F)
  expr_fpkm <- read.table("../blca_fpkm2.xls", check.names = F)
  fpkm_tumor <- subset(expr_fpkm, select = dataSmTP)
  fpkm_tumor <- fpkm_tumor[, dataSmTP]
  all(colnames(fpkm_tumor) == dataSmTP)
  colnames(fpkm_tumor) <- dataSmTP_sub
  fpkm_tumor <- subset(fpkm_tumor, select = tumor_clinical_sample)
  write.table(fpkm_tumor,
              file = "blca_fpkm_tumor.xls",
              quote = F)
}

pyroptosis_gene <- read.csv("/Users/tielinwu/Projects/BLCA/data/pyroptosis.csv")
pyroptosis_gene <- pyroptosis_gene$gene
write.table(dataPrep,
            file = "blca_count_raw.xls",
            quote = F,
            sep = "\t")

write.table(dataPrep3,
            file = "blca_fpkm_raw.xls",
            quote = F,
            sep = "\t")
# 01 DEG ------------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./01_DEG")){
  dir.create("./01_DEG")
}
setwd("./01_DEG")

expr <- expr_count

rownames(group_list) <- group_list$sample
group_list$group <- factor(group_list$group, levels = c("Normal", "Tumor"))

all(rownames(group_list) %in% colnames(expr))
all(rownames(group_list) == colnames(expr))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = group_list,
                              design =~ group)
dds = dds[rownames(counts(dds)) > 1,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("group","Tumor","Normal"))
resOrdered <- res[order(res$pvalue),]
DEG <- as.data.frame(resOrdered)
DEG <- na.omit(DEG)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff ,'Up','Down'),'Not')
)


png(filename = "MA.png", height = 450, width = 600)
plotMA(res, ylim=c(-2,2))
dev.off()
pdf(file = "MA.pdf", height = 5)
plotMA(res, ylim=c(-2,2))
dev.off()

sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= 1)
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = F)

up_DEG <- DEG[which(DEG$change == "Up"),]
up_DEG <- up_DEG[order(-up_DEG$log2FoldChange, up_DEG$padj),]
down_DEG <- DEG[which(DEG$change == "Down"),]
down_DEG <- down_DEG[order(down_DEG$log2FoldChange, down_DEG$padj),]
sig_diff_write <- rbind(up_DEG, down_DEG)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = F)

data_repel <- rbind(up_DEG[1:10,], down_DEG[1:10,])
pyroptosis_sig <- sig_diff[which(rownames(sig_diff) %in% pyroptosis_gene),]
pyroptosis_sig <- pyroptosis_sig[order(pyroptosis_sig$log2FoldChange),]

library(VennDiagram)

venn.list2 = list("DEGs"=rownames(sig_diff),
                  "Pyroptosis"=pyroptosis_gene)
venn.diagram(venn.list2, 
             filename = 'DEG_pyroptosis_venn.png',
             imagetype = 'png',
             fill = c("blue", "red"), 
             alpha = 0.3, 
             fontface = "bold",
             # label.col = c("white", "white", "white","white",
             #               "white", "white", "white"),
             cat.col = c("blue", "red"),
             cat.cex = 1, 
             cat.pos = c(-20,40),
             cat.fontfamily = 'serif',
             col = c("blue", "red"),
             cex = 0.8, 
             ext.text = F,
             fontfamily = 'serif',
             lwd = 0.1,
             margin = 0.05,
             height = 2000, width = 2000)

# Volcano
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(patchwork)
library(ggrepel)
DEG$change <- factor(DEG$change, levels =  c("Up", "Not", "Down"))
{
  dss_con_plot <- ggplot(data = DEG,
                         aes(x = log2FoldChange,
                             y = -log10(padj), 
                             colour = change)) +
    scale_color_manual(values = c("red", "darkgray","blue")) +
    scale_x_continuous(breaks = c(-5,-3,-1,0,1,3,5)) +
    # scale_y_continuous(trans = "log1p") +
    geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
    ylim(0,55) +
    geom_label_repel(
      data = pyroptosis_sig[which(pyroptosis_sig$change=="Up"),],
      aes(label = rownames(pyroptosis_sig[which(pyroptosis_sig$change=="Up"),])),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black",
      show.legend = F,
      direction = "y",
      hjust = 0,
      xlim = c(12,20),
      ylim = c(1,50)
      )+
    geom_label_repel(
      data = pyroptosis_sig[which(pyroptosis_sig$change=="Down"),],
      aes(label = rownames(pyroptosis_sig[which(pyroptosis_sig$change=="Down"),])),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black",
      show.legend = F,
      direction = "y",
      hjust = 0,
      xlim = c(-5,-10),
      ylim = c(1,50)
    )+    theme_bw(base_size = 12, base_family = "Times") +
    geom_hline(yintercept = -log10(0.05),
               lty = 4,
               col = "darkgray",
               lwd = 0.6) +
    geom_vline(xintercept = c(-1,1),
               lty = 4,
               col = "darkgray",
               lwd = 0.6) +
    theme(
      legend.justification = c(1,1),
          legend.position = c(1,1),
          legend.background = element_rect(fill = "white", color = "black", size = 0.2),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face="bold",
                                     color="black",
                                     family = "Times",
                                     size=12),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face = "bold",
                                     color = "black",
                                     size = 12),
          axis.text.y = element_text(face = "bold",
                                     color = "black",
                                     size = 12),
          axis.title.x = element_text(face = "bold",
                                      color = "black",
                                      family = "Times",
                                      size = 18),
          axis.title.y = element_text(face = "bold",
                                      color = "black",
                                      family = "Times",
                                      size = 18)) +
    theme(text = element_text(family = "Times")) +
    labs(x = "log2 (Fold Change)",
         y = "-log10 (p-adjust)",
         title = "")
  
  dss_con_plot
  write_fig(dss_con_plot,
            file = "volcano_DEG_0.05_1.pdf",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(dss_con_plot,
            file = "volcano_DEG_0.05_1.png",
            width = 10,
            height = 8,
            devices = NULL,
            res = 300,
            show = F)
}
dss_con_plot <- ggplot(data = DEG,
                       aes(x = log2FoldChange,
                           y = -log10(padj), 
                           colour = change)) +
  scale_color_manual(values = c("red", "darkgray","blue")) +
  scale_x_continuous(breaks = c(-5,-3,-1,0,1,3,5),
                     limits = c(-6,10)) +
  # scale_y_continuous(trans = "log1p") +
  geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
  geom_label_repel(
    data = pyroptosis_sig[which(pyroptosis_sig$change=="Up"),],
    aes(label = rownames(pyroptosis_sig[which(pyroptosis_sig$change=="Up"),])),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black",
    show.legend = F,
    direction = "y",
    hjust = 0,
    xlim = c(8,15),
    ylim = c(1,50)
  )+
  geom_label_repel(
    data = pyroptosis_sig[which(pyroptosis_sig$change=="Down"),],
    aes(label = rownames(pyroptosis_sig[which(pyroptosis_sig$change=="Down"),])),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black",
    show.legend = F,
    direction = "y",
    hjust = 0,
    xlim = c(-5,-10),
    ylim = c(1,50)
  )+
  theme_bw(base_size = 12, base_family = "Times") +
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6) +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6) +
  theme(legend.justification = c(1,1),
    legend.position = c(1,1),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    family = "Times",
                                    size = 18),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    family = "Times",
                                    size = 18)) +
  theme(text = element_text(family = "Times")) +
  labs(x = "log2 (Fold Change)",
       y = "-log10 (p-adjust)",
       title = "")

dss_con_plot
write_fig(dss_con_plot,
          file = "volcano_DEG_0.05_1(2).pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(dss_con_plot,
          file = "volcano_DEG_0.05_1(2).png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)

# PCA
dds_vst <- vst(dds, blind = FALSE)
plotPCA(dds_vst, intgroup="group") +
  theme_bw() +
  geom_point(size = 5) +
  #scale_y_continuous(limits = c(-10, 10)) +
  ggtitle(label = "Principal Component Analysis(PCA)",
          subtitle = "")

# pheatmap
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))

mat <- assay(dds_vst[row.names(pyroptosis_sig)])
annotation_col <- data.frame(
  Group = factor(colData(dds_vst)$group),
  row.names = rownames(colData(dds_vst))
)
ann_colors <- list(
  Group = c(Normal="darkgreen", Tumor="darkorange"),
  Expression = c( Up="red",Down="blue")
)

mat <- mat[rownames(pyroptosis_sig[order(pyroptosis_sig$change, decreasing = T),]),]
gene_col <- data.frame(Expression=c(rep("Up",8), rep("Down", 1)))
rownames(gene_col) <- rownames(mat)

newnames <- lapply(
  rownames(mat),
  function(x) bquote(italic(.(x))))

pheatplot <- pheatmap(mat=scale(mat),
         annotation_col = annotation_col,
         annotation_row = gene_col,
         annotation_colors = ann_colors,
         fontsize = 12,
         labels_row = as.expression(newnames),
         show_colnames = FALSE,
         cluster_cols = F,
         cluster_rows = F,
         annotation_names_row = F,
         gaps_col = 19,
         gaps_row = 8)
ggsave(file="pheatmap.png", pheatplot)
ggsave(file="pheatmap.pdf", pheatplot)


# 02 gene KM --------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./02_gene_KM/")){
  dir.create("./02_gene_KM/")
}
setwd("./02_gene_KM/")

test <- subset(clinical_data, 
               select = c("days_to_last_follow_up","age_at_diagnosis","year_of_diagnosis","vital_status","age_at_index","days_to_birth","year_of_birth","days_to_death","year_of_death"))
View(test)
View(clinical_data)
clinical_data2 <- data.frame(sample = clinical_data$submitter_id,
                             OS = clinical_data$vital_status,
                             OS.time = ifelse(clinical_data$vital_status != "Dead",
                                              ifelse(!is.na(clinical_data$year_of_diagnosis),
                                                     (2019 - clinical_data$year_of_diagnosis) * 365,
                                                     clinical_data$days_to_last_follow_up),
                                              ifelse(!is.na(clinical_data$days_to_death),
                                                     clinical_data$days_to_death,
                                                     ifelse(!is.na(clinical_data$days_to_last_follow_up),
                                                            clinical_data$days_to_last_follow_up,
                                                            0))),
                             gender = clinical_data$gender,
                             age = clinical_data$age_at_index,
                             stage = clinical_data$ajcc_pathologic_stage,
                             T_stage = clinical_data$ajcc_pathologic_t,
                             N_stage = clinical_data$ajcc_pathologic_n,
                             M_stage = clinical_data$ajcc_pathologic_m)
clinical_data2$OS <- ifelse(clinical_data2$OS != "Dead",
                            0, 1)
View(clinical_data2)
km_data <- as.data.frame(t(pyroptosis_sig_count))
km_data$sample = rownames(km_data)
km_data <- merge(km_data, clinical_data2, by = "sample")

library(survival)
library(survminer)
library(Ipaper)
for (gene in rownames(pyroptosis_sig)){
  group <- paste(gene,"group",sep = "_")
  
  km_data[,group] <- ifelse(km_data[,gene]>median(km_data[,gene]), "High", "Low")
  
  gene_fit <- survfit(Surv(OS.time, OS)~get(group), data=km_data)
  assign(paste(gene,"surv", sep = "_")  ,ggsurvplot(gene_fit,
                          pval = TRUE, 
                          conf.int = F,
                          legend.labs=c("High Expr","Low Expr"),
                          legend.title="",
                          title=paste(gene, sep=""),
                          risk.table = TRUE, 
                          risk.table.col = "strata", 
                          linetype = "strata", 
                          # surv.median.line = "hv",
                          risk.table.fontsize = 5,
                          ggtheme = theme_bw(), 
                          palette = c("#A73030FF", "#0073C2FF"),
                          font.family = "Times"
                          )
         )
}


for (i in c("IL6","GSDMA","CASP6","ZBP1","CHMP4C","IL1A","AIM2","NLRP2","NLRP7")){
  s_surv <- get(paste(i,"surv",sep="_"))
  
  s_surv$table <- s_surv$table + 
    labs(x = "Overall Survival (days)") +
    theme(panel.grid = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          text = element_text(size = 15, face = "bold"))
  s_surv$plot <- s_surv$plot + 
    labs(x = "Overall Survival (days)") +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(1,1),
          legend.position = c(1,1),
          legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(fill = "white", color = "black", size = 0.2),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 13, face = "bold"),
          text = element_text(size = 20, face = "bold"))
  
  write_fig(s_surv,
            file = paste(i,"km.png",sep="_"),
            width = 5,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(s_surv,
            file = paste(i,"km.pdf",sep="_"),
            width = 5,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}

library(patchwork)
surv_all <- IL6_surv$plot + GSDMA_surv$plot + CASP6_surv$plot + ZBP1_surv$plot + 
CHMP4C_surv$plot + IL1A_surv$plot + AIM2_surv$plot + NLRP2_surv$plot + 
  NLRP7_surv$plot +plot_layout(ncol =3) &
  labs(x = "Overall Survival (days)") &
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10, face = "bold"),
        text = element_text(size = 20, face = "bold"))
surv_all
write_fig(surv_all,
          file = "all_gene_km.pdf",
          width = 10,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
write_fig(surv_all,
          file = "all_gene_km.png",
          width = 10,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)



# 03 cluster --------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./03_cluster")){
  dir.create("./03_cluster")
}
setwd("./03_cluster")

library(ConsensusClusterPlus)

sig_km_gene <- c("CHMP4C", "GSDMA")
pyroptosis_sig_km_expr <- expr_tumor[which(rownames(expr_tumor) %in% sig_km_gene),]
d <- sweep(pyroptosis_sig_km_expr, 1, apply(pyroptosis_sig_km_expr, 1, median, na.rm=T))
cluster_result <- ConsensusClusterPlus(as.matrix(d),
                                       maxK = 8,
                                       reps = 1000,
                                       pItem = 0.8,
                                       pFeature = 1,
                                       title = "ConsensusClusterPlus",
                                       clusterAlg = "hc",
                                       distance = "spearman",
                                       seed = 123,
                                       plot = "pdf")
new_group <- cluster_result[[2]][["consensusClass"]]
cluster1 <- names(new_group[new_group == 1])
cluster2 <- names(new_group[new_group == 2])
write.table(cluster_list,
            file = "cluster_list.txt",
            quote = F,
            sep = "\t",
            row.names = F)
# 04 tSNE -----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./04_tSNE")){
  dir.create("./04_tSNE")
}
setwd("./04_tSNE")

{
  # t-SNE 
  d = pyroptosis_sig_km_expr
  require(Rtsne)
  dtsne = sweep(pyroptosis_sig_km_expr, 1, apply(pyroptosis_sig_km_expr, 1, median, na.rm=T))
  set.seed(18)
  tsne_out <- Rtsne(t(dtsne),
                dims = 2, 
                perplexity=10, 
                verbose=TRUE, 
                max_iter = 100,
                check_duplicates = F,
                theta = 0)
  tsne_data_group <- as.data.frame(t(dtsne))
  tsne_data_group$Cluster <- ifelse(rownames(tsne_data_group) %in% cluster1,
                                  "Cluster1", "Cluster2")
  table(tsne_data_group$Cluster)
    library(ggplot2)
    tsne_res <- as.data.frame(tsne_out$Y)
    colnames(tsne_res) <- c("tSNE1","tSNE2")
    save(tsne_res, file = "tsne_res_seed18.RData")
    tsne_plot <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=tsne_data_group$Cluster)) + 
      geom_point(size = 1) + theme_bw() + 
      geom_hline(yintercept = 0,lty=2,col="grey") + 
      geom_vline(xintercept = 0,lty=2,col="grey") +
      stat_ellipse(data = tsne_res,
                   geom="polygon",
                   aes(x=tSNE1,y=tSNE2,
                       group=tsne_data_group$Cluster,
                       fill=tsne_data_group$Cluster),
                   alpha = 0.2,
                   lty="dashed",
                   # color = "black",
                   key_glyph = "blank") +
      guides(fill="none")+
      scale_color_discrete(labels = c("Cluster1(n=209)", "Cluster2(n=199)")) +
      theme(legend.position = "top",
            legend.text = element_text(size = 15, family = "Times"),
            axis.text = element_text(size = 12, family = "Times"),
            axis.title = element_text(size = 16, family = "Times"),
            panel.grid = element_blank()) + 
      labs(title = "",color="")
    tsne_plot
    ggsave(filename = "tsne_plot.png", tsne_plot, width = 7, height = 5)
    ggsave(filename = "tsne_plot.pdf", tsne_plot, width = 7, height = 5)

}

write.table(tsne_data_group,
            file = "tsne_data_group.xls",
            sep = "\t",
            quote = F)
write.table(tsne_res,
            file = "tsne_result.xls",
            sep = "\t",
            quote =F)


# 05 cluster survival analysis -----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./05_cluster_KM")){
  dir.create("./05_cluster_KM")
}
setwd("./05_cluster_KM")

km_data$Cluster <- factor(ifelse(km_data$sample %in% cluster1,
                          "Cluster1", "Cluster2"))
write.table(km_data,
            file = "km_data.xls",
            sep = "\t",
            quote = F)
library(survival)
library(survminer)
library(Ipaper)
kmfit <- survfit(Surv(OS.time, OS)~Cluster, data=km_data)
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE, 
                     # pval.method = T,
                     conf.int = F,
                     legend.labs=c("Cluster1","Cluster2"),
                     legend.title="Cluster",
                     # title="Cluster KM",
                     risk.table = TRUE, 
                     risk.table.col = "strata",
                     linetype = "strata",
                     # surv.median.line = "hv",
                     risk.table.fontsize = 5,
                     ggtheme = theme_bw(),
                     palette = c("#f8776e", "#00c4c6"),
                     font.family = "Times"
                     )
s_surv

s_surv$table <- s_surv$table + 
  labs(x = "Overall Survival (days)") +
  theme(panel.grid = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 15, face = "bold"))
s_surv$plot <- s_surv$plot + 
  labs(x = "Overall Survival (days)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 13, face = "bold"),
        text = element_text(size = 20, face = "bold"))

s_surv

write_fig(s_surv,
          file = "cluster_km.pdf",
          width = 5,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(s_surv,
          file = "cluster_km.png",
          width = 5,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)

# 06 Clinical analysis -----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./06_clinical_index")){
  dir.create("./06_clinical_index")
}
setwd("./06_clinical_index")

clinical_index <- clinical_data2
clinical_index <- merge(clinical_index,
                        cluster_list,
                        by = "sample")
clinical_index$T_stage <- gsub("^T0$", "TX", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^TX$", "TX", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T1[a-c]*$", "T1", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T2[a-c]*$", "T2", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T3[a-c]*$", "T3", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T4[a-c]*$", "T4", clinical_index$T_stage)

clinical_index$N_stage <- gsub("^NX$", "NX", clinical_index$N_stage)
clinical_index$M_stage <- gsub("^MX$", "MX", clinical_index$M_stage)

clinical_index$age_stage <- ifelse(clinical_index$age <= 60,
                                   "30~60", ifelse(clinical_index$age <= 70 & clinical_index$age > 60,
                                                   "60~70", ">70"))

# pheatmap
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))

mat <- expr_tumor[rownames(data_repel2),c(cluster1,cluster2)]
dim(mat)
rownames(clinical_index) <- clinical_index$sample
clinical_index <- clinical_index[rownames(colData(dds_vst)),]
all(rownames(clinical_index) == rownames(colData(dds_vst)))
annotation_col <- subset(clinical_index, 
                         select = c("cluster","T_stage","N_stage","M_stage",
                                    "gender","age_stage"))
annotation_col$age_stage <- factor(annotation_col$age_stage)
names(annotation_col) = c("Cluster","T**","N","M**",
                          "Gender","Age")
annotation_col <- annotation_col[,c("Age", "Gender","M**","N","T**", "Cluster")]
ann_colors <- list(
  Cluster = c(Cluster1="#f8776e", Cluster2="#00c4c6"),
  "T**" = c(T1="#FEF0D9", T2="#FDCC8A", T3="#FC8D59", T4="#D7301F", TX="grey"),
  N = c(N0="#EDF8FB", N1="#B3CDE3", N2="#8C96C6", N3="#88419D", NX="grey"),
  "M**" = c(M0="#FBB4B9", M1="#AE017E", MX="grey"),
  Gender = c(female="#c71585", male="#1e90ff"),
  Age = c("30~60"="#33A02C", "60~70" = "#1F78B4", ">70"="#A6CEE3")
)


newnames <- lapply(
  rownames(mat),
  function(x) bquote(italic(.(x))))

pheatplot <- pheatmap::pheatmap(log(mat+1),
                                # scale = "row",
                                annotation_col = annotation_col,
                                annotation_colors = ann_colors,
                                fontsize = 12,
                                # labels_row = as.expression(newnames),
                                show_colnames = FALSE,
                                cluster_cols = F,
                                cluster_rows = T,
                                na_col = "grey",
                                annotation_names_row = F,
                                show_rownames = F
                                )
ggsave(file="pheatmap.png", pheatplot, height = 10)
ggsave(file="pheatmap.pdf", pheatplot, height = 10)


{
clinical_index$stage <- gsub("^Stage I$", "Stage I/II", clinical_index$stage)
clinical_index$stage <- gsub("^Stage II$", "Stage I/II", clinical_index$stage)
clinical_index$stage <- gsub("^Stage III$", "Stage III/IV", clinical_index$stage)
clinical_index$stage <- gsub("^Stage IV$", "Stage III/IV", clinical_index$stage)

clinical_index$T_stage <- gsub("T0", NA, clinical_index$T_stage)
clinical_index$T_stage <- gsub("TX", NA, clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T1$", "T1/T2", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T2[a-c]*$", "T1/T2", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T3[a-c]*$", "T3/T4", clinical_index$T_stage)
clinical_index$T_stage <- gsub("^T4[a-c]*$", "T3/T4", clinical_index$T_stage)

clinical_index$N_stage <- gsub("NX", NA, clinical_index$N_stage)
clinical_index$N_stage <- gsub("^N0$", "N0/1", clinical_index$N_stage)
clinical_index$N_stage <- gsub("^N1$", "N0/1", clinical_index$N_stage)
clinical_index$N_stage <- gsub("^N2$", "N2/3", clinical_index$N_stage)
clinical_index$N_stage <- gsub("^N3$", "N2/3", clinical_index$N_stage)


library(ggpubr)
library(Ipaper)
library(ggthemes)
my_comparisons <- list(c("Cluster1", "Cluster2"))

male_clinical <- clinical_index[which(clinical_index$gender=="male"),]
famale_clinical <- clinical_index[which(clinical_index$gender=="female"),]
gender_data <- read.table("gender_data.xls", header = T)
p_gender <- ggplot(gender_data,
                   aes(x=Gender,
                       y=Count,
                       fill=Cluster)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#f8776e", "#00c4c6")) +
  theme_base() +
  labs(title = "", x = "") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        text = element_text(size = 20, face = "bold"))

p_gender
write_fig(p_gender,
          file = "gender.pdf",
          width = 4,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p_gender,
          file = "gender.png",
          width = 4,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)

library(ggplot2)
library(ggridges)
library(RColorBrewer)

p_age <- ggplot(clinical_index,aes(x = age,
                          y = cluster,
                          fill = cluster)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height=0.00, size = 0.3, alpha = 0.3) +
  # scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)) +
  scale_fill_manual(values = c("#f8776e", "#00c4c6")) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        text = element_text(size = 20, face = "bold"))
p_age

write_fig(p_age,
          file = "age.pdf",
          width = 8,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p_age,
          file = "age.png",
          width = 8,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
#T
{
  t1_t2 <- clinical_index[which(clinical_index$T_stage=="T1/T2"),]
  t3_t4 <- clinical_index[which(clinical_index$T_stage=="T3/T4"),]
  t_stage_data <- read.table("t_stage_data.xls", header = T)
  p_t_stage <- ggplot(t_stage_data,
                      aes(x=T_stage,
                          y=Count,
                          fill=Cluster)) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6")) +
    theme_base() +
    labs(title = "", x = "") +
    theme(legend.title = element_blank(),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 13, face = "bold"),
          text = element_text(size = 20, face = "bold"))
  
  p_t_stage
  write_fig(p_t_stage,
            file = "t_stage.pdf",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p_t_stage,
            file = "t_stage.png",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}
#M
{
  m0 <- clinical_index[which(clinical_index$M_stage=="M0"),]
  m1 <- clinical_index[which(clinical_index$M_stage=="M1"),]
  mx <- clinical_index[which(clinical_index$M_stage=="MX"),]
  m_stage_data <- read.table("m_stage_data.xls", header = T)
  p_m_stage <- ggplot(m_stage_data,
                      aes(x=M_stage,
                          y=Count,
                          fill=Cluster)) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6")) +
    theme_base() +
    labs(title = "", x = "") +
    theme(legend.title = element_blank(),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 13, face = "bold"),
          text = element_text(size = 20, face = "bold"))
  
  p_m_stage
  write_fig(p_m_stage,
            file = "m_stage.pdf",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p_m_stage,
            file = "m_stage.png",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}
#N
{
  n0_n1 <- clinical_index[which(clinical_index$N_stage=="N0/1"),]
  n2_n3 <- clinical_index[which(clinical_index$N_stage=="N2/3"),]
  n_stage_data <- read.table("n_stage_data.xls", header = T)
  p_n_stage <- ggplot(n_stage_data,
                      aes(x=N_stage,
                          y=Count,
                          fill=Cluster)) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6")) +
    theme_base() +
    labs(title = "", x = "") +
    theme(legend.title = element_blank(),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 13, face = "bold"),
          text = element_text(size = 20, face = "bold"))
  
  p_n_stage
  write_fig(p_n_stage,
            file = "n_stage.pdf",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p_n_stage,
            file = "n_stage.png",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}
#stage
{
  stage1_2 <- clinical_index[which(clinical_index$stage =="Stage I/II"),]
  stage3_4 <- clinical_index[which(clinical_index$stage =="Stage III/IV"),]
  stage_data <- read.table("stage_data.xls", header = T, sep = "\t")
  p_stage <- ggplot(stage_data,
                      aes(x=Stage,
                          y=Count,
                          fill=Cluster)) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6")) +
    theme_base() +
    labs(title = "", x = "") +
    theme(legend.title = element_blank(),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 13, face = "bold"),
          text = element_text(size = 20, face = "bold"))
  
  p_stage
  write_fig(p_stage,
            file = "stage.pdf",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p_stage,
            file = "stage.png",
            width = 4,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}

library(patchwork)
p_all_clinical <- p_age /
  (p_gender + p_stage & theme(legend.position = "none", axis.title.y = element_blank()))  /
  (p_t_stage + p_m_stage + p_n_stage& theme(legend.position="none", axis.title.y = element_blank())) 
p_all_clinical

write_fig(p_all_clinical,
          file = "all_clinical_index.pdf",
          width = 9,
          height = 12,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p_all_clinical,
          file = "all_clinical_index.png",
          width = 9,
          height = 12,
          devices = NULL,
          res = 600,
          show = F)

}
# 07 TIICs ----------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./07_TIICs")){
  dir.create("./07_TIICs")
}
setwd("./07_TIICs")

fpkm_tumor2 <- dplyr::filter(as.data.frame(fpkm_tumor), !duplicated(rownames(as.data.frame(fpkm_tumor))))

write.table(fpkm_tumor2,
            file = "expr_fpkm_tumor.xls",
            sep = "\t",
            quote = F)
# CIBERSORT
{
  source("/Users/tielinwu/Projects/BLCA/data/CIBERSORT.R")
  result <- CIBERSORT('/Users/tielinwu/Projects/BLCA/data/LM22.txt',
                      'expr_fpkm_tumor.xls', 
                      perm = 1000, 
                      QN = F)
  cibersort_raw <- read.table("CIBERSORT-Results.txt",
                              header = T,
                              sep = "\t",
                              row.names = 1,
                              check.names = F)
  cibersort_result <- t(cibersort_raw[,-c(23,24,25)])
}

{
  tiics_result <- cibersort_result
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)){
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, cluster1],
                                         tiics_result[i, cluster2])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, cluster1]) - 
      mean(tiics_result[i, cluster2])
  }
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <- data.frame(log2FoldChange, 
                       pvalue, 
                       padj,
                       row.names = rownames(tiics_result))
  cluster1_res <- signif(apply(tiics_result[rownames(rTable), cluster1], 
                       1,
                       mean), 4)
  cluster2_res <- signif(apply(tiics_result[rownames(rTable), cluster2], 
                      1, 
                      mean), 4)
  rTable <- data.frame(cluster1_res, 
                       cluster2_res,
                       rTable[, c("padj", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(rTable$padj < 0.05,
                       ifelse(rTable$padj < 0.01, 
                              ifelse(rTable$padj < 0.001,
                                     ifelse(rTable$padj < 0.0001,
                                            paste(rTable$immune_cell, "****",  sep = ""),
                                            paste(rTable$immune_cell, "***", sep = "")),
                                     paste(rTable$immune_cell, "**", sep = "")),
                              paste(rTable$immune_cell, "*",  sep = "")), 
                       rTable$immune_cell)
  
  write.table(rTable,
              file = "tiics_wilcox_test.xls",
              quote = F,
              row.names = F)
}

library(ggradar)
cibersort_result2 <- cibersort_result
cibersort_result2 <- cibersort_result2[rTable$immune_cell,]
rownames(cibersort_result2) <- rTable$sig
cibersort_result_t <- as.data.frame(t(cibersort_result2))
cibersort_result_t <- cbind(group = ifelse(rownames(cibersort_result_t) %in% cluster1,
                                           "Cluster1", "Cluster2"),
                            cibersort_result_t)     
cibersort_result_t <- aggregate(cibersort_result_t[,2:ncol(cibersort_result_t)],
                                by = list(group = cibersort_result_t$group), mean)

tiics_radar <- ggradar(cibersort_result_t,
                       values.radar = c("0", "0.07", "0.14"),
                       grid.min = 0,
                       grid.mid = 0.07,
                       grid.max = 0.145,
                       group.line.width = 0.8, 
                       group.point.size = 1.5,
                       plot.extent.x.sf = 2,
                       plot.extent.y.sf = 1.25,
                       axis.label.size = 4.5,
                       grid.label.size = 5,
                       background.circle.colour = "white",
                       gridline.mid.colour = "grey",
                       gridline.min.colour = "grey",
                       axis.line.colour = "grey",
                       group.colours = c("#f8776e", "#00c4c6"),
                       plot.title = "TIICs Distribution") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = c(0.8,1))
tiics_radar
write_fig(tiics_radar,
          file = "tiics_radar.pdf",
          width = 14,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(tiics_radar,
          file = "tiics_radar.png",
          width = 16,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)
{
  cibersort_result <- cibersort_result[,c(cluster1,cluster2)]
  pdf("tiics_barplot.pdf",height=10,width=18)
  col=rainbow(nrow(cibersort_result),s=0.7,v=0.7)
  par(las=1,mar=c(8,5,4,20),mgp=c(3,0.1,0),cex.axis=1.5)
  a1=barplot(cibersort_result,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
  a2=axis(2,tick=F,labels=F)
  axis(2,a2,paste0(a2*100,"%"))
  par(srt=0,xpd=T)
  rect(xleft = a1[1], ybottom = -0.01, xright = a1[209], ytop = -0.06,col="#f8776e")
  text(a1[209]/2,-0.035,"Cluster1",cex=2)
  rect(xleft = a1[209], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="#00c4c6")
  text((a1[length(a1)]+a1[209])/2,-0.035,"Cluster2",cex=2)
  ytick2 = cumsum(cibersort_result[,ncol(cibersort_result)])
  ytick1 = c(0,ytick2[-length(ytick2)])
  legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(cibersort_result),col=col,pch=15,bty="n",cex=1.6)
  dev.off()
  
  png("tiics_barplot.png",height=10,width=18, units = "in", res = 300)
  col=rainbow(nrow(cibersort_result),s=0.7,v=0.7)
  par(las=1,mar=c(8,5,4,20),mgp=c(3,0.1,0),cex.axis=1.5)
  a1=barplot(cibersort_result,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
  a2=axis(2,tick=F,labels=F)
  axis(2,a2,paste0(a2*100,"%"))
  par(srt=0,xpd=T)
  rect(xleft = a1[1], ybottom = -0.01, xright = a1[209], ytop = -0.06,col="#f8776e")
  text(a1[209]/2,-0.035,"Cluster1",cex=2)
  rect(xleft = a1[209], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="#00c4c6")
  text((a1[length(a1)]+a1[209])/2,-0.035,"Cluster2",cex=2)
  ytick2 = cumsum(cibersort_result[,ncol(cibersort_result)])
  ytick1 = c(0,ytick2[-length(ytick2)])
  legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(cibersort_result),col=col,pch=15,bty="n",cex=1.6)
  dev.off()
}

# 08 ESTIMATE analysis ------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./08_ESTIMATE")){
  dir.create("./08_ESTIMATE")
}
setwd("./08_ESTIMATE")

library(estimate)

expr_train <- log2(fpkm_tumor2+1)
filterCommonGenes(input.f = './expr_sample408_log2.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
estimateScore('expr_train.gct', 'train_purity.gct', platform="illumina")
es_score <- read.table('train_purity.gct', skip = 2, header = T, check.names = F)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME

write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)

violin_dat <- data.frame(t(immu_score))
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% gsub("-",".",cluster1),
                           "Cluster1", "Cluster2")

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + # violin plot
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6"), name = "") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ 
    ylim(-2600,2200) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(legend.position = "none",
          axis.text.y=element_text(family="Times",size=12,face="bold"), 
          axis.text.x = element_text(family="Times",size = 15,face="bold"),
          axis.title.y=element_text(family="Times",size = 15,face="bold"), 
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold", size = 10))+
    labs(title="Violin plot of StromalScore", x="", y="StromalScore")
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + 
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6"), name = "Cluster") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ 
    ylim(-2000, 3200) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(legend.position = "none",
          axis.text.y=element_text(family="Times",size=12,face="bold"), 
          axis.text.x = element_text(family="Times",size = 15,face="bold"),
          axis.title.y=element_text(family="Times",size = 15,face="bold"), 
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold", size = 10))+
    labs(title="Violin plot of ImmuneScore", x="", y="ImmuneScore")
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = c("#f8776e", "#00c4c6"), name = "Cluster") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ 
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(-4500, 5500) +
    theme(legend.position = "none",
          axis.text.y=element_text(family="Times",size=12,face="bold"),
          axis.text.x = element_text(family="Times",size = 15,face="bold"),
          axis.title.y=element_text(family="Times",size = 15,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold", size = 10))+
    labs(title="Violin plot of ESTIMATEScore", x="", y="ESTIMATEScore")
  p3
  
  p4 <- cowplot::plot_grid(p1,p2,p3,
                           nrow = 1, 
                           align = 'h', 
                           vjust = -0.3)
  p4
}
{
  write_fig(p4,
            file = "estimate_plot.pdf",
            width = 20,
            height = 6,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p4,
            file = "estimate_plot.png",
            width = 20,
            height = 6,
            devices = NULL,
            res = 600,
            show = F)
  
  
  write_fig(p1,
            file = "violin_StromalScore.pdf",
            width = 6,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p1,
            file = "violin_StromalScore.png",
            width = 6,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  
  write_fig(p2,
            file = "violin_ImmuneScore.pdf",
            width = 6,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p2,
            file = "violin_ImmuneScore.png",
            width = 6,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  
  write_fig(p3,
            file = "violin_ESTIMATEScore.pdf",
            width = 6,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(p3,
            file = "violin_ESTIMATEScore.png",
            width = 6,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}


# 09 TMB ------------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./09_TMB")){
  dir.create("./09_TMB")
}
setwd("./09_TMB")

library(TCGAmutations)
library(stringr)

tmp=as.data.frame(tcga_available())
dt_blca <- TCGAmutations::tcga_load(study = "BLCA")
dt_blca2 <- dt_blca@data
dt1 <- as.data.frame( table(dt_blca2$Tumor_Sample_Barcode))
names(dt1) <- c('sample', 'Freq')
dt1$tmb <- dt1$Freq/38
dt1$sample <- str_sub(dt1$sample, 1, 15)
names(dt1)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
violin_dat <- dt1
violin_dat$group <- ifelse(violin_dat$sample %in% gsub("$", "-01", cluster1),
                            "Cluster1", "Cluster2") 
head(violin_dat)
my_comprasons <- c("Cluster1", "Cluster2")
violin_plot <- ggplot(violin_dat, aes(x=group,
                                      y=tmb,
                                      fill=group)) +
  geom_violin() +
  geom_boxplot(fill = "white", width = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = c("#f8776e", "#00c4c6"), name = "Cluster") +
  stat_compare_means(aes(group = group),method = 'wilcox')+ 
  theme_bw()+ 
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank()) + 
  labs(y="TMB", x="") +
  theme(legend.position = "none",
    axis.text.x=element_text(family="Times",size=12,face="bold"), 
        axis.text.y=element_text(family="Times",size=12,face="bold"), 
        axis.title.y=element_text(family="Times",size = 15,face="bold"), 
        axis.title.x=element_text(family="Times",size = 10,face="bold"))
violin_plot

write_fig(violin_plot,
          file = "violin_TMB.pdf",
          width = 4,
          height = 4,
          devices = NULL,
          res = 600,
          show = F)
write_fig(violin_plot,
          file = "violin_TMB.png",
          width = 4,
          height = 4,
          devices = NULL,
          res = 600,
          show = F)


write.table(dt1, file = 'TMB.csv',
            quote = F,
            row.names = F)

# 10 checkpoint -------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./10_checkpoint")){
  dir.create("./10_checkpoint")
}
setwd("./10_checkpoint")

checkpoint <- read.table("../data/checkpoint.txt",
                         header = F)
checkpoint <- checkpoint$V1
common_checkpoint <- intersect(checkpoint, rownames(fpkm_tumor2))
length(common_checkpoint)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

common_checkpoint_expr <- fpkm_tumor2[which(rownames(fpkm_tumor2) %in% common_checkpoint),]
dim(common_checkpoint_expr)
common_checkpoint_expr <- log2(common_checkpoint_expr + 1)
common_checkpoint_expr$gene <- rownames(common_checkpoint_expr)

violin_dat <- gather(common_checkpoint_expr, key=sample, value='log2(FPKM+1)', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% cluster1,
                           "Cluster1", "Cluster2") 
head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=`log2(FPKM+1)`,
                                      fill=group)) +
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#f8776e", "#00c4c6"), name = "Cluster")+
  labs(title="", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot
write_fig(violin_plot,
          file = "checkpoint_diff.pdf",
          width = 15,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(violin_plot,
          file = "checkpoint_diff.png",
          width = 15,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)


# 11 cluster DEG ----------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./11_Cluster_DEG")){
  dir.create("./11_Cluster_DEG")
}
setwd("./11_Cluster_DEG")

cluster <- read.table("../03_cluster/cluster_list.txt", header = T)
cluster1 <- cluster$sample[which(cluster$cluster=="Cluster1")]
cluster2 <- cluster$sample[which(cluster$cluster=="Cluster2")]
cluster_sample <- c(cluster1, cluster2)
cluster_list <- data.frame(sample=cluster_sample)
rownames(cluster_list) <- cluster_list$sample
cluster_list$cluster <- ifelse(rownames(cluster_list) %in% cluster1, "Cluster1", "Cluster2") 
cluster_list$cluster <- factor(cluster_list$cluster, levels = c("Cluster1", "Cluster2"))

expr_tumor2 <- expr_tumor[, rownames(cluster_list)]
all(rownames(cluster_list) %in% colnames(expr_tumor2))
all(rownames(cluster_list) == colnames(expr_tumor2))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = expr_tumor2,
                              colData = cluster_list,
                              design =~ cluster)
dds = dds[rownames(counts(dds)) > 1,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("cluster","Cluster1","Cluster2"))
resOrdered <- res[order(res$pvalue),]
DEG2 <- as.data.frame(resOrdered)
DEG2 <- na.omit(DEG2)
logFC_cutoff <- 1
DEG2$change = as.factor(
  ifelse(DEG2$padj < 0.05 & abs(DEG2$log2FoldChange) > logFC_cutoff,
         ifelse(DEG2$log2FoldChange > logFC_cutoff ,'Up','Down'),'Not')
)

png(filename = "MA.png", height = 450, width = 600)
plotMA(res, ylim=c(-2,2))
dev.off()
pdf(file = "MA.pdf", height = 5)
plotMA(res, ylim=c(-2,2))
dev.off()


sig_diff2 <- subset(DEG2,
                   DEG2$padj < 0.05 & abs(DEG2$log2FoldChange) >= 1)
cluster_diff_gene <- rownames(sig_diff2)

DEG_write2 <- cbind(GeneSymbol=rownames(DEG2), DEG2)
write.table(DEG_write2, file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = F)

up_DEG2 <- DEG2[which(DEG2$change == "Up"),]
up_DEG2 <- up_DEG2[order(-up_DEG2$log2FoldChange, up_DEG2$padj),]
down_DEG2 <- DEG2[which(DEG2$change == "Down"),]
down_DEG2 <- down_DEG2[order(down_DEG2$log2FoldChange, down_DEG2$padj),]
sig_diff_write2 <- rbind(up_DEG2, down_DEG2)
sig_diff_write2 <- cbind(GeneSymbol=rownames(sig_diff2), sig_diff2)
write.table(sig_diff_write2, file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = F)

data_repel2 <- rbind(up_DEG2[1:20,], down_DEG2[1:20,])
diff_cluster_pygene <- intersect(pyroptosis_gene, rownames(sig_diff2))
diff_cluster_pygene_res <- sig_diff2[diff_cluster_pygene,]
diff_cluster_pygene_res <- diff_cluster_pygene_res[order(diff_cluster_pygene_res$log2FoldChange, decreasing = T),]

library(VennDiagram)
venn.list2 = list("DEGs"=rownames(sig_diff2),
                  "Pyroptosis"=pyroptosis_gene)
venn.diagram(venn.list2, 
             filename = 'DEG_pyroptosis_venn.png',
             imagetype = 'png',
             fill = c("blue", "red"), 
             alpha = 0.3, 
             fontface = "bold",
             # label.col = c("white", "white", "white","white",
             #               "white", "white", "white"),
             cat.col = c("blue", "red"),
             cat.cex = 1, 
             cat.pos = c(-20,30),
             cat.fontfamily = 'serif',
             col = c("blue", "red"),
             cex = 0.8, 
             ext.text = F,
             fontfamily = 'serif',
             lwd = 0.1,
             margin = 0.05,
             height = 2000, width = 2000)

# Volcano
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(patchwork)
library(ggrepel)
DEG2$change <- factor(DEG2$change, levels =  c("Up", "Not", "Down"))
{
  dss_con_plot2 <- ggplot(data = DEG2,
                         aes(x = log2FoldChange,
                             y = -log10(padj), 
                             colour = change)) +
    scale_color_manual(values = c("red", "darkgray","blue")) +
    scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) +
    # scale_y_continuous(trans = "log1p") +
    geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
    theme_bw(base_size = 12, base_family = "Times") +
    geom_hline(yintercept = -log10(0.05),
               lty = 4,
               col = "darkgray",
               lwd = 0.6) +
    geom_vline(xintercept = c(-1,1),
               lty = 4,
               col = "darkgray",
               lwd = 0.6) +
    theme(legend.justification = c(1,1),
          legend.position = c(1,1),
          legend.background = element_rect(fill = "white", color = "black", size = 0.2),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face="bold",
                                     color="black",
                                     family = "Times",
                                     size=12),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face = "bold",
                                     color = "black",
                                     size = 12),
          axis.text.y = element_text(face = "bold",
                                     color = "black",
                                     size = 12),
          axis.title.x = element_text(face = "bold",
                                      color = "black",
                                      size = 18),
          axis.title.y = element_text(face = "bold",
                                      color = "black",
                                      size = 18)) +
    labs(x = "log2 (Fold Change)",
         y = "-log10 (p-adjust)",
         title = "")
  
  dss_con_plot2
  write_fig(dss_con_plot2,
            file = "volcano_DEG_0.05_1.pdf",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(dss_con_plot2,
            file = "volcano_DEG_0.05_1.png",
            width = 10,
            height = 8,
            devices = NULL,
            res = 300,
            show = F)
}


# PCA
dds_vst <- vst(dds, blind = FALSE)
plotPCA(dds_vst, intgroup="cluster") +
  theme_bw() +
  geom_point(size = 5) +
  #scale_y_continuous(limits = c(-10, 10)) +
  ggtitle(label = "Principal Component Analysis(PCA)",
          subtitle = "")


# pheatmap
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))

mat <- assay(dds_vst[row.names(diff_cluster_pygene_res)])
annotation_col <- data.frame(
  Cluster = factor(colData(dds_vst)$cluster),
  row.names = rownames(colData(dds_vst))
)
ann_colors <- list(
  Cluster = c(Cluster1="#f8776e", Cluster2="#00c4c6"),
  Expression = c(Up="red", Down="blue")
)

mat <- mat[rownames(diff_cluster_pygene_res),]
gene_col <- data.frame(Expression=c(rep("Up", 5), rep("Down",1)))
rownames(gene_col) <- rownames(mat)

newnames <- lapply(
  rownames(mat),
  function(x) bquote(italic(.(x))))

pheatplot <- pheatmap::pheatmap(mat=scale(mat),
                      annotation_col = annotation_col,
                      annotation_row = gene_col,
                      annotation_colors = ann_colors,
                      fontsize = 12,
                      labels_row = as.expression(newnames),
                      show_colnames = FALSE,
                      cluster_cols = F,
                      cluster_rows = F,
                      annotation_names_row = F,
                      gaps_row = 5,
                      gaps_col = 209)
ggsave(file="pheatmap.png", pheatplot)
ggsave(file="pheatmap.pdf", pheatplot)


# 12 enrichment -----------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./12_Enrichment")){
  dir.create("./12_Enrichment")
}
setwd("./12_Enrichment")

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(Ipaper)

gene_transform <- bitr(cluster_diff_gene,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID", "ENSEMBL", "REFSEQ"),
                       OrgDb = "org.Hs.eg.db")

# GO enrichment analysis
ego <- enrichGO(gene = gene_transform$ENTREZID, 
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", 
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)


write.table(ego@result,file = "GO.xls",sep = "\t",quote = F,row.names = F)
go_bar <- barplot(ego, showCategory=10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
  theme(axis.text.y = element_text(size  = 15))
go_bar
go_dot <- dotplot(ego, showCategory=10, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
  theme(axis.text.y = element_text(size  = 15))
go_dot

write_fig(go_bar,
          file = "GO_bar.pdf",
          width = 13,
          height = 12,
          devices = NULL,
          res = 600,
          show = F)
write_fig(go_bar,
          file = "GO_bar.png",
          width = 13,
          height = 12,
          devices = NULL,
          res = 600,
          show = F)
write_fig(go_dot,
          file = "GO_dot.pdf",
          width = 13,
          height = 12,
          devices = NULL,
          res = 600,
          show = F)
write_fig(go_dot,
          file = "GO_dot.png",
          width = 13,
          height = 12,
          devices = NULL,
          res = 600,
          show = F)
# KEGG enrichment analysis
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)

kk_bar <- barplot(kk, showCategory=20) +
  theme(axis.text.y = element_text(size  = 15))
kk_bar
kk_dot <- dotplot(kk, showCategory=20) +
  theme(axis.text.y = element_text(size  = 15))
kk_dot
write_fig(kk_bar,
          file = "KEGG_bar.pdf",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)
write_fig(kk_bar,
          file = "KEGG_bar.png",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)
write_fig(kk_dot,
          file = "KEGG_dot.pdf",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)
write_fig(kk_dot,
          file = "KEGG_dot.png",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)

# 13 GSEA -----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./13_GSEA")){
  dir.create("./13_GSEA")
}
setwd("./13_GSEA")

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

gene <- rownames(DEG2)
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
DEG2$SYMBOL <- rownames(DEG2)
data_all <- DEG2 %>%
  dplyr::inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>%
  arrange(desc(log2FoldChange))
geneList = data_all_sort$log2FoldChange
names(geneList) <- data_all_sort$ENTREZID
head(geneList)
write.table(geneList,
            file = "cluster_diff_all_gene_sort.txt",
            quote = F,
            col.names = F)

enrichplot::gseaplot2(gobp_gsea,
          geneSetID = 1:5, 
          # color="preranked", 
          base_size = 15, 
          # subplots = 1:2, 
          # pvalue_table = T
          ) 
ggsave(filename = "gobp_baseplot.pdf", height = 13, width = 15)
ggsave(filename = "gobp_baseplot.png", height = 13, width = 15)

enrichplot::gseaplot2(gocc_gsea,
                     geneSetID = 1:5,
                     base_size = 15) 
ggsave(filename = "gocc_baseplot.pdf", height = 13, width = 15)
ggsave(filename = "gocc_baseplot.png", height = 13, width = 15)

enrichplot::gseaplot2(gomf_gsea,
                      geneSetID = 1:5,
                      base_size = 15)
ggsave(filename = "gomf_baseplot.pdf", height = 13, width = 15)
ggsave(filename = "gomf_baseplot.png", height = 13, width = 15)

enrichplot::gseaplot2(kegg_gsea,
                      geneSetID = 1:5,
                      base_size = 15)
ggsave(filename = "kegg_baseplot.pdf", height = 13, width = 15)
ggsave(filename = "kegg_baseplot.png", height  =13, width = 15)

{
#GO-BP
gobp_set <- read.gmt("../data/c5.go.bp.v7.4.symbols.gmt")
gobp_gsea <- GSEA(geneList, TERM2GENE=gobp_set, verbose=FALSE, eps = 1e-100,
                  nPermSimple = 100000)
gobp_result <- gobp_gsea@result
NES = 1
gobp_result$change = as.factor(
  ifelse(gobp_result$p.adjust < 0.001 & abs(gobp_result$NES) > NES,
         ifelse(gobp_result$NES > NES,'UP','DOWN'),'NOT')
)
DEGeneSets <- subset(gobp_result,
                     gobp_result$p.adjust < 0.001 & abs(gobp_result$NES) > 1)
write.table(gobp_result,
            file = "gsea_gobp_result.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$p.adjust),]
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  data_repel <- rbind(DEGeneSets[1:10,], DEGeneSets[1:10,])
  data_repel$gobp <- rownames(data_repel)
  
  # Volcano
  library(ggplot2)
  library(ggthemes)
  library(Ipaper)
  library(ggrepel)
  gobp_result$change <- factor(gobp_result$change, 
                               levels = c("UP", "NOT", "DOWN"),
                               labels = c("Enriched in Cluster1", 
                                          "Not Significantly Enriched",
                                          "Enriched in Cluster2"))
  volcano_gobp <- ggplot(data = gobp_result,
                         aes(y = NES,
                             x = -log10(p.adjust), 
                             colour = change)) +
    scale_color_manual(values = c("#f8776e", "darkgray","#00c4c6"), name = "Cluster") +
    geom_point(aes(size = setSize/22880),alpha = 0.5, na.rm=T) +
    scale_size(range = c(0.3, 2), name = "GeneRatio") +
    xlim(0, 42) +
    geom_label_repel(
      data = up_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", 
      show.legend = FALSE,
      direction = "y",
      hjust = 0,
      force = 0.5,
      force_pull = 0,
      # nudge_x = 5,
      box.padding = 0.2,
      max.overlaps = 20,
      xlim = c(25,40)) +
    geom_label_repel(
      data = down_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", show.legend = FALSE,
      direction = "y",
      # nudge_x = 5,
      hjust = 0,
      force = 1,
      force_pull = 0,
      xlim = c(25,40),
      max.overlaps = 20) +
    theme_bw(base_size = 12, base_family = "Times") +
    geom_vline(xintercept = -log10(0.001),
               lty = 4,
               col = "darkgray",
               lwd = 0.6) +
    labs(y = "Normalized Enrichment Score",
         x = "-log10 (p-adjust)",
         title = "GSEA for Biological Process in GO Database") +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          legend.title = element_text(face = "bold", size = 17),
          legend.text = element_text(
            color="black",
            family = "Times",
            size=17),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
          axis.text.x = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.title.x = element_text(face = "bold",
                                      color = "black",
                                      size = 20),
          axis.title.y = element_text(face = "bold",
                                      color = "black",
                                      size = 20)) 
  
  volcano_gobp
  write_fig(volcano_gobp,
            file = "volcano_gobp.pdf",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(volcano_gobp,
            file = "volcano_gobp.png",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
}

#GO-CC
gocc_set <- read.gmt("../data/c5.go.cc.v7.4.symbols.gmt")
gocc_gsea <- GSEA(geneList, TERM2GENE = gocc_set, verbose = FALSE, eps = 1e-50, nPermSimple = 100000)
gocc_result <- gocc_gsea@result
gocc_result$change <- as.factor(ifelse(gocc_result$p.adjust < 0.001 & abs(gocc_result$NES) > NES,
                                       ifelse(gocc_result$NES > NES, "UP", "DOWN"), "NOT"))
DEGeneSets <- subset(gocc_result,
                     gocc_result$p.adjust < 0.001 & abs(gocc_result$NES) > 1)
write.table(gocc_result,
            file = "gsea_gocc_result.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$p.adjust),]
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  # Volcano
  library(ggplot2)
  library(ggthemes)
  library(Ipaper)
  library(ggrepel)
  gocc_result$change <- factor(gocc_result$change, levels = c("UP", "NOT", "DOWN"),
                               labels = c("Enriched in Cluster1", 
                                          "Not Significantly Enriched",
                                          "Enriched in Cluster2"))
  volcano_gocc <- ggplot(data = gocc_result,
                         aes(y = NES,
                             x = -log10(p.adjust), 
                             colour = change)) +
    scale_color_manual(values = c("#f8776e", "darkgray","#00c4c6"), name = "Cluster") +
    geom_point(aes(size = setSize/22880), alpha = 0.5, na.rm=T) +
    scale_size(range = c(0.3, 2), name = "GeneRatio") +
    xlim(0, 40) +
    geom_label_repel(
      data = up_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", 
      show.legend = FALSE,
      direction = "y",
      hjust = 0,
      force = 0.5,
      force_pull = 0,
      # nudge_x = 5,
      box.padding = 0.2,
      max.overlaps = 20,
      xlim = c(30,40)) +
    geom_label_repel(
      data = down_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", show.legend = FALSE,
      direction = "y",
      # nudge_x = 5,
      hjust = 0,
      force = 1,
      force_pull = 0,
      xlim = c(30,40),
      max.overlaps = 20) +
    theme_bw(base_size = 12, base_family = "Times") +
    geom_vline(xintercept = -log10(0.001),
               lty = 4,
               col = "darkgray",
               lwd = 0.6)+
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          legend.title = element_text(face = "bold", size = 17),
          legend.text = element_text(
            color="black",
            family = "Times",
            size=17),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
          axis.text.x = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.title.x = element_text(face = "bold",
                                      color = "black",
                                      size = 20),
          axis.title.y = element_text(face = "bold",
                                      color = "black",
                                      size = 20))  +
    labs(y = "Normalized Enrichment Score",
         x = "-log10 (p-adjust)",
         title = "GSEA for Cellular Component in GO Database")
  
  volcano_gocc
  write_fig(volcano_gocc,
            file = "volcano_gocc.pdf",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(volcano_gocc,
            file = "volcano_gocc.png",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
}

#GO-MF
gomf_set <- read.gmt("../data/c5.go.mf.v7.4.symbols.gmt")
gomf_gsea <- GSEA(geneList, TERM2GENE = gomf_set, verbose = FALSE, eps = 1e-50, nPermSimple = 100000)
gomf_result <- gomf_gsea@result
gomf_result$change <- ifelse(gomf_result$p.adjust < 0.001 & abs(gomf_result$NES) > NES,
                             ifelse(gomf_result$NES > NES, "UP", "DOWN"), "NOT") %>% as.factor()
DEGeneSets <- subset(gomf_result,
                     gomf_result$p.adjust < 0.001 & abs(gomf_result$NES) > 1)
write.table(gomf_result,
            file = "gsea_gomf_result.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$p.adjust),]
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  
  # Volcano
  library(ggplot2)
  library(ggthemes)
  library(Ipaper)
  library(ggrepel)
  gomf_result$change <- factor(gomf_result$change, levels = c("UP", "NOT", "DOWN"),
                               labels = c("Enriched in Cluster1", 
                                          "Not Significantly Enriched",
                                          "Enriched in Cluster2"))
  volcano_gomf <- ggplot(data = gomf_result,
                         aes(y = NES,
                             x = -log10(p.adjust), 
                             colour = change)) +
    scale_color_manual(values = c("#f8776e", "darkgray","#00c4c6"), name = "Cluster") +
    geom_point(aes(size = setSize/22880), alpha = 0.5, na.rm=T) +
    scale_size(range = c(0.3, 2), name = "GeneRatio") +
    xlim(0, 31) +
    geom_label_repel(
      data = up_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", 
      show.legend = FALSE,
      direction = "y",
      hjust = 0,
      force = 0.5,
      force_pull = 0,
      # nudge_x = 5,
      box.padding = 0.2,
      max.overlaps = 20,
      xlim = c(20,30)) +
    geom_label_repel(
      data = down_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", show.legend = FALSE,
      direction = "y",
      # nudge_x = 5,
      hjust = 0,
      force = 1,
      force_pull = 0,
      xlim = c(20,30),
      max.overlaps = 20) +
    theme_bw(base_size = 12, base_family = "Times") +
    geom_vline(xintercept = -log10(0.001),
               lty = 4,
               col = "darkgray",
               lwd = 0.6)+
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          legend.title = element_text(face = "bold", size = 17),
          legend.text = element_text(
            color="black",
            family = "Times",
            size=17),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
          axis.text.x = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.title.x = element_text(face = "bold",
                                      color = "black",
                                      size = 20),
          axis.title.y = element_text(face = "bold",
                                      color = "black",
                                      size = 20)) +
    labs(y = "Normalized Enrichment Score",
         x = "-log10 (p-adjust)",
         title = "GSEA for Molecular Function in GO Database")
  
  volcano_gomf
  write_fig(volcano_gomf,
            file = "volcano_gomf.pdf",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(volcano_gomf,
            file = "volcano_gomf.png",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
}

# KEGG
kegg_set <- read.gmt("../data/c2.cp.kegg.v7.4.symbols.gmt")
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, verbose = FALSE, eps = 1e-50)
kegg_result <- kegg_gsea@result

kegg_result$change <- ifelse(kegg_result$p.adjust < 0.001 & abs(kegg_result$NES) > NES,
                             ifelse(kegg_result$NES > NES, "UP", "DOWN"), "NOT") %>% as.factor()
DEGeneSets <- subset(kegg_result,
                     kegg_result$p.adjust < 0.001 & abs(kegg_result$NES) > 1)
write.table(kegg_result,
            file = "gsea_kegg_result.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$p.adjust),]
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  
  # Volcano
  library(ggplot2)
  library(ggthemes)
  library(Ipaper)
  library(ggrepel)
  kegg_result$change <- factor(kegg_result$change, levels = c("UP", "NOT", "DOWN"),
                               labels = c("Enriched in Cluster1", 
                                          "Not Significantly Enriched",
                                          "Enriched in Cluster2"))
  volcano_kegg <- ggplot(data = kegg_result,
                         aes(y = NES,
                             x = -log10(p.adjust), 
                             colour = change)) +
    scale_color_manual(values = c("#f8776e", "darkgray","#00c4c6"), name = "Cluster") +
    geom_point(aes(size = setSize/22880), alpha = 0.5, na.rm=T) +
    scale_size(range = c(0.3, 2), name = "GeneRatio") +
    xlim(0, 20) +
    geom_label_repel(
      data = up_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", 
      show.legend = FALSE,
      direction = "y",
      hjust = 0,
      force = 0.5,
      force_pull = 0,
      # nudge_x = 5,
      box.padding = 0.2,
      max.overlaps = 20,
      xlim = c(13,18)) +
    geom_label_repel(
      data = down_DEG[1:10,],
      aes(label = gobp),
      fontface = "italic",
      size = 3.5,
      color = "black",
      segment.color = "black", show.legend = FALSE,
      direction = "y",
      # nudge_x = 5,
      hjust = 0,
      force = 1,
      force_pull = 0,
      xlim = c(13,18),
      max.overlaps = 20) +
    theme_bw(base_size = 12, base_family = "Times") +
    geom_vline(xintercept = -log10(0.001),
               lty = 4,
               col = "darkgray",
               lwd = 0.6)+
    guides(color = guide_legend(order = 1)) +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          legend.title = element_text(face = "bold", size = 17),
          legend.text = element_text(
            color="black",
            family = "Times",
            size=17),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
          axis.text.x = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     color = "black",
                                     size = 15),
          axis.title.x = element_text(face = "bold",
                                      color = "black",
                                      size = 20),
          axis.title.y = element_text(face = "bold",
                                      color = "black",
                                      size = 20))  +
    labs(y = "Normalized Enrichment Score",
         x = "-log10 (p-adjust)",
         title = "GSEA for Pathways in KEGG Database") 
  volcano_kegg
  write_fig(volcano_kegg,
            file = "volcano_kegg.pdf",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(volcano_kegg,
            file = "volcano_kegg.png",
            width = 20,
            height = 9,
            devices = NULL,
            res = 600,
            show = F)
}
}


# 13.1 GSVA ---------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./GSVA")){
  dir.create("./GSVA")
}
setwd("./GSVA")
library(GSVA)
library(GSEABase)
library(limma)

gsva_exp <- fpkm_tumor[, cluster_list$sample]
# 分组
group_list2 <- cluster_list$cluster %>% as.factor()
design_risk <- model.matrix(~0 + group_list2)
rownames(design_risk) <- colnames(gsva_exp)
colnames(design_risk) <- levels(group_list2)

compare_risk <- makeContrasts("Cluster1-Cluster2", levels = design_risk)


# GO-BP
GOBP_ref <- getGmt("data/c5.go.bp.v7.4.symbols.gmt")
es_GOBP <- gsva(as.matrix(gsva_exp), GOBP_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOBP <- as.data.frame(es_GOBP)

fit <- lmFit(es_GOBP, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GOBP.xls",
            quote = F,
            sep = "\t",
            row.names = T)

{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  
  up_DEG <- up_DEG[order(-up_DEG$logFC, up_DEG$adj.P.Val),]
  down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$adj.P.Val),]
  
  data_repel <- rbind(up_DEG[1:20,], down_DEG[1:20,])
  
  mat <- es_GOBP[rownames(data_repel),]
  mat <- mat[, cluster_list$sample]
  annotation_col <- data.frame(
    Cluster = factor(cluster_list$cluster),
    row.names = cluster_list$sample
  )
  ann_colors <- list(
    Cluster = c(Cluster1="#f8776e", Cluster2="#00c4c6"),
    GSVA_GO_BP = c( Up="red",Down="blue")
  )
  
  head(data_repel[order(data_repel$change, decreasing = T),]$change)
  
  mat <- mat[rownames(data_repel[order(data_repel$change, decreasing = T),]),]
  gene_col <- data.frame(GSVA_GO_BP=c(rep("Up",20), rep("Down", 20)))
  rownames(gene_col) <- rownames(mat)
  library(pheatmap)
  pheatplot <- pheatmap(mat=mat,
                        annotation_col = annotation_col,
                        # annotation_row = gene_col,
                        annotation_colors = ann_colors,
                        fontsize = 12,
                        show_colnames = FALSE,
                        cluster_cols = F,
                        cluster_rows = T,
                        annotation_names_row = F)
  ggsave(file="gsva_gobp_pheatmap.png", pheatplot, width = 20, height = 10)
  ggsave(file="gsva_gobp_pheatmap.pdf", pheatplot, width = 20, height = 10)
 
}

# GO-CC
GOCC_ref <- getGmt("data/c5.go.cc.v7.4.symbols.gmt")
es_GOCC <- gsva(as.matrix(gsva_exp), GOCC_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOCC <- as.data.frame(es_GOCC)
fit <- lmFit(es_GOCC, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GOCC.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
  dim(DEGeneSets)
  # [1] 369   7
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  
  up_DEG <- up_DEG[order(-up_DEG$logFC, up_DEG$adj.P.Val),]
  down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$adj.P.Val),]
  
  data_repel <- rbind(up_DEG[1:20,], down_DEG[1:20,])
  
  mat <- es_GOCC[rownames(data_repel),]
  mat <- mat[, cluster_list$sample]
  annotation_col <- data.frame(
    Cluster = factor(cluster_list$cluster),
    row.names = cluster_list$sample
  )
  ann_colors <- list(
    Cluster = c(Cluster1="#f8776e", Cluster2="#00c4c6"),
    GSVA_GO_CC = c( Up="red",Down="blue")
  )
  
  head(data_repel[order(data_repel$change, decreasing = T),]$change)
  
  mat <- mat[rownames(data_repel[order(data_repel$change, decreasing = T),]),]
  gene_col <- data.frame(GSVA_GO_CC=c(rep("Up",20), rep("Down", 20)))
  rownames(gene_col) <- rownames(mat)
  library(pheatmap)
  pheatplot <- pheatmap(mat=mat,
                        annotation_col = annotation_col,
                        # annotation_row = gene_col,
                        annotation_colors = ann_colors,
                        fontsize = 12,
                        show_colnames = FALSE,
                        cluster_cols = F,
                        cluster_rows = T,
                        annotation_names_row = F)
  ggsave(file="gsva_gocc_pheatmap.png", pheatplot, width = 20, height = 10)
  ggsave(file="gsva_gocc_pheatmap.pdf", pheatplot, width = 20, height = 10)
}

# GO-MF
GOMF_ref <- getGmt("data/c5.go.mf.v7.4.symbols.gmt")
es_GOMF <- gsva(as.matrix(gsva_exp), GOMF_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOMF <- as.data.frame(es_GOMF)
fit <- lmFit(es_GOMF, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GOMF.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  
  up_DEG <- up_DEG[order(-up_DEG$logFC, up_DEG$adj.P.Val),]
  down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$adj.P.Val),]
  
  data_repel <- rbind(up_DEG[1:20,], down_DEG[1:20,])
  
  mat <- es_GOMF[rownames(data_repel),]
  mat <- mat[, cluster_list$sample]
  annotation_col <- data.frame(
    Cluster = factor(cluster_list$cluster),
    row.names = cluster_list$sample
  )
  ann_colors <- list(
    Cluster = c(Cluster1="#f8776e", Cluster2="#00c4c6"),
    GSVA_GO_MF = c( Up="red",Down="blue")
  )
  
  head(data_repel[order(data_repel$change, decreasing = T),]$change)
  
  mat <- mat[rownames(data_repel[order(data_repel$change, decreasing = T),]),]
  gene_col <- data.frame(GSVA_GO_MF=c(rep("Up",20), rep("Down", 20)))
  rownames(gene_col) <- rownames(mat)
  library(pheatmap)
  pheatplot <- pheatmap(mat=mat,
                        annotation_col = annotation_col,
                        # annotation_row = gene_col,
                        annotation_colors = ann_colors,
                        fontsize = 15,
                        show_colnames = FALSE,
                        cluster_cols = F,
                        cluster_rows = T,
                        annotation_names_row = F)
  ggsave(file="gsva_gomf_pheatmap.png", pheatplot, width = 40, height = 18)
  ggsave(file="gsva_gomf_pheatmap.pdf", pheatplot, width = 40, height = 18)
}

# KEGG
KEGG_ref <- getGmt("data/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), KEGG_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
{
  DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
  dim(DEGeneSets)
  # [1] 118   7
  up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
  down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
  
  up_DEG$gobp <- rownames(up_DEG)
  down_DEG$gobp <- rownames(down_DEG)
  
  up_DEG <- up_DEG[order(-up_DEG$logFC, up_DEG$adj.P.Val),]
  down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$adj.P.Val),]
  
  data_repel <- rbind(up_DEG[1:20,], down_DEG[1:20,])
  
  mat <- es_KEGG[rownames(data_repel),]
  mat <- mat[, cluster_list$sample]
  annotation_col <- data.frame(
    Cluster = factor(cluster_list$cluster),
    row.names = cluster_list$sample
  )
  ann_colors <- list(
    Cluster = c(Cluster1="#f8776e", Cluster2="#00c4c6"),
    GSVA_KEGG = c( Up="red",Down="blue")
  )
  
  head(data_repel[order(data_repel$change, decreasing = T),]$change)
  
  mat <- mat[rownames(data_repel[order(data_repel$change, decreasing = T),]),]
  gene_col <- data.frame(GSVA_KEGG=c(rep("Up",20), rep("Down", 20)))
  rownames(gene_col) <- rownames(mat)
  library(pheatmap)
  pheatplot <- pheatmap(mat=mat,
                        annotation_col = annotation_col,
                        # annotation_row = gene_col,
                        annotation_colors = ann_colors,
                        fontsize = 12,
                        show_colnames = FALSE,
                        cluster_cols = F,
                        cluster_rows = T,
                        annotation_names_row = F)
  ggsave(file="gsva_kegg_pheatmap.png", pheatplot, width = 20, height = 10)
  ggsave(file="gsva_kegg_pheatmap.pdf", pheatplot, width = 20, height = 10)

}


# 14 Univariate cox analysis---------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./14_Unix_cox")){
  dir.create("./14_Unix_cox")
}
setwd("./14_Unix_cox")


diff_expr <- fpkm_tumor[cluster_diff_gene,]
write.table(diff_expr,
            file = "diff_expr.xls",
            quote = F,
            sep = "\t")
diff_expr_t <- as.data.frame(t(log2(diff_expr+1)))
diff_expr_t$sample <- rownames(diff_expr_t)
diff_expr_clinical <- merge(diff_expr_t,
                            clinical_data2,
                            by = "sample")
rownames(diff_expr_clinical) <- diff_expr_clinical$sample
diff_expr_clinical <- diff_expr_clinical[,c(rownames(diff_expr), "OS", "OS.time")]
diff_expr_clinical <- diff_expr_clinical[-which(diff_expr_clinical$OS.time == 0),]
dim(diff_expr_clinical)
# [1] 405 956
write.table(diff_expr_clinical,
            file = "diff_expr_clinical.xls",
            row.names = T,
            quote = F,
            sep = "\t")

library(survival)
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))

univ_models <- lapply(univ_formulas,
                      function(x) coxph(x, data = diff_expr_clinical))
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #p value
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #HR
                         HR <-signif(x$coef[2], digits=3);
                         #95% CI
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
write.table(res,
            file = "univariate_cox_result.csv",
            quote = F,
            row.names = T)
res_results_0.05 <- res[which(as.numeric(res$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.csv",
            quote = F,
            row.names = T)
dim(res_results_0.05)

res_results_0.01 <- res[which(as.numeric(res$p.value) < 0.01),]
res_results_0.01 <- na.omit(res_results_0.01)
write.table(res_results_0.01,
            file = "univariate_cox_result_0.01.csv",
            quote = F,
            row.names = T)
dim(res_results_0.01)

res_results_0.001 <- res[which(as.numeric(res$p.value) < 0.001),]
res_results_0.001 <- na.omit(res_results_0.001)
write.table(res_results_0.001,
            file = "univariate_cox_result_0.001.csv",
            quote = F,
            row.names = T)
dim(res_results_0.001)

res_results_0.0001 <- res[which(as.numeric(res$p.value) < 0.0001),]
res_results_0.0001 <- na.omit(res_results_0.0001)
write.table(res_results_0.0001,
            file = "univariate_cox_result_0.0001.csv",
            quote = F,
            row.names = T)
dim(res_results_0.0001)

library(tidyr)
res_results_0.05_2 <- separate(res_results_0.0001, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]
hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.0001,
                                         "< 0.0001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = "univariate_cox_forest.pdf", height = 16, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95% CI lower
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95% CI upper
           boxsize=0.2,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,    
           colgap=unit(5,"mm"),    
           xticks = c(0, 1, 2, 3), 
           lwd.xaxis=2,            
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.5,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1.3, fontface = "bold"),
                          title=gpar(cex = 1.5, fontface = "bold")),
           xlab="Hazard Ratio",
           cex.lab=2,
           grid = T) 

dev.off()
png(filename = "univariate_cox_forest.png", height = 1100, width = 800)
forestplot(labeltext=tabletext, 
           graph.pos=4,  
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), 
           upper=c(NA,NA,res_results_0.05_2$HR.95H), 
           boxsize=0.2,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,   
           colgap=unit(5,"mm"),   
           xticks = c(0, 1, 2, 3), 
           lwd.xaxis=2,            
           lineheight = unit(0.8,"cm"),
           graphwidth = unit(.5,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"),
           #                 "59" = gpar(lwd=2, col="black")),
           # mar=unit(rep(0.01, times = 4), "cm"),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1.3, fontface = "bold"),
                          title=gpar(cex = 1.5, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) 

dev.off()



# 15 lasso ----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./15_Lasso")){
  dir.create("./15_Lasso")
}
setwd("./15_Lasso")

library(glmnet)
x_all <- subset(diff_expr_clinical, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.0001)]
y_all <- subset(diff_expr_clinical, select = c(OS, OS.time))

# Model fitting
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
plot(fit, xvar = "lambda",label = TRUE, las=1)

png(filename = "lasso_model.png", height = 450, width = 600)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()


# Cross validation
set.seed(33)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),
                  nfold=50,
                  family = "cox") 
plot(cvfit, las =1)

png(filename = "lasso_verify.png", height = 450, width = 600)
plot(cvfit, las =1)
dev.off()
pdf(file = "lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()

# get coefficients
coef.min = coef(cvfit, s = "lambda.min")
cvfit$lambda.min
active.min = which(coef.min@i != 0)

# get geneids
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
write(lasso_geneids, "lasso_genes.csv")
write.csv(x_all,file = "Lasso_x.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)


# 16 multivariate analysis ------------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./16_Mul_cox")){
  dir.create("./16_Mul_cox")
}
setwd("./16_Mul_cox")

library(survival)
library(survminer)
library(Ipaper)
gene_list <- lasso_geneids
cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                              paste(gene_list,
                                    sep = '',
                                    collapse = '+')))
cox_more <- coxph(cox_data,
                  data = diff_expr_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]

cox_formula <- as.formula(paste("Surv(OS.time, OS)~",
                                paste(rownames(cox_table)[cox_table[,3]>0.05],
                                      collapse = "+")))
cox_more_2 <- coxph(cox_formula, data = diff_expr_clinical)
mul_cox_result <- summary(cox_more_2)$coefficients
mul_cox_gene <- rownames(mul_cox_result)
mul_cox_gene
library(survminer)
train_HRforest <- ggforest(model = cox_more_2,
                           data = diff_expr_clinical,
                           main = "Hazard ratio of Candidate Genes",
                           fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
train_HRforest
write_fig(train_HRforest,
          file = "train_HRforest.pdf",
          width = 10,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(train_HRforest,
          file = "train_HRforest.png",
          width = 10,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)

riskScore=data.frame(riskScore=predict(cox_more_2,type="lp",newdata=diff_expr_clinical))
coxGene=rownames(as.data.frame(cox_more_2$coefficients))
outCol=c("OS","OS.time",coxGene)
risk <- subset(diff_expr_clinical, select = outCol)
risk <- merge(riskScore, risk, by = "row.names") %>% tibble::column_to_rownames(var="Row.names")
risk$risk <- ifelse(risk$riskScore>median(risk$riskScore), "High", "Low")
{
  library(ggplot2)
  library(ggthemes)
  library(Ipaper)
  risk_dis <- ggplot(risk, aes(x=reorder(id, riskScore), 
                               y=riskScore, 
                               color = factor(risk, 
                                              levels = c(0, 1), 
                                              labels = c("High Risk", "Low Risk")))) +
    geom_point() +
    scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
    scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400)],
                     labels = c(1,100,200,300,400),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
               lty = 2) +
    geom_hline(yintercept = median(riskScore),
               lty =2) +
    labs(x = "Patients(increasing risk score)",
         y = "Risk Score",
         title = "Train Risk Score Distribution") + 
    theme_base() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          plot.title = element_text(size = 15, hjust = 0.5))
  risk_dis
  write_fig(risk_dis,
            file = "train_riskscore_dis.pdf",
            width = 10.5,
            height = 5.4,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(risk_dis,
            file = "train_riskscore_dis.png",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  
  surv_stat <- ggplot(risk, aes(x=reorder(id, riskScore),
                                y=OS.time,
                                color = factor(OS,
                                               levels = c(0,1),
                                               labels = c("Alive", "Dead")))) +
    geom_point() +
    scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
    scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400)],
                     labels = c(1,100,200,300,400),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
               lty = 2) +
    labs(x = "Patients(increasing risk score)",
         y = "Survival time (days)",
         title = "Train Survival State Distribution") + 
    theme_base() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          plot.title = element_text(size = 15, hjust = 0.5))
  surv_stat
  write_fig(surv_stat,
            file = "train_survstat_dis.pdf",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(surv_stat,
            file = "train_survstat_dis.png",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}


kmfit<-survfit(Surv(OS.time, OS) ~ risk, data =  risk)
train_survival_median <- ggsurvplot(kmfit,
                                    pval = TRUE, 
                                    conf.int = F,
                                    legend.labs=c("High risk","Low risk" ),
                                    legend.title="Risk score",
                                    title="Train KM",
                                    font.main = c(15,"bold"),
                                    risk.table = TRUE, 
                                    risk.table.col = "strata", 
                                    linetype = "strata", 
                                    # surv.median.line = "hv", 
                                    ggtheme = theme_bw(), 
                                    palette = c("#A73030FF", "#0073C2FF"))
train_survival_median

train_survival_median$table <- train_survival_median$table + 
  labs(x = "Overall Survival (days)") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 15, face = "bold"))
train_survival_median$plot <- train_survival_median$plot + 
  labs(x = "Overall Survival (days)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 13, face = "bold"),
        text = element_text(size = 20, face = "bold"))
train_survival_median

write_fig(train_survival_median,
          file = "train_survival_median.pdf",
          width = 5,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(train_survival_median,
          file = "train_survival_median.png",
          width = 5,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)



# ROC
riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}
candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(diff_expr_clinical,
                                         candidate_genes_for_cox2,
                                         cox_more_2)

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# ROC curves 
library(scales)
library(geomROC)
library(plotROC)

auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC1 <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#BC392F", "#076DB2", "#DD822D")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.position = c(1,0),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5),) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))),
           size = 5)
ROC1
write_fig(ROC1,
          file = "train_ROC.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(ROC1,
          file = "train_ROC.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)

{
  library(survivalROC)
  library(tidyverse)
  rt = subset(risk, select = c(OS, OS.time, riskScore))
  rt$OS.time <- rt$OS.time / 365
  survivalROC_helper <- function(t) {
    
    survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =t, method="KM")
  }
  survivalROC_data <- data_frame(t = c(1,3,5)) %>%
    
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
}
survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.2f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)
survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2))

survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.2f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

year =factor(survivalROC_data1$year)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  scale_color_manual(
                     values = c("#BC392F", "#076DB2", "#DD822D")) +
  geom_path(aes(color= year))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  facet_wrap( ~ year) +
  theme_bw() +
  labs(x = "False positive fraction", y = "True positive fraction", title = "Train ROC") +
  theme(axis.text.x = element_text(vjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        legend.position = "none",
        panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15))

ggsave(filename = "train_ROC2.pdf", height = 6, width = 15)
ggsave(filename = "train_ROC2.png", height = 6, width = 15)

# 17 Model validation -----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./17_risk_model_verify")){
  dir.create("./17_risk_model_verify")
}
setwd("./17_risk_model_verify")

# library(GEOquery)
# eSet_out <- getGEO('GSE13507', destdir='./', getGPL=F)
# expr_out <- exprs(eSet_out[[1]])
# dim(expr_out)
# # [1] 43148   256
# phenotype_out <- pData(eSet_out[[1]])
# clinical_out <- data.frame(OS.time=as.numeric(phenotype_out$`survival month:ch1`)*30,
#                            OS=phenotype_out$`overall survival:ch1`)
# 
# View(phenotype_out)

expr_out2 <- read.table("../data/GSE13507_exp.txt",header = T)
phenotype_out2 <- read.table("../data/pheno.txt", header = T, sep="\t")
clinical_out2 <- data.frame(sample=phenotype_out2$Sample.name,
                            OS.time=round(phenotype_out2$survivalMonth*31),
                            OS=ifelse(phenotype_out2$overall.survival==1,0,1))
clinical_out_sample <- clinical_out2$sample
expr_out2 <- expr_out2[,clinical_out_sample]
gse_all_gene <- rownames(expr_out2)
expr_out3 <- expr_out2[mul_cox_gene, clinical_out_sample]
verify_data <- as.data.frame(t((log2(expr_out3+1))))
verify_data$sample <- rownames(verify_data)
verify_expr_clinical <- merge(verify_data,
                              clinical_out2,
                              by = "sample")
rownames(verify_expr_clinical) <- verify_expr_clinical$sample
colnames_sum <- colnames(verify_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(verify_expr_clinical) <- colnames_sum

riskScore_out = predict(cox_more_2,type="lp",newdata=verify_expr_clinical)
risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
risk_out <- as.data.frame(c(cbind(id=rownames(cbind(verify_expr_clinical[,outCol],
                                                    riskScore_out,
                                                    risk_out)),
                                  cbind(verify_expr_clinical[,outCol],
                                        riskScore_out,
                                        risk_out))))
{
  library(ggplot2)
  library(ggthemes)
  median(riskScore_out)
  # [1] 1.678871
  risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                       y=riskScore_out, 
                                       color = factor(risk_out, 
                                                      levels = c(0, 1), 
                                                      labels = c("High Risk", "Low Risk")))) +
    geom_point() +
    scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
    scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,20,40,60,80, 100,120,140,160)],
                     labels = c(1,20,40,60,80, 100,120,140,160),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
               lty = 2) +
    geom_hline(yintercept = median(riskScore_out),
               lty =2) +
    labs(x = "Patients(increasing risk score)",
         y = "Risk Score",
         title = "Validation Risk Score Distribution") + 
    theme_base() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          plot.title = element_text(size = 15, hjust = 0.5))
  risk_dis_out
  write_fig(risk_dis_out,
            file = "verify_riskscore_dis.pdf",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(risk_dis_out,
            file = "verify_riskscore_dis.png",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  
  surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                        y=OS.time,
                                        color = factor(OS,
                                                       levels = c(0,1),
                                                       labels = c("Alive", "Dead")))) +
    geom_point() +
    scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
    scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,20,40,60,80, 100,120,140,160)],
                     labels = c(1,20,40,60,80, 100,120,140,160),
                     expand = c(0.02,0)) +
    ylim(x=c(0,5000)) +
    geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
               lty = 2) +
    labs(x = "Patients(increasing risk score)",
         y = "Survival time (days)",
         title = "Validation Survival State Distribution") + 
    theme_base() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          plot.title = element_text(size = 15, hjust = 0.5))
  surv_stat_out
  write_fig(surv_stat_out,
            file = "verify_survstat_dis.pdf",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
  write_fig(surv_stat_out,
            file = "verify_survstat_dis.png",
            width = 8,
            height = 5,
            devices = NULL,
            res = 600,
            show = F)
}


library(survival)
library(survminer)
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
verify_survival_median <- ggsurvplot(kmfit_out,
                                     pval = TRUE, 
                                     conf.int = F,
                                     legend.labs=c("High risk","Low risk" ),
                                     legend.title="Risk score",
                                     title="Validition KM",
                                     font.main = c(15,"bold"),
                                     risk.table = TRUE, 
                                     risk.table.col = "strata", 
                                     linetype = "strata", 
                                     # surv.median.line = "hv", 
                                     ggtheme = theme_bw(), 
                                     palette = c("#A73030FF", "#0073C2FF"))
verify_survival_median
verify_survival_median$table <- verify_survival_median$table + 
  labs(x = "Overall Survival (days)") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 15, face = "bold"))
verify_survival_median$plot <- verify_survival_median$plot + 
  labs(x = "Overall Survival (days)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 13, face = "bold"),
        text = element_text(size = 20, face = "bold"))
verify_survival_median
write_fig(verify_survival_median,
          file = "verify_survival_median.pdf",
          width = 5,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(verify_survival_median,
          file = "verify_survival_median.png",
          width = 5,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)

riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}
candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(verify_expr_clinical,
                                         candidate_genes_for_cox2,
                                         cox_more_2)

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_out,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk_out)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# ROC curves
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC3 <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#BC392F", "#076DB2", "#DD822D")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Validation ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.position = c(1,0),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5),) +
          annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))),
           size = 5)

ROC3

write_fig(ROC3,
          file = "verify_ROC.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(ROC3,
          file = "verify_ROC.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
{
  library(survivalROC)
  library(tidyverse)
  rt = subset(risk_out, select = c(OS, OS.time, riskScore_out))
  rt$OS.time <- rt$OS.time / 365
  survivalROC_helper <- function(t) {
    
    survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore_out, 
                predict.time =t, method="KM")
    
  }
  survivalROC_data <- data_frame(t = c(1,3,5)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
}
survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.2f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2))

survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.2f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")


year =factor(survivalROC_data1$year)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  scale_color_manual(
    values = c("#BC392F", "#076DB2", "#DD822D")) +
  geom_path(aes(color= year))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  facet_wrap( ~ year) +
  theme_bw() +
  labs(x = "False positive fraction", y = "True positive fraction", title = "Validation ROC") +
  theme(axis.text.x = element_text(vjust = 0.5),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        legend.position = "none",
        panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15))

ggsave(filename = "verify_ROC2.pdf", height = 6, width = 15)
ggsave(filename = "verify_ROC2.png", height = 6, width = 15)

# 18 Nomogram model -----------------------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./18_prog_model")){
  dir.create("./18_prog_model")
}
setwd("./18_prog_model")


clinical_data3 <- clinical_data2
clinical_data3$T_stage <- gsub("T0", NA, clinical_data3$T_stage)
clinical_data3$T_stage <- gsub("TX", NA, clinical_data3$T_stage)
clinical_data3$T_stage <- gsub("^T1[a-c]*$", "T1/T2", clinical_data3$T_stage)
clinical_data3$T_stage <- gsub("^T2[a-c]*$", "T1/T2", clinical_data3$T_stage)
clinical_data3$T_stage <- gsub("^T3[a-c]*$", "T3/T4", clinical_data3$T_stage)
clinical_data3$T_stage <- gsub("^T4[a-c]*$", "T3/T4", clinical_data3$T_stage)

clinical_data3$N_stage <- gsub("NX", NA, clinical_data3$N_stage)
clinical_data3$M_stage <- gsub("MX", NA, clinical_data3$M_stage)
sub_risk <- subset(risk, select = c(id, riskScore))
names(sub_risk) <- c("sample", "riskScore")
train_risk_clinical <- merge(clinical_data3,
                             sub_risk,
                             by = "sample")
rownames(train_risk_clinical) <- train_risk_clinical$sample
train_risk_clinical = subset(train_risk_clinical, select = -c(sample))
dim(train_risk_clinical)

train_risk_clinical$gender <- factor(train_risk_clinical$gender)
train_risk_clinical$stage <- factor(train_risk_clinical$stage)
train_risk_clinical$T_stage <- factor(train_risk_clinical$T_stage)
train_risk_clinical$N_stage <- factor(train_risk_clinical$N_stage)
train_risk_clinical$M_stage <- factor(train_risk_clinical$M_stage)

library(survival)
library(survminer)
library(Ipaper)

colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]
covariates_train

covariates_train <- c(
                      # "gender",
                      "age",
                      # "stage",
                      # "T_stage",
                      "N_stage",
                      # "M_stage",
                      "riskScore")
cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                              paste(covariates_train,
                                    sep = '',
                                    collapse = '+')))
cox_more <- coxph(cox_data,
                  data = train_risk_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]

cox_formula <- as.formula(paste("Surv(OS.time, OS)~",
                                paste(rownames(cox_table)[cox_table[,3]>0.05],
                                      collapse = "+")))
cox_more_2 <- coxph(cox_formula, data = train_risk_clinical)

library(survminer)
train_HRforest_prog <- ggforest(model = cox_more_2,
                           data = train_risk_clinical,
                           main = "",
                           fontsize = 1.3) +
  theme(plot.title = element_text(face = "bold", size = 20))
train_HRforest_prog

write_fig(train_HRforest_prog,
          file = "train_HRforest_prog2.pdf",
          width = 13,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(train_HRforest_prog,
          file = "train_HRforest_prog2.png",
          width = 13,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)

write_fig(train_HRforest_prog,
          file = "train_HRforest_prog.pdf",
          width = 13,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
write_fig(train_HRforest_prog,
          file = "train_HRforest_prog.png",
          width = 13,
          height = 5,
          devices = NULL,
          res = 600,
          show = F)
# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# Cox model
res.cox <- psm(cox_data,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) 
function(x) surv(365, x) # 1-year
function(x) surv(1095, x) # 3-year
function(x) surv(1825, x) # 5-year

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99))
plot(nom.cox, cex.axis = 1.5)
png(filename = "nomogram_line_points.png", height = 650, width = 900)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "nomogram_line_points.pdf", height = 9, width = 11)
plot(nom.cox, cex.axis = 1.5, cex.var = 1.7)
dev.off()

coxm_1 <- cph(cox_data,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=100,B=1000)
coxm_3 <- cph(cox_data,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=100,B=1000)
coxm_5 <- cph(cox_data,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=100,B=1000)


png(filename = "nomogram_predicted_1-5_2021-11-22.png", height = 550, width = 800)
par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',
     ylab='Actual 1-5 year Overall Survival(proportion)',
     col="#BC392F",
     xlim = c(0,1),ylim = c(0,1)) 
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',
     ylab='Actual 1-5 year Overall Survival(proportion)',
     col="#076DB2",
     xlim = c(0,1),ylim = c(0,1))
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlab='Nomogram-Predicted Probability of 1-year Overall Survival',
     ylab='Actual 1-year Overall Survival(proportion)',
     col="#DD822D",
     xlim = c(0,1),ylim = c(0,1))


# Legend
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#BC392F", "#076DB2", "#DD822D"), 
       lwd=2)
# Diagonal
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

pdf(file = "nomogram_predicted_1-5_2021-11-22.pdf", height = 7, width = 9)
par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',
     ylab='Actual 1-5 year Overall Survival(proportion)',
     col="#BC392F",
     xlim = c(0,1),ylim = c(0,1)) 
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',
     ylab='Actual 1-5 year Overall Survival(proportion)',
     col="#076DB2",
     xlim = c(0,1),ylim = c(0,1)) 
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlab='Nomogram-Predicted Probability of 1-year Overall Survival',
     ylab='Actual 1-year Overall Survival(proportion)',
     col="#DD822D",
     xlim = c(0,1),ylim = c(0,1)) 


#Legend
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#BC392F", "#076DB2", "#DD822D"), 
       lwd=2)
#Diagonal
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()


# 19 Risk model and Pyroptosis-related genes --------------------------------------------------

setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./19_gene_corr")){
  dir.create("./19_gene_corr")
}
setwd("./19_gene_corr")

sig_km_gene_expr <- fpkm_tumor[which(rownames(fpkm_tumor)%in%sig_km_gene),]
dim(sig_km_gene_expr)

mul_cox_gene_expr <- fpkm_tumor[which(rownames(fpkm_tumor)%in%mul_cox_gene),]
dim(mul_cox_gene_expr)

CHMP4C_expr <- fpkm_tumor["CHMP4C",]
GSDMA_expr <- fpkm_tumor[which(rownames(fpkm_tumor)=="GSDMA"),]

CHMP4C_cor_df <- fpkm_tumor[which(rownames(fpkm_tumor)%in%c(mul_cox_gene,"CHMP4C")),]
GSDMA_cor_df <- fpkm_tumor[which(rownames(fpkm_tumor)%in%c(mul_cox_gene,"GSDMA")),]
CHMP4C_cor_df <- log2(CHMP4C_cor_df+1)
GSDMA_cor_df <- log2(GSDMA_cor_df+1)

y <- as.numeric(CHMP4C_cor_df["CHMP4C",])
rownames <- rownames(CHMP4C_cor_df)
CHMP4C_cor_res <- data.frame(rownames)
for (i in 1:length(rownames)){
  test <- cor.test(as.numeric(CHMP4C_cor_df[i,]),y,type="spearman")
  CHMP4C_cor_res[i,2] <- test$estimate
  CHMP4C_cor_res[i,3] <- test$p.value
}
names(CHMP4C_cor_res) <- c("symbol","correlation","pvalue")
View(CHMP4C_cor_res)

y <- as.numeric(GSDMA_cor_df["GSDMA",])
rownames <- rownames(GSDMA_cor_df)
GSDMA_cor_res <- data.frame(rownames)
for (i in 1:length(rownames)){
  test <- cor.test(as.numeric(GSDMA_cor_df[i,]),y,type="spearman")
  GSDMA_cor_res[i,2] <- test$estimate
  GSDMA_cor_res[i,3] <- test$p.value
}
names(GSDMA_cor_res) <- c("symbol","correlation","pvalue")
View(GSDMA_cor_res)

library(ggstatsplot)

for (i in 1:length(mul_cox_gene)){
  gene <- mul_cox_gene[i]
  print(gene)
  ggscatterstats(data = as.data.frame(t(CHMP4C_cor_df)),
                 y = get(gene),
                 x = CHMP4C,
                 centrality.para = "mean",
                 margins = "both",
                 xfill = "#CC79A7",
                 yfill = "#009E73",
                 type = "spearman",
                 marginal.type = "histogram",
                 title = paste("Relationship between CHMP4C and",i)
                 ) 
}

plot_grid(p1,p2,nrow = 1,labels = LETTERS[1:2])
ggscatterstats(data = as.data.frame(t(CHMP4C_cor_df)),
               y = SCEL,
               x = CHMP4C,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               type = "spearman",
               marginal.type = "histogram",
               title = paste("Relationship between CHMP4C and",i)
)


ggscatterstats(data = as.data.frame(t(GSDMA_cor_df)),
               y = SCEL,
               x = GSDMA,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               type = "spearman",
               marginal.type = "histogram",
               title = paste("Relationship between CHMP4C and",i)
)


# 20 Chemotherapy sensitivity--------------------------------------------------------------
setwd("/Users/tielinwu/Projects/BLCA")
if (!dir.exists("./20_Medicinal_Sensity")){
  dir.create("./20_Medicinal_Sensity")
}
setwd("./20_Medicinal_Sensity")

library("pRRophetic")

library(ggplot2)

set.seed(12345)

model_expr<-fpkm_tumor[mul_cox_gene, rownames(risk)]

riskscore<-data.frame(rownames(risk),risk$risk)

colnames(riskscore)<-c('sample','risk')

riskscore$risk[which(riskscore$risk=="Low")] <-'Low risk'

riskscore$risk[which(riskscore$risk=="High")] <-'High risk'
head(riskscore)
drug<-read.table('../data/drugs.txt',sep='\t',header=F)

ic50<-data.frame(riskscore$sample)
a<-data.frame(row.names=riskscore$sample,riskscore$risk)

colnames(a)<-'risk'

cnt<-1

while (cnt < 139) {
  
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  
  Tipifarnib<-data.frame(predictedPtype)
  
  colnames(Tipifarnib)<-drug[cnt,]
  
  a<-cbind(a,Tipifarnib)
  
  cnt = cnt + 1
}

write.table(a,'IC50_2.xls',sep='\t',quote=F)

b<-a
b[b<0]<-NA
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}

c<-removeColsAllNa(b)

na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
medicinal_result <- t(subset(x, select = -risk)) 
high_group <- rownames(risk)[which(risk$risk=="High")]
low_group <- rownames(risk)[which(risk$risk=="Low")]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                             1,
                             mean), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                             1, 
                             mean), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test2.xls",
            quote = F,
            row.names = F)

rTable2 <- rTable[which(rTable$high_group_res > rTable$low_group_res &
                          rTable$padj < 0.05),]
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, padj=rTable$padj)
drugs_res <- drugs_res[which(rownames(drugs_res) %in% 
                               c("Cisplatin", "Vinblastine",
                                 "Methotrexate","Gemcitabine","Doxorubicin")),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","padj"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High risk", "Low risk") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High risk", "Low risk"))
violin_dat$score <- log(violin_dat$score, base = exp(1))

drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                              color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                              add = "jitter",
                              short.panel.labs = T,
                              ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                         facet.by = "drugs",
                         short.panel.labs = T,
                         panel.labs.background = list(fill = "white"),
                         ncol = 3,
                         scales = "free_y") + xlab("") + ylab("Estimated ln(IC50)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "drugs_Cisplatin_boxplot.pdf", height = 5, width = 4)
ggsave(filename = "drugs_Cisplatin_boxplot.png", height = 5, width = 4)


library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
rTable3 <- rTable[which(rTable$high_group_res < rTable$low_group_res &
                          rTable$padj < 0.05),]
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, padj=rTable$padj)
drugs_res <- drugs_res[rownames(rTable3),]

violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","padj"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High risk", "Low risk") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High risk", "Low risk"))
head(violin_dat)
violin_dat$score <- log(violin_dat$score, base = exp(1))

drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 3,
                           scales = "free_y") + xlab("") + ylab("Estimated ln(IC50)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "drugs_Cisplatin_boxplot.pdf", height = 6, width = 5)
ggsave(filename = "drugs_Cisplatin_boxplot.png", height = 6, width = 5)


drugs_box <- ggboxplot(violin_dat,
                      x = "indivs",
                      y = "score",
                      fill = "indivs",
                      palette =c("#A73030FF", "#0073C2FF")) +
  stat_compare_means() +
  theme_bw() +
  labs(title = "", x = "", y = "Cisplatin Estimated ln(IC50)") +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(colour="black",face="bold",size=15), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.y=element_text(size=18,face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
drugs_box
ggsave(filename = "Cisplatin_boxplot_20220527.png", height = 5, width = 4,drugs_box)
ggsave(filename = "Cisplatin_boxplot_20220527.pdf", height = 5, width = 4,drugs_box)
