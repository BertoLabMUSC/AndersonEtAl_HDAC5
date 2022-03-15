# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(future.apply)
library(DESeq2)
library(pheatmap)
library(limma)
})

dir.create("d1_dge_onebatch_noout")

load("input/Expression_Input_D1.RData")

pd <- metadata_d1 %>% filter(Batch == "Batch1")
p <- exp_d1[,colnames(exp_d1) %in% rownames(pd)]

pdf("d1_dge_onebatch_noout/PCA_D1_onebatch.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette=c("black","red"), 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() 
dev.off()

p$HDAC5_D1_6 <- NULL
pd <- pd[-4,] %>% droplevels()

# DGE 
design <- model.matrix(~Genotype + sex,pd) # describe model to be fit

fitLM = lmFit(p,design,method="robust");
fitEb = eBayes(fitLM);
D1_FullTab_onebatch = topTable(fitEb, coef = "GenotypeWT",number=nrow(p)) %>%
             rownames_to_column("Gene") %>%
             mutate(logFC = -1*logFC)


D1_DGE_onebatch <- D1_FullTab_onebatch %>%
                  mutate(Abs = abs(logFC)) %>%
                  filter(adj.P.Val < 0.05 & Abs > 0.3) %>%
                  arrange(desc(Abs))


save(D1_FullTab_onebatch,D1_DGE_onebatch, file = "d1_dge_onebatch_noout/D1_Dge_Data_onebatch_noout.RData")


openxlsx::write.xlsx(D1_FullTab_onebatch, 
                     file = "d1_dge_onebatch_noout/D1_DGE_FullStats_onebatch_noout.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite = TRUE)

openxlsx::write.xlsx(D1_DGE_onebatch, 
                     file = "d1_dge_onebatch_noout/D1_DGE_onebatch_noout.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite = TRUE)


# Regression Expression
pd_sva <- pd %>%
            dplyr::select(-Genotype,-Batch) %>% #Removing ID, Region, Diagnosis from the model
            droplevels()

betas<-future_lapply(1:nrow(p), function(x)
            {
                lm(unlist(p[x,])~., data = pd_sva)
            })

residuals<-future_lapply(betas, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
p_regressed <- residuals+matrix(future_apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(p_regressed)<-rownames(p)

write.table(p_regressed,"d1_dge_onebatch_noout/expression_regressed_onebatch.txt",sep="\t",quote=F)

pdf("d1_dge_onebatch_noout/PCA_D1_AdjustedForConfound_onebatch_noout.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p_regressed))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette=c("black","red"), 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() 
dev.off()

# Input for viz Plot
df <- D1_FullTab_onebatch %>% 
        mutate(LOG = -log10(adj.P.Val), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(adj.P.Val < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & adj.P.Val < 0.05 ~ "UpReg", logFC < -0.3 & adj.P.Val < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(adj.P.Val) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
mat <- p_regressed[rownames(p_regressed)%in% top_labelled$Gene,] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("d1_dge_onebatch_noout/Boxplots_TopGenes_D1_onebatch_onebatch.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("red", "black")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("d1_dge_onebatch_noout/Vulcano_Plot_D1_onebatch_onebatch.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-5,+5)
dev.off()

# heatmap
mat <- p_regressed[rownames(p_regressed)%in% D1_DGE_onebatch$Gene,]
anno <- pd
Genotype        <- c("black", "red")
names(Genotype) <- c("WT", "HDAC5")
anno_colors <- list(Genotype = Genotype)
pdf("d1_dge_onebatch_noout/Heatmap_D1_onebatch_onebatch.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors,cluster_cols=F)
dev.off()










