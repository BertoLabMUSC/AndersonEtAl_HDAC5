suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(UpSetR)
library(RColorBrewer)
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(clusterProfiler)
library(RRHO2)
})

dir.create("Combined_Data")

d1 <- read.table("d1_dge/D1_Dge_MouseID.txt",header=T)
d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)

df <- rbind(d1,d2)

df <- df %>% 
		filter(Direction %in% c("D1_Upreg","D1_Downreg","D2_Upreg","D2_Downreg"))


pdf("Combined_Data/Intersection_D1_D2.pdf",width=6,height=4,useDingbats=FALSE)
l <- split(as.character(df$Gene),df$Direction)
Waves <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Waves, ToTGene))
names(metadata) <- c("Genotype", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), list(type = "matrix_rows", 
    column = "sets", colors = c(D1_Upreg = "blue", D1_Downreg = "red", D2_Upreg = "orange",D2_Downreg = "purple"), 
    alpha = 0.5))))
dev.off()

# database
load("d1_dge/D1_Dge_Data.RData")
load("d2_dge/D2_Dge_Data.RData")

d1 <- D1_FullTab %>%
	  select(Gene,logFC,adj.P.Val) %>%
	  rename(D1_logFC = logFC, D1_FDR = adj.P.Val)


d2 <- D2_FullTab %>%
	  select(Gene,logFC,adj.P.Val) %>%
	  rename(D2_logFC = logFC, D2_FDR = adj.P.Val)

d1_d2_comb <- Reduce(dplyr::full_join, list(d1,d2)) %>%
		 		as.data.frame()


d1_d2_comb_sign <- d1_d2_comb %>%
					filter(D1_FDR < 0.05, D2_FDR < 0.05, abs(D1_logFC) > 0.3, abs(D2_logFC) > 0.3, sign(D2_FDR) == sign(D1_FDR))

save(d1_d2_comb,d1_d2_comb_sign,file = "Combined_Data/D1_D2_Combined_Stats.RData")

openxlsx::write.xlsx(d1_d2_comb, 
                     file = "Combined_Data/D1_D2_Combined_Stats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=TRUE)

openxlsx::write.xlsx(d1_d2_comb_sign, 
                     file = "Combined_Data/D1_D2_Combined_DGE.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=TRUE)


pdf("Combined_Data/Scatter_D1_D2_Combined_DGE.pdf",width=4,height=4,useDingbats=FALSE)

Genes <- d1_d2_comb_sign %>%
		filter(Gene %in% c("Slc17a7","Akap5","Kcnv1"))


ggscatter(d1_d2_comb_sign, 
            x = "D1_logFC", 
            y = "D2_logFC",
            color = "red",
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("D1 log2(Fold Change)")+ 
      ylab("D2 log2(Fold Change)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
      ylim(-2,+2) + xlim(-2,+2) + 
      geom_text_repel(data = Genes, 
                      mapping = aes(label = Gene), 
                      size = 5,
    				  nudge_x = .15,
    				  nudge_y = .5,
   					  segment.linetype = 6,
    				  segment.curvature = -1e-20,
    				  arrow = arrow(length = unit(0.015, "npc")),
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))
dev.off()

# Venn by sign

d1 <- read.table("d1_dge/D1_Dge_MouseID.txt",header=T)
d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)
df <- rbind(d1,d2)

df <- df %>% 
		filter(Direction %in% c("D1_Upreg","D1_Downreg","D2_Upreg","D2_Downreg"))

list <- split(df, df$Direction)

x <- list(D1_Downreg=as.character(list[[1]]$Gene),D1_Upreg=as.character(list[[2]]$Gene),D2_Downreg=as.character(list[[3]]$Gene),D2_Upreg=as.character(list[[4]]$Gene))

ven <- venn.diagram(x= list(D1_Downreg=as.character(list[[1]]$Gene),
							D1_Upreg=as.character(list[[2]]$Gene),
							D2_Downreg=as.character(list[[3]]$Gene),
							D2_Upreg=as.character(list[[4]]$Gene)),
filename = NULL,
fill = c("red","grey60","magenta","cyan"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2)

pdf(file="Combined_Data/Venn_DGEs.pdf",width=5,height=5)
grid.draw(ven)
dev.off()

# All
d1 <- read.table("d1_dge/D1_Dge_MouseID.txt",header=T)
d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)
df <- rbind(d1,d2)

df <- df %>% 
		filter(Direction %in% c("D1_All","D2_All"))

list <- split(df, df$Direction)

x <- list(D1_All=as.character(list[[1]]$Gene),D2_All=as.character(list[[2]]$Gene))

ven <- venn.diagram(x= list(D1=as.character(list[[1]]$Gene),
							D2=as.character(list[[2]]$Gene)),
filename = NULL,
fill = c("red","cyan"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2)

pdf(file="Combined_Data/Venn_DGEs_All.pdf",width=5,height=5)
grid.draw(ven)
dev.off()

# Gene Ontology
load("Combined_Data/D1_D2_Combined_Stats.RData")
GOI <- bitr(as.character(d1_d2_comb_sign$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("Combined_Data/GeneOnto_D1_D2_Combined.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(GeneOnto), 
                     file = "Combined_Data/GeneOnto_D1_D2_Combined.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

pdf("Combined_Data/GeneOnto_D1_D2_Combined_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(GeneOnto, categorySize="pvalue")
dev.off()

# For D1 and D2

d1 <- read.table("d1_dge/D1_Dge_MouseID.txt",header=T)

d1 <- d1 %>% filter(Direction == "D1_All")

GOI <- bitr(as.character(d1$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D1_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("d1_dge/GeneOnto_D1.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D1_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D1_GeneOnto), 
                     file = "d1_dge/GeneOnto_D1.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

pdf("d1_dge/GeneOnto_D1_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(D1_GeneOnto, categorySize="pvalue")
dev.off()


d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "D2_All")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("d2_dge/GeneOnto_D2.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "d2_dge/GeneOnto_D2.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

pdf("d2_dge/GeneOnto_D2_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(D2_GeneOnto, categorySize="pvalue")
dev.off()

# Split by signature
d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "D2_Upreg")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("d2_dge/GeneOnto_D2_Upreg.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "d2_dge/GeneOnto_D2_Upreg.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

pdf("d2_dge/GeneOnto_D2_Upreg_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(D2_GeneOnto, categorySize="pvalue")
dev.off()


d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "D2_Downreg")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("d2_dge/GeneOnto_D2_Downreg.pdf",width=8,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "d2_dge/GeneOnto_D2_Downreg.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

pdf("d2_dge/GeneOnto_D2_Downreg_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(D2_GeneOnto, categorySize="pvalue")
dev.off()

# RRHO2
load("d1_dge/D1_Dge_Data.RData")
load("d2_dge/D2_Dge_Data.RData")

d1_up <- D1_DGE %>% filter(logFC > 0) %>% select(Gene, P.Value) %>% mutate(DDE = -log10(P.Value)*1) %>% select(Gene, DDE)
d1_down <- D1_DGE %>% filter(logFC < 0) %>% select(Gene, P.Value) %>% mutate(DDE = -log10(P.Value) * -1) %>% select(Gene, DDE)
d1 <- rbind(d1_up, d1_down)

d2_up <- D2_DGE %>% filter(logFC > 0) %>% select(Gene, P.Value) %>% mutate(DDE = -log10(P.Value)* 1) %>% select(Gene, DDE)
d2_down <- D2_DGE %>% filter(logFC < 0) %>% select(Gene, P.Value) %>% mutate(DDE = -log10(P.Value) * -1) %>% select(Gene, DDE)
d2 <- rbind(d2_up, d2_down)

d1 <- d1 %>% filter(Gene %in% intersect(d1$Gene, d2$Gene))
d2 <- d2 %>% filter(Gene %in% intersect(d1$Gene, d2$Gene))

RRHO_obj <-  RRHO2_initialize(d1, d2, labels = c("D1", "D2"), log10.ind=TRUE,method="fisher")

# Combined
d1 <- read.table("d1_dge/D1_Dge_MouseID.txt",header=T)
d2 <- read.table("d2_dge/D2_Dge_MouseID.txt",header=T)
df <- rbind(d1,d2)

df <- df %>% 
      filter(Direction %in% c("D1_All","D2_All"))

list <- split(df, df$Direction)

a <- intersect(list[[1]]$Gene, list[[2]]$Gene)
b <- list[[1]][!(list[[1]]$Gene %in% a),]
c <- list[[2]][!(list[[2]]$Gene %in% a),]

a <- data.frame(Gene = a, Direction = rep("Overlap",length(a)))
b <- data.frame(Gene = b$Gene, Direction = rep("D1_Specific",nrow(b)))
c <- data.frame(Gene = c$Gene, Direction = rep("D2_Specific",nrow(c))) 

Comb <- rbind(a,b,c)

write.table(Comb,"Combined_Data/Combined_Genes.txt",sep="\t",quote=F,row.names=F)

openxlsx::write.xlsx(Comb, 
                     file = "Combined_Data/Combined_Genes.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats", 
                     overwrite = TRUE)


human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = Comb$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)


signComb <- merge(Comb,MGI,by.x="Gene",by.y="MGI.symbol",all=F) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(Direction)

write.table(signComb,"Combined_Data/Combined_Genes_HumanID.txt",sep="\t",quote=F,row.names=F)

openxlsx::write.xlsx(signComb, 
                     file = "Combined_Data/Combined_Genes_HumanID.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats", 
                     overwrite = TRUE)


