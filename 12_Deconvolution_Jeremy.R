suppressPackageStartupMessages({
library(xbioc)
library(MuSiC)
library(reshape2)
library(edgeR)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(Seurat)
})

dir.create("deconvolution")

# Load CMC data and create ExpSet
load("input/Expression_Input.RData")
load("input/NAc_Rat_JeremyDay.robj")

NAc <- subset(x = All_Groups_log, subset = (Stim == "Saline"))

# Exp Set for Allen
phenoSC <- NAc@meta.data %>% as.data.frame()

phenoSC$CellID <- as.factor(NAc@active.ident)
phenoSC$Sample_ID <- rownames(phenoSC)
phenoDataSC <- new("AnnotatedDataFrame",data=phenoSC)

exprSC <- abs(as.matrix(GetAssayData(object = NAc[["RNA"]],slot = "data")))
exprSC <- log2(exprSC+1)

scSet <- ExpressionSet(assayData=exprSC, phenoData=phenoDataSC)

# Getting D1 D2 ready
pdD1D2 <- metadata %>%
			rownames_to_column("ID") %>%
			filter(Genotype == "WT") %>%
			select(ID, Cell) %>%
			column_to_rownames("ID")


expD1D2 <- log2(as.matrix(exp[,colnames(exp)%in%rownames(pdD1D2)])+1)



phenoData <- new("AnnotatedDataFrame",data=pdD1D2)
D1D2data <- ExpressionSet(assayData=expD1D2, phenoData=phenoData)

D1D2_Deconv = music_prop(
	bulk.eset = D1D2data, 
	sc.eset = scSet,
	clusters = 'CellID', 
	samples = 'Sample_ID',
	iter.max = 100,
	nu = 1e-10)

save(D1D2_Deconv,file = "deconvolution/D1D2_Deconv_Jeremy.RData")
Prob <- as.data.frame(D1D2_Deconv$Est.prop.weighted)
Prob$Rows <- rownames(Prob)
tmp <- melt(Prob)
tmp$variable <- as.character(tmp$variable)
df <- merge(tmp,pdD1D2,by.x="Rows",by.y="row.names",all=T)
df <- df %>%
		arrange(variable) 

pdf("deconvolution/Deconvolution_JeremyDay.pdf",width=8,height=5)
ggboxplot(df, "variable", "value",
color = "Cell")+
ggtitle("Cell Type Proportion")+
xlab("")+ 
ylab("Proportion")+
theme(legend.position="right")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
ylim(0,1)
dev.off()

