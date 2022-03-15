suppressPackageStartupMessages({
library(xbioc)
library(MuSiC)
library(reshape2)
library(edgeR)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
})

dir.create("deconvolution")

# Load CMC data and create ExpSet
load("input/Expression_Input.RData")
#load("input/NAc_Rat_JeremyDay.robj")

load("input/SCE_NAc-n8_tran-etal.rda")
load("utils/geneset/GeneSets_TranEtAl_NAc.RData")
genes <- do.call(rbind,GeneSets)


sce <- subset(sce.nac.tran, , !(cellType%in%c("drop.doublet_A","drop.doublet_B", "drop.doublet_C","drop.doublet_D","drop.lowNTx")))

phenoSC <- colData(sce) %>%
        as.data.frame()

exprSC <- logcounts(sce) %>%
        as.data.frame()

exprSC_2 <- exprSC[rownames(exprSC)%in%genes$Gene,]

human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = biomaRt::getLDS(attributes = c("hgnc_symbol"), 
						filters = "hgnc_symbol", 
						values = rownames(exprSC_2), 
						mart = human, 
						attributesL = c("mgi_symbol","gene_biotype"), 
						martL = mouse, 
						uniqueRows=T)

exprSC_filt <- full_join(exprSC_2 %>%
                           mutate(Symbol = rownames(exprSC_2)),MGI, by = c("Symbol" = "HGNC.symbol")) %>%
				filter(Gene.type == "protein_coding") %>%
                dplyr::select(-Symbol,-Gene.type) %>%
                distinct(MGI.symbol, .keep_all = TRUE) %>%
                column_to_rownames("MGI.symbol")

exprSC_filt <- exprSC_filt[-which(rownames(exprSC_filt) == ""), ]


# Exp Set for Allen
phenoSC$CellID <- as.factor(phenoSC$cellType)
phenoSC$Sample_ID <- rownames(phenoSC)
phenoDataSC <- new("AnnotatedDataFrame",data=phenoSC)

scSet <- ExpressionSet(assayData=as.matrix(exprSC_filt), phenoData=phenoDataSC)

# Getting D1 D2 ready
pdD1D2 <- metadata %>%
			rownames_to_column("ID") %>%
			filter(Genotype == "WT") %>%
			select(ID, Cell) %>%
			column_to_rownames("ID")


expD1D2 <- log2(as.matrix(exp[,colnames(exp)%in%rownames(pdD1D2)])+1)

phenoData <- new("AnnotatedDataFrame",data=pdD1D2)
D1D2Set <- ExpressionSet(assayData=expD1D2, phenoData=phenoData)

save(scSet,D1D2Set, file = "deconvolution/Deconvolution_Jaffe_Markers_Input.RData")

D1D2_Deconv = music_prop(
	bulk.eset = D1D2Set, 
	sc.eset = scSet,
	clusters = 'CellID', 
	samples = 'Sample_ID',
	iter.max = 100,
	nu = 1e-10)

save(D1D2_Deconv,file = "deconvolution/D1D2_Jaffe_Markers_Deconv.RData")
Prob <- as.data.frame(D1D2_Deconv$Est.prop.weighted)
Prob$Rows <- rownames(Prob)
tmp <- melt(Prob)
tmp$variable <- as.character(tmp$variable)
df <- merge(tmp,pdD1D2,by.x="Rows",by.y="row.names",all=T)
df <- df %>%
		arrange(variable) 

pdf("deconvolution/Deconvolution_Jaffe_Markers.pdf",width=8,height=4)
ggbarplot(df, x = "variable", y = "value", 
          add = c("mean_se"),
          color = "Cell", palette = "jco",
          position = position_dodge(0.8)) +
xlab("")+ 
ylab("Proportion")+
theme(legend.position="right")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
ylim(0,0.3)
dev.off()

