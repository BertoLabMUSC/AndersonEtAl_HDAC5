library(oligo)
library(tidyverse)
library(pd.rta.1.0)
library(affycoretools)



celFiles <- list.celfiles('cels', full.names=TRUE)

rawData <- read.celfiles(celFiles)


names <- sampleNames(rawData)

meta <- read.table("cels/metadata.txt",header=T)
meta <- meta[match(names,meta$FileName),]

sampleNames(rawData) <- meta$NewFileName

rownames(meta) <- meta$NewFileName

pd <- meta %>% 
		select(Genotype,Cell,Batch,sex)

pheno <- data.frame(pData(rawData), pd)

pData(rawData) <- pheno

# Background correction (RMA)
Ethan_BkgCorrect <- backgroundCorrect(rawData)
Ethan_Normalized <- normalize(Ethan_BkgCorrect)

Ethan_Normalized <- rma(rawData, target="core")
Ethan_Normalized_Anno <- annotateEset(Ethan_Normalized, pd.rta.1.0)

# Remove genes and filter for protein coding
geneid <- pData(featureData(Ethan_Normalized_Anno))["SYMBOL"]

exp_g <- exprs(Ethan_Normalized_Anno)
rownames(exp_g) <- geneid$SYMBOL
exp_g2 <- exp_g[!is.na(rownames(exp_g)), ]

get_mean=apply(exp_g2, 1, mean)
sel_mean=tapply(1:nrow(exp_g2), rownames(exp_g2), function(x){x[which.max(get_mean[x])]})
exp_g2_2 <- exp_g2[sel_mean, ]

ensembl = biomaRt::useMart("ensembl")
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attr = c('mgi_symbol','gene_biotype')

Genes <- biomaRt::getBM(attributes = attr, filters = "mgi_symbol", values = rownames(exp_g2_2), uniqueRows = TRUE, mart = mouse)

exp_g2_2=merge(exp_g2_2,Genes, by.x="row.names", by.y="mgi_symbol", all=F)

exp_g2_3 <- exp_g2_2 %>%
			filter(gene_biotype == "protein_coding") %>%
			select(-gene_biotype) %>%
			column_to_rownames("Row.names")

pheno <- pData(Ethan_Normalized_Anno)
metadata <- pheno %>% 
			as.data.frame() %>%
			select(-index) %>%
			arrange(Genotype)

exp <- exp_g2_3[,match(colnames(exp_g2_3),rownames(metadata))]

dir.create("input")
save(exp,metadata, file = "input/Expression_Input.RData")


