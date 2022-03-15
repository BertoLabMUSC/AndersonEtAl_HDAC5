suppressPackageStartupMessages({
library(biomaRt)
library(tidyverse)
})

# Load Data
load("d1_dge/D1_Dge_Data.RData")

# Convert to human 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = D1_DGE$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(D1_DGE,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(logFC > 0 ~ "D1_Upreg", logFC < 0  ~ "D1_Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("D1_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df,"d1_dge/D1_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df <- D1_DGE %>%
                mutate(Direction = case_when(logFC > 0 ~ "D1_Upreg", logFC < 0  ~ "D1_Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp <- data.frame(Gene = df$Gene, Direction = rep("D1_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df,"d1_dge/D1_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)

# D2
load("d2_dge/D2_Dge_Data.RData")

# Convert to human 
MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = D2_DGE$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(D2_DGE,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(logFC > 0 ~ "D2_Upreg", logFC < 0  ~ "D2_Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("D2_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df,"d2_dge/D2_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df <- D2_DGE %>%
                mutate(Direction = case_when(logFC > 0 ~ "D2_Upreg", logFC < 0  ~ "D2_Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp <- data.frame(Gene = df$Gene, Direction = rep("D2_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df,"d2_dge/D2_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)

