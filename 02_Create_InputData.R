library(tidyverse)
library(data.table)
library(future.apply)
library(janitor)
plan(multicore)

setwd("input")

load("Expression_Input.RData")

metadata <- metadata %>% 
			as.data.frame() %>%
			arrange(Genotype,Cell)

exp <- exp[,match(rownames(metadata),colnames(exp))]

save(exp,metadata, file = "Expression_Input.RData")

# Filter by region
exp_d1 <- exp[,grep("D1_",colnames(exp))] 

metadata_d1 <- metadata %>% 
				filter(Cell == "D1") %>%
				select(-Cell) %>%
				arrange(Genotype)

exp_d1 <- exp_d1[,match(rownames(metadata_d1),colnames(exp_d1))]


exp_d2 <- exp[,grep("D2_",colnames(exp))] 


metadata_d2 <- metadata %>% 
				filter(Cell == "D2") %>%
				select(-Cell) %>%
				arrange(Genotype)

exp_d2 <- exp_d2[,match(rownames(metadata_d2),colnames(exp_d2))]


save(exp_d1,metadata_d1, file = "Expression_Input_D1.RData")
save(exp_d2,metadata_d2, file = "Expression_Input_D2.RData")

