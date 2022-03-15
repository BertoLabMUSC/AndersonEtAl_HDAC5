# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir Combined_Data/enrichments_combined/

cp utils/geneset/*.RData Combined_Data/enrichments_combined/
cp utils/Enrichment.r Combined_Data/enrichments_combined/
cp Combined_Data/Combined_Genes_HumanID.txt Combined_Data/enrichments_combined/
cp Combined_Data/Combined_Genes.txt Combined_Data/enrichments_combined/

cd Combined_Data/enrichments_combined/

mkdir STATS/

# scRNA healty
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_Markers -W 4 -H 4
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l GeneSets_BBLake_2018.RData -p -b 15585 -o STATS/BBLake_Markers -W 4 -H 4
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l GeneSets_TranEtAl_NAc.RData -p -b 15585 -o STATS/Tran_Markers -W 4 -H 6

Rscript Enrichment.r -g Combined_Genes.txt -l GeneSets_scMouse.RData -p -b 15585 -o STATS/scMouse -W 4 -H 5
Rscript Enrichment.r -g Combined_Genes.txt -l GeneSets_NAc_CellReport.RData -p -b 15585 -o STATS/scMouse_NAc -W 5 -H 6
Rscript Enrichment.r -g Combined_Genes.txt -l GeneSets_IEG.RData -p -b 15585 -o STATS/IEGs -W 5 -H 4

# scRNA Disorders
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell -W 4 -H 5
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell -W 4 -H 5

# Neuropsy
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 4 -H 3
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 4 -H 5
Rscript Enrichment.r -g Combined_Genes_HumanID.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari -W 4 -H 2

rm *.RData

