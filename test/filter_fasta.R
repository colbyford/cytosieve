library(readr)
library(dplyr)
library(Biostrings)

pv01_genes <- readDNAStringSet("data/PVP01_Genes.fasta")

profile <- read_delim("data/profile_41467_2019_8312_MOESM7_ESM.txt", "\t")

profile_filtered <- profile %>% filter(Schizonts > 1) %>% select(GeneID)
# profile_filtered <- profile %>% filter(Schizonts == 1) %>% select(GeneID)

pv01_genes_filtered <- pv01_genes[profile_filtered$GeneID]

writeXStringSet(pv01_genes_filtered, "data/PVP01_Genes_Schizont_gt1.fasta", format="fasta")
