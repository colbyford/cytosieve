library(cytosieve)
library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)

## Genes that are NOT expressed in Schizonts
genes <- readDNAStringSet("data/Pvivax_schizont_filter_all29.fasta")

gene_name_pattern = "PVP01_[0-9]+"
gene_names <- names(genes) %>% str_extract(gene_name_pattern)


## Append all sequence files as a single fasta
# # txtfiles <- list.files(path = "201019EN-032/", pattern="*.txt")
# txtfiles <- paste0("test/schizont_test/reference_filtering/", list.files(path = "test/schizont_test/reference_filtering/", pattern="*.fasta"))
#
# # outFile <- file("output/all_reads.fasta", "w")
# outFile <- file("test/schizont_test/reference_filtering/Pvivax_ASM241v2.fasta", "w")
#
# for (i in txtfiles){
#   x <- readLines(i)
#   writeLines(x, outFile)
# }
#
# close(outFile)

# input_fasta_path <- "test/schizont_test/reference_filtering/Pvivax_ASM241v2.fasta"
# output_fasta_path <- paste0(str_remove(input_fasta_path, ".fasta"), "_filtered.fasta")
#
# fasta_orig <- ShortRead::readFasta(input_fasta_path)

input_gff_path <- "test/schizont_test/reference_filtering/PlasmoDB-48_PvivaxP01.gff"
output_gff_path <- paste0(str_remove(input_gff_path, ".gff"), "_filtered.gff")

gff <- read_gff(input_gff_path)

gff_filtered <- filter_reference(gff, gene_names)

write_gff(gff_filtered, path = output_gff_path, append = FALSE)
