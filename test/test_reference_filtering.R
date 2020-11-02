library(cytosieve)
library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)

## Genes that are NOT expressed in Schizonts
genes <- readDNAStringSet("data/Pvivax_schizont_filter_all29.fasta")[1:5]

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

gff_header <- readLines(input_gff_path)
gff_header <- gff_header[startsWith(gff_header, "##")]


gff_body <- ape::read.gff(input_gff_path)

gff_body_filtered <- gff_body %>%
  mutate(ID = stringr::str_extract(attributes, "(?<=ID=)([A-Za-z0-9\\-\\_\\.]+)(?=;)")) %>%
  filter(ID %in% gene_names) %>%
  select(-ID)


read_gff <- function(x, path){

}


write_gff <- function(x, path){

}
