library(cytosieve)
library(Biostrings)
library(stringr)

## Genes that are NOT expressed in Schizonts
# genes <- readDNAStringSet("sample_data/Pvivax_schizont_filter_all29.fasta")[1:5]
# genes <- readDNAStringSet("sample_data/PVP01_Genes_Schizont_gt1.fasta")[1:5]
genes <- readDNAStringSet("sample_data/PVP01_Genes_Schizont_gt1_concatenated.fasta")


# samples <- read.table("sample_summary.txt", header = TRUE, stringsAsFactors = FALSE)
# input_R1_fastq_path <- samples$FileR1[sample]
# output_R1_fastq_path <- paste0(str_remove(input_R1_fastq_path, ".fastq"), "_filtered.fastq")
# input_R2_fastq_path <- samples$FileR2[sample]
# output_R2_fastq_path <- paste0(str_remove(input_R2_fastq_path, ".fastq"), "_filtered.fastq")

## An example paired set of P. vivax samples
# input_R1_fastq_path <- "sample_data/Pvivax_sample_r1_1000reads.fastq"
input_R1_fastq_path <- "sample_data/GHC-006_R1_001_4Msample.fastq"
output_R1_fastq_path <- paste0(str_remove(input_R1_fastq_path, ".fastq"), "_filtered.fastq")
# input_R2_fastq_path <- "sample_data/Pvivax_sample_r2_1000reads.fastq"
input_R2_fastq_path <- "sample_data/GHC-006_R2_001_4Msample.fastq"
output_R2_fastq_path <- paste0(str_remove(input_R2_fastq_path, ".fastq"), "_filtered.fastq")

# strsplit(basename(input_R1_fastq_path), split="\\.")[[1]][2]

library(doParallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)


filter_reads(input_path = input_R1_fastq_path,
             output_path = output_R1_fastq_path,
             genes_to_find = genes,
             # eliminate_matches = TRUE,
             eliminate_matches = FALSE,
             pct_variability = 0.10,
             paired = TRUE,
             input_r2_path = input_R2_fastq_path,
             output_r2_path = output_R2_fastq_path,
             # par_method = "multicore",
             par_method = "foreach",
             verbose = TRUE)

stopCluster(cl)
