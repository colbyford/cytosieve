# *cytosieve* - An R package for _in silico_ cytometric filtering of bulk sequencing data using transcriptional profiles.
<h3 align = "right">Colby T. Ford, Ph.D., Daniel Kepple, Daniel Janies, Ph.D., and Eugenia Lo, Ph.D.</h3>

<img align="right" src="https://raw.githubusercontent.com/colbyford/cytosieve/main/img/cytosieve_hex.png" alt="cytosieve icon" width="200">

## Motivation
Bulk sequencing, as compared to single cell sequencing, is a more cost-effective method of sequencing cells and offers useful insights into comparative transcriptomics. However, bulk sequencing is insufficient for studies that need the information of a single cell type among heterogeneous assortments of cells in a single sample.

Using transcriptional profiles, like those used in tools like [CIBERSORTx](https://cibersortx.stanford.edu/), a list of exclusion or inclusion genes by which to filter can be defined. Using this list of a genes, we can remove a given read if the read is part of a gene that is not expressed in the cell type of interest.

<img align="right" src="https://raw.githubusercontent.com/colbyford/cytosieve/main/img/process.png" width=500>

--------------------------

## Installation
```r
# install.packages("devtools")
devtools::install_github("colbyford/cytosieve")
```

## Usage
This package will look at each read in the input FASTQ file (R1) and will find any matches across the genes in your filtering set. If a match is found, the read is removed as it is determined to be from an undesired cell type. The following code will output filtered versions of the input FASTQ files, better purifying the reads to only those from the desired types of cells.

```r
library(cytosieve)
library(Biostrings)
library(stringr)

## Define genes that are NOT expressed in the cell type of interest
genes <- readDNAStringSet("data/filter_genes.fasta") # The exclusion list of sequences

## A paired set of sample reads (can also be run on non-paired samples)
input_R1_fastq_path <- "data/sample_r1_reads.fastq"
output_R1_fastq_path <- paste0(str_remove(input_R1_fastq_path, ".fastq"), "_filtered.fastq")
input_R2_fastq_path <- "data/sample_r2_reads.fastq"
output_R2_fastq_path <- paste0(str_remove(input_R2_fastq_path, ".fastq"), "_filtered.fastq")

## Filter out reads from the input R1 (and R2) sample files
filter_reads(input_path = input_R1_fastq_path,
             output_path = output_R1_fastq_path,
             genes_to_find = genes,                 # Genes that do not occur in the desired cell type
             eliminate_matches = TRUE,              # Remove genes that match from the exclusion list
             pct_variability = 0.10,
             paired = TRUE,                         # Optional
             input_r2_path = input_R2_fastq_path,   # Optional
             output_r2_path = output_R2_fastq_path, # Optional
             verbose = FALSE)
```
