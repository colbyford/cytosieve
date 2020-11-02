#############################
#' @title Filter Reference Files
#' @name filter_reference
#' @param input_path String. Path to the input .fastq file.
#' @param output_path String. Desired path to the output filtered .fastq file.
#' @param genes_to_find XString set of gene sequences for which to search. Usually longer sequences than the individual reads. Usually from a FASTA file.
#' @param eliminate_matches Logical. Should the matches that are found be filtered out (TRUE) or should non-matched be filtered out (FALSE)? (Default: TRUE)
#' @param verbose Logical. Should the process print out which read and gene it is currently searching? (Default: TRUE)
#'
#' @description Filter out reads that have genes that are not expressed in a particular stage/cell type of interest.
#' Generates a FASTQ file by filtering out reads based on a reference set of gene sequences.
#' @export filter_reference

filter_reference <- function(input_path,
                             output_path,
                             genes_to_find = NULL,
                             eliminate_matches = TRUE,
                             verbose = TRUE){

  ext <- strsplit(basename(input_path), split="\\.")[[1]][2]
  if(ext != "fasta"){
    stop("Your input file is not in the FASTQ format.\nCurrently, this function only supports .fastqainput files.")
  }

  fasta_orig <- ShortRead::readFasta(input_path)
  ## GFF


}
