#############################
#' @title Filter Reads
#' @name filter_reads
#' @param input_path String. Path to the input .fastq file.
#' @param output_path String. Desired path to the output filtered .fastq file.
#' @param genes_to_find XString set of gene sequences for which to search. Usually longer sequences than the individual reads. Usually from a FASTA file.
#' @param eliminate_matches Logical. Should the matches that are found be filtered out (TRUE) or should non-matched be filtered out (FALSE)? (Default: TRUE)
#' @param pct_variability Numerical. The amount of variability that will be allowed in considering a match. (Default: 0.10)
#' @param paired Logical. Should the filtering occur on both a paired forward (R1) and reverse (R2) read?
#' If TRUE, then reads that are filtered from the R1 file will also be removed from the R2 file.  (Default: FALSE)
#' @param input_r2_path String. Path to the input (reverse) .fastq file. (Optional, only used when paired = TRUE)
#' @param output_r2_path String. Desired path to the output (reverse) filtered .fastq file. (Optional, only used when paired = TRUE)
#' @param par_method String. Either "multicore" for single node parallelism or "foreach" for foreach-based distribution.
#' If "foreach" is used, then this tasks will be distributed according to the doParallel, doMPI, doSnow, etc. specifications. (Default = "multicore")
#' @param verbose Logical. Should the process print out which read and gene it is currently searching? (Default: TRUE)
#'
#' @description Filter out reads that have genes that are not expressed in a particular stage/cell type of interest.
#' Generates a FASTQ file by filtering out reads based on a reference set of gene sequences.
#' @export filter_reads

filter_reads <- function(input_path,
                         output_path,
                         genes_to_find = NULL,
                         eliminate_matches = TRUE,
                         pct_variability = 0.10,
                         paired = FALSE,
                         input_r2_path = NULL,
                         output_r2_path = NULL,
                         par_method = "multicore",
                         verbose = TRUE){

  ext <- strsplit(basename(input_path), split="\\.")[[1]][2]
  if(ext != "fastq"){
    stop("Your input file is not in the FASTQ format.\nCurrently, this function only supports .fastq input files.")
  }

  fastq_orig <- ShortRead::readFastq(input_path)

  ## Intialize list for which reads to keep (FALSEs) or remove (TRUEs)
  match_list <- c(rep(FALSE, length(fastq_orig)))

  ## Find the read segment in the reference gene(s)
  search_read <- function(x, genes = genes_to_find, percentage = pct_variability){

    if(verbose){cat("Searching for matches in read: ", x, "\n")}

    read <- fastq_orig@sread[x] %>% as.character()

    mismatch_threshold <- (nchar(read) * percentage) %>% as.integer()

    for (g in seq_along(genes)){
      gene <- genes[g] %>% as.character()

      if(verbose){cat("\tSearching filter gene: ", names(gene), "\n")}

      read_matches <- Biostrings::matchPattern(read,
                                               gene %>%
                                                 Biostrings::DNAString(),
                                               with.indels=TRUE,
                                               max.mismatch = mismatch_threshold)

      if (length(read_matches) > 0){
        cat("\tFound a match!: [ Read:", x, ", Gene:", names(gene),"]\n")
        # match_list[x] = TRUE
        # return(TRUE)
        match_found <- TRUE
        ## Cease the gene search on this read and go to the next
        break
      } else {
        match_found <- FALSE
        # return(FALSE)
      }

      rm(read_matches)

    }

    if(match_found){
      return(TRUE)
    } else {
      return(FALSE)
    }

  }

  if (par_method == "multicore"){
    match_list <- mclapply(seq_along(fastq_orig@sread),
                           search_read,
                           mc.cores = if (.Platform$OS.type == "windows"){
                             1
                           } else {
                             parallel::detectCores()-1
                           }) %>% unlist()

  } else if (par_method == "foreach"){
    match_list <- foreach (row = seq_along(fastq_orig@sread), .combine = rbind, .packages = "cytosieve") %dopar% {
      search_read(row,
                  genes = genes_to_find,
                  percentage = pct_variability)
    }

  } else {
    stop("One of `multicore` or `foreach` parallelization methods must be selected.")
  }


  # match_list <- (match_list + match_list_iter) %>%
  #   as.logical()

  ## Loop through genes
  # for (g in seq_along(genes_to_find)){
  #   gene <- genes_to_find[g] %>% as.character()
  #
  #   cat("Current filter gene: ", names(gene), "\n")
  #
  #   match_list_iter <- mclapply(seq_along(fastq_orig@sread),
  #                               search_read,
  #                               mc.cores = if (.Platform$OS.type == "windows"){
  #                                 1
  #                               } else {
  #                                 parallel::detectCores()-1
  #                               }) %>% unlist()
  #
  #   match_list <- (match_list + match_list_iter) %>%
  #     as.logical()
  #
  # }

  ## Fix any skipped reads at the end
  # length(match_list) <- length(fastq_orig)
  # match_list[is.na(match_list)] <- FALSE

  which(match_list)

  ## Define filtering function
  if(eliminate_matches){
    ##  Filter out matching reads
    fun <- function(x) {
      x[-match_list]
    }
  } else {
    ##  Filter out non-matching reads
    fun <- function(x) {
      x[match_list]
    }
  }

  cat("***Found",
      length(which(match_list)), "match(es) out of",
      length(fastq_orig@sread), "reads by searching across",
      length(genes_to_find), "genes.***\n")

  ## Filter FASTQ(s)
  fastq_filtered <- fastq_orig[-which(match_list)]
  ShortRead::writeFastq(fastq_filtered, file=output_path, compress=FALSE)
  # filterFastq(input_path, output_path, filter=fun, compress=FALSE)

  if(paired){
    fastq_r2_orig <- ShortRead::readFastq(input_r2_path)
    fastq_r2_filtered <- fastq_r2_orig[-which(match_list)]
    ShortRead::writeFastq(fastq_r2_filtered, file=output_r2_path, compress=FALSE)
    # filterFastq(input_r2_path, output_r2_path, filter=fun, compress=FALSE)
    }

  # return(which(match_list))
}
