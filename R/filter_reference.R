#############################
#' @title Filter GFF Reference Files
#' @name filter_reference
#' @param gff gffObject. Object from the read_gff function.
#' @param genes_names Vector of strings. List of gene IDs by which to filter.
#' @param eliminate_matches Logical. Should the matches that are found be filtered out (TRUE) or should non-matched be filtered out (FALSE)? (Default: TRUE)
#' @param verbose Logical. Should the process print out which read and gene it is currently searching? (Default: TRUE)
#'
#' @description Filter out genes from a GFF file that are not expressed in a particular stage/cell type of interest.
#' Works by filtering out lines with genes that match based on a reference set of gene names.
#' @export filter_reference
#' @examples
#' gff_filtered <- filter_reference(gff, gene_names = c("geneA", "geneB"))

filter_reference <- function(gff,
                             genes_names = NULL,
                             eliminate_matches = TRUE,
                             verbose = TRUE){

  # ext <- strsplit(basename(input_path), split="\\.")[[1]][2]
  # if(ext != "gff"){
  #   stop("Your input file is not in the GFF format.\nCurrently, this function only supports .gff input files.")
  # }

  if(class(gff)!= "gffObject"){
    stop("Your input object is not in the class `gffObject`. Hint: Use the `read_gff` function to make this object from a .gff file.")
  }

  if (eliminate_matches){
    gff$body <- gff$body %>%
      mutate(ID = stringr::str_extract(attributes, "(?<=ID=)([A-Za-z0-9\\-\\_\\.]+)(?=;)")) %>%
      filter(!ID %in% gene_names) %>%
      select(-ID)
  } else {
    gff$body <- gff$body %>%
      mutate(ID = stringr::str_extract(attributes, "(?<=ID=)([A-Za-z0-9\\-\\_\\.]+)(?=;)")) %>%
      filter(ID %in% gene_names) %>%
      select(-ID)
  }


  return(gff)

}
