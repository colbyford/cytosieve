#############################
#' @title Reading and Writing GFF Files
#' @name read_gff
#' @param path String. Path to the input/output .gff file.
#' @param na_strings Vector of strings.
#' Which strings in the GFF file are missing values and should be converted to NAs. (Default: `c(".", "?")`)
#' @param append Logical. Should the process overwrite the existing file? (Default: FALSE)
#'
#' @description Read and write General Feature Format (a.k. gene-finding format, generic feature format, GFF) is a file format for describing genes and other genetic features.
#'
#' @export read_gff
#' @export write_gff
#'
#' @examples
#' gff <- read_gff("folder/myfile.gff")
#' write_gff(gff_filtered, path = "folder/myfile_filtered.gff")

read_gff <- function(path, na_strings = c(".", "?")){
  gff_lines <- readLines(path)

  gff_header <- gff_lines[startsWith(gff_lines, "##")] %>%
    data.frame(stringsAsFactors = FALSE)

  colnames(gff_header) <- c("header")

  gff_body <- gff_lines[!startsWith(gff_lines, "##")] %>%
    strsplit("\t") %>%
    do.call(rbind.data.frame, .) %>%
    mutate_all(~ replace(., . %in% na_strings, NA)) ## Replace missing data with NAs

  colnames(gff_body) <- c("seqid", "source", "type", "start", "end",
                          "score", "strand", "phase", "attributes")


  gff <- list(header = gff_header,
              body = gff_body)

  class(gff) <- c("gffObject")

  return(gff)
}


write_gff <- function(x, path, append = FALSE){

  ## Delete file if it exists
  if (!append){
    if (file.exists(path)){
      file.remove(path)
    }
  }

  invisible(lapply(x, function(y) write.table(data.frame(y),
                                              file = path,
                                              quote = FALSE,
                                              col.names = FALSE,
                                              row.names = FALSE,
                                              append = TRUE,
                                              sep='\t')))

}
