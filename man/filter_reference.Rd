% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_reference.R
\name{filter_reference}
\alias{filter_reference}
\title{Filter GFF Reference Files}
\usage{
filter_reference(gff, genes_names = NULL, eliminate_matches = TRUE,
  verbose = TRUE)
}
\arguments{
\item{gff}{gffObject. Object from the read_gff function.}

\item{genes_names}{Vector of strings. List of gene IDs by which to filter.}

\item{eliminate_matches}{Logical. Should the matches that are found be filtered out (TRUE) or should non-matched be filtered out (FALSE)? (Default: TRUE)}

\item{verbose}{Logical. Should the process print out which read and gene it is currently searching? (Default: TRUE)}
}
\description{
Filter out genes from a GFF file that are not expressed in a particular stage/cell type of interest.
Works by filtering out lines with genes that match based on a reference set of gene names.
}
\examples{
gff_filtered <- filter_reference(gff, gene_names = c("geneA", "geneB"))
}
