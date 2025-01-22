#' @importFrom magrittr `%<>%` `%>%` `%$%`
#' @import Matrix igraph
NULL
library(magrittr)
library(Matrix)
library(igraph)

#' @export
transcripts_to_count_matrix <- function(transcripts) {
  cm <- transcripts %$% table(gene, cell)
  class(cm) <- "matrix"
  return(cm)
}


#' @export
run_pagoda_de <- function(cm, groups) {
  p2 <- pagoda2::Pagoda2$new(cm, n.cores=1, min.cells.per.gene=0, verbose=FALSE)
  de <- p2$getDifferentialGenes(
      groups=groups, z.threshold = 0,
      upregulated.only=FALSE, append.specificity.metrics=FALSE
  )

  return(de)
}