#' Get vector of MHC genes using biomaRt
#'
#' @return A vector of unique genes that fall within MHC region of chr6.
#' @export
#'
#' @examples
get_mhc_genes <- function() {

  mart <- biomaRt::useMart("ensembl")
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
  mhc_genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                     filters = c("chromosome_name","start","end"),
                     values = list(chromosome = "6", start = "28510120", end = "33480577"),
                     mart = mart)
  mhc_genes_uniq <- stringi::stri_remove_empty(unique(mhc_genes$hgnc_symbol), na_empty = FALSE)
  cat('\n\nMHC genes:', length(mhc_genes_uniq), '\n')

  return(mhc_genes_uniq)

}
