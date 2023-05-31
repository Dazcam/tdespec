#' Get vector of hg38 genes overlapping MHC region using biomaRt
#'
#' @return A vector of unique genes that fall within MHC region of chr6.
#' @export
#'
#' @param start Start position of MHC region.
#' @param end End position of MHC region.
#' @param server Ensembl server to use to call genes.
#'
#' @examples get_mhc_genes(start = 28510120, end = 33480577, server = "ensembl")
get_mhc_genes <- function(

    start = "28510120",
    end = "33480577",
    server = "ensembl"

    ) {

  cat('\n\nConnecting to Ensembl server ...')
  mart <- biomaRt::useMart(server)
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)

  cat('\nCalling MHC region genes ...\n')
  mhc_genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                     filters = c("chromosome_name","start","end"),
                     values = list(chromosome = "6", start = start, end = end,
                     mart = mart))
  mhc_genes_uniq <- stringi::stri_remove_empty(unique(mhc_genes$hgnc_symbol), na_empty = FALSE)
  cat('\n\nMHC genes detected:', length(mhc_genes_uniq), '\n\n')

  return(mhc_genes_uniq)

}
