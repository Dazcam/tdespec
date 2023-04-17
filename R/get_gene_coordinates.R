#' Get gene co-ordinates from MAGMA hg19 file with genes overlapping MHC region removed
#'
#' @return A set of gene co-ordinates
#' @export
#'
#' @examples
get_gene_coordinates <- function() {

  gene_coordinates <- readr::read_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.gene.loc.txt'),
                                      col_names = FALSE, col_types = 'cciicc') %>%
    dplyr::rename(entrez = "X1", chr = "X2", start = 'X3',
                  end = 'X4', strand = "X5", hgnc = 'X6') %>%
    dplyr::filter(!hgnc %in% mhc_genes_uniq) %>%
    #write_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.MHCremoved.gene.loc.txt'), col_names = FALSE) %>%
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::select(chr, start, end, entrez, hgnc)

  return(gene_coordinates)

}