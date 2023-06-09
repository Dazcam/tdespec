#' Get gene co-ordinates from MAGMA hg19 file with genes overlapping MHC region removed
#'
#' @return A set of gene co-ordinates
#' @export
#'
#' @importFrom dplyr %>% relocate select mutate filter rename
#' @importFrom readr read_tsv write_tsv
#' @importFrom utils download.file
#' @importFrom rlang .data
#'
#' @param mhc_genes A vector of genes overlapping the MHC region.
#'
#' @examples get_gene_coordinates(c('SPI1', 'MEF2C', 'GAD1', 'TREM2', 'NEUROD1'))
get_gene_coordinates <- function(

  mhc_genes = NULL) {

  message('\n\nDowloading MAGMA hg19 reference file ...')
  temp <- tempfile()
  utils::download.file("https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip", temp)

  gene_coordinates <- readr::read_tsv(unz(temp, "NCBI37.3.gene.loc"),
                                      col_names = FALSE, col_types = 'cciicc') %>%
    dplyr::rename(entrez = "X1", chr = "X2", start = 'X3',
                  end = 'X4', strand = "X5", hgnc = 'X6') %>%
  dplyr::filter(!.data$hgnc %in% mhc_genes) %>%
  #write_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.MHCremoved.gene.loc.txt'), col_names = FALSE) %>%
  dplyr::mutate(chr = paste0("chr", .data$chr)) %>%
  dplyr::select(.data$chr, .data$start, .data$end, .data$entrez, .data$hgnc)

  unlink(temp)

  return(gene_coordinates)

}
