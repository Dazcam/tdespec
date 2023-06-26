#' Get gene co-ordinates from MAGMA hg19 file with genes overlapping MHC region removed
#'
#' @return A set of gene co-ordinates
#' @export
#'
#' @importFrom dplyr relocate select mutate filter rename
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv write_tsv
#' @importFrom utils download.file
#' @importFrom rlang .data
#'
#' @param mhc_genes A vector of genes overlapping the MHC region.
#' @param write_ref Save a copy of the MAGMA reference to file.
#' @param ref_dir A directory to save the reference file to.
#'
#' @examples
#' data('mhc_genes')
#' get_gene_coordinates(mhc_genes)
get_gene_coordinates <- function(

  mhc_genes = NULL,
  write_ref = FALSE,
  ref_dir = NULL

  ) {

  message('\n\nDowloading MAGMA hg19 reference file ...')
  temp <- tempfile()
  utils::download.file("https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip", temp)

  if (write_ref) {

    gene_coordinates <- readr::read_tsv(unz(temp, "NCBI37.3.gene.loc"),
                                        col_names = FALSE, col_types = 'cciicc') %>%
      dplyr::rename(entrez = "X1", chr = "X2", start = 'X3',
                    end = 'X4', strand = "X5", hgnc = 'X6') %>%
    dplyr::filter(!.data$hgnc %in% mhc_genes) %>%
    readr::write_tsv(paste0(ref_dir, 'NCBI37.3.MHCremoved.gene.loc.tsv'), col_names = FALSE) %>%
    dplyr::mutate(chr = paste0("chr", .data$chr)) %>%
    dplyr::select(.data$chr, .data$start, .data$end, .data$entrez, .data$hgnc)

  } else {

    gene_coordinates <- readr::read_tsv(unz(temp, "NCBI37.3.gene.loc"),
                                        col_names = FALSE, col_types = 'cciicc') %>%
      dplyr::rename(entrez = "X1", chr = "X2", start = 'X3',
                    end = 'X4', strand = "X5", hgnc = 'X6') %>%
      dplyr::filter(!.data$hgnc %in% mhc_genes) %>%
      dplyr::mutate(chr = paste0("chr", .data$chr)) %>%
      dplyr::select(.data$chr, .data$start, .data$end, .data$entrez, .data$hgnc)

  }

  unlink(temp)

  return(gene_coordinates)

}
