#' A small scRNAseq matrix in Seurat format
#'
#' A Seurat object of 1452 cells and 23310 features (genes). It contains
#' cell type annotations at two levels of granularity that can be accessed
#' using seurat_small$cluster_level_1 and seurat_small$cluster_level_2.
#'
#' @docType data
#'
#' @usage data(seurat_small)
#'
#' @format An object of class \code{"SeuratObject"}; see \code{\link[Seurat vingette]{https://satijalab.org/seurat/articles/get_started.html}}.
#'
#' @keywords datasets
#'
#' @references Add later when paper published
#'
#' @source Add later when paper published
#'
#' @examples
#' data(seurat_small)
"seurat_small"

#' A Dataframe of genes and their genomic coordinates overlapping the MHC
#'
#' A data frame of genes extracted from ENSEMBL databses using REST API or BiomaRt.
#'
#' @docType data
#'
#' @usage data(mhc_genes)
#'
#' @format An object of class \code{"data.frame"};
#'
#' @keywords datasets
#'
#' @references Add later when paper published
#'
#' @source Add later when paper published
#'
#' @examples
#' data(mhc_genes)
"mhc_genes"

#' A Dataframe of genes and their hg19 co-ordinates which have had MHC
#' overlapping genes removed.
#'
#' A Dataframe of genes and their genomic co-ordinates taken from a gene reference
#' file which have had MHC overlapping genes removed. Entrez and HGNC IDs are
#' included.
#'
#' @docType data
#'
#' @usage data(gene_coord)
#'
#' @format An object of class \code{"data.frame"};
#'
#' @keywords datasets
#'
#' @references Downloaded from
#'
#' @source \href{https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip}{MAGMA hg19 ref}
#'
#' @examples
#' data(gene_coord)
"gene_coord"

#' A CTD object containing cell specificty scores for an snRNAseq dtatset
#'
#' A list
#'
#' @docType data
#'
#' @usage data(ctd_obj)
#'
#' @format An object of class \code{"list"};
#'
#' @keywords datasets
#'
#' @references Downloaded from
#'
#' @source Generated using the EWCE::generate_celltype_data
#'
#' @examples
#' data(ctd_obj)
"ctd_obj"
