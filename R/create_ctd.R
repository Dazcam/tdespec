#' Create gene expression specificity scores for scRNAseq data
#'
#' @return A CTD object containing
#'
#' @export
#'
#' @importFrom EWCE generate_celltype_data drop_uninformative_genes
#' @importFrom edgeR cpm
#'
#' @param gex_matrix A raw count gene expression matrix - needs to be genes x cell.
#' @param meta_lev_1 A vector of cell-type labels for level 1 cell types.
#' @param meta_lev_2  A vector of cell-type labels for level 2 cell types.
#' @param outdir The output directory to store the ctd object.
#' @param mhc_genes A vector of HGNC ID genes overlapping the MHC to be removed.
#' @param group_name An unique identifier (string) for the dataset.
#' @param threads The number of cores to assign to process. Default is NULL,
#' which will auto detect total cores and estimate number to use.
#'
#' @returns A path to the location where where the ctd object was stored.
#'
#' @examples
#' data('seurat_small', 'mhc_genes')
#' create_ctd(gex_matrix = seurat_small@assays$RNA@data,
#'                     meta_lev_1 = seurat_small$cluster_level_1,
#'                     meta_lev_2 = seurat_small$cluster_level_2,
#'                     outdir = 'test',
#'                     mhc_genes = mhc_genes,
#'                     group_name = 'brain_study',
#'                     threads = 1)

create_ctd <- function(gex_matrix = NULL, meta_lev_1 = NULL, meta_lev_2 = NULL,
                       outdir = NULL, mhc_genes = NULL, group_name = NULL,
                       threads = NULL) {

  message('\n\nCreating CTD object ... \n')

  # Remove MHC genes from gex_matrix
  gex_matrix_no_mhc <- gex_matrix[!(rownames(gex_matrix) %in% mhc_genes), ]
  message('MHC genes removed: ', dim(gex_matrix)[1] - dim(gex_matrix_no_mhc)[1])

  # Can also be done directly to seurat object!!
  # gex_matrix_no_mhc <- seurat_obj[!(rownames(seurat_obj) %in% mhc_genes), ]

  # Create annotations
  annotations <- as.data.frame(cbind(as.vector(colnames(gex_matrix)),
                                     as.vector(meta_lev_1),
                                     as.vector(meta_lev_2)))
  colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
  rownames(annotations) <- NULL
  annotLevels <- list(level1class = annotations$level1class,
                      level2class = annotations$level2class)

  # Normalize - this is optional, was not used in the original EWCE publication
  # cat('\nRunning SCT ... ', '\n\n')
  # COUNTS_SCT <- EWCE::sct_normalize(RAW_COUNTS_NO_MHC)

  message('\nRunning CPM ...')
  gex_matrix_cpm <- edgeR::cpm(gex_matrix_no_mhc)

  # cat('\nDropping uninformative genes raw ... ', '\n\n')
  # DROP_GENES_RAW <- EWCE::drop_uninformative_genes(
  #   exp = RAW_COUNTS_NO_MHC,
  #   input_species = "human",
  #   output_species = "human",
  #   level2annot = annotLevels$level2class)
  #
  # cat('\nDropping uninformative genes sct norm ... ', '\n\n')
  # DROP_GENES_SCT <- EWCE::drop_uninformative_genes(
  #   exp = COUNTS_SCT,
  #   input_species = "human",
  #   output_species = "human",
  #   level2annot = annotLevels$level2class)

  message('\nDropping uninformative genes cpm norm ...')
  gex_matrix_drop_genes <- EWCE::drop_uninformative_genes(
    exp = gex_matrix_cpm,
    input_species = "human",
    output_species = "human",
    level2annot = annotLevels$level2class,
    no_cores = threads)

  message('\nGene counts:',
      '\n\nRAW:', dim(gex_matrix)[1],
      '\nRAW_NO_MHC:', dim(gex_matrix_cpm)[1],
      '\nCPM_DROP_GENES:', dim(gex_matrix_drop_genes)[1])

  # cat('\nGene counts:',
  #     '\n\nRAW:', dim(RAW_COUNTS)[1],
  #     '\nRAW_NO_MHC:', dim(RAW_COUNTS_NO_MHC)[1],
  #     '\nNO_NORM_DROP_GENES:', dim(DROP_GENES_RAW)[1],
  #     '\nSCT_DROP_GENES:', dim(DROP_GENES_SCT)[1],
  #     '\nCPM_DROP_GENES:', dim(DROP_GENES_CPM)[1])

  # Create object - saves ctd obj to folder
  message('\nCreating CTD object ... \n\n')
  dir.create(outdir, showWarnings = FALSE)
  ctd <- EWCE::generate_celltype_data(exp = gex_matrix_drop_genes,
                                      annotLevels = annotLevels,
                                      groupName = group_name,
                                      savePath = outdir,
                                      numberOfBins = 10,
                                      no_cores = threads)

  return(ctd)


}



