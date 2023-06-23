#' Create gene expression specificity scores for scRNAseq data
#'
#' @return A CTD object containing
#'
#' @export
#'
#' @importFrom EWCE generate_celltype_data drop_uninformative_genes
#' @importFrom edgeR cpm
#'
#' @param matrix Raw count gene matrix - needs to be genes x cell and annotation data
#' @param meta_lev_1 A vector of cell IDs for level 1 cell types
#' @param meta_lev_2  A vector of cell IDs for level 1 cell types
#' @param outdir The output directory to store the ctd object
#' @param mhc_genes A set of MHC genes to be removed
#' @param group_name An string identifier for the dataset
#'
#' @returns directory where ctd object was stored
#'
#' @examples

create_ctd <- function(matrix = NULL, meta_lev_1 = NULL, meta_lev_2 = NULL,
                       outdir = NULL, mhc_genes = NULL, group_name = NULL) {

  message('\n\nCreating CTD object ... \n')

  # Remove MHC genes from matrix
  matrix_no_mhc <- matrix[!(rownames(matrix) %in% mhc_genes), ]
  message('MHC genes removed: ', dim(matrix)[1] - dim(matrix_no_mhc)[1])

  # Create annotations
  annotations <- as.data.frame(cbind(as.vector(colnames(matrix)),
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
  matrix_cpm <- edgeR::cpm(matrix_no_mhc)

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
  matrix_drop_genes <- EWCE::drop_uninformative_genes(
    exp = matrix_cpm,
    input_species = "human",
    output_species = "human",
    level2annot = annotLevels$level2class)

  message('\nGene counts:',
      '\n\nRAW:', dim(matrix)[1],
      '\nRAW_NO_MHC:', dim(matrix_cpm)[1],
      '\nCPM_DROP_GENES:', dim(matrix_drop_genes)[1])

  # cat('\nGene counts:',
  #     '\n\nRAW:', dim(RAW_COUNTS)[1],
  #     '\nRAW_NO_MHC:', dim(RAW_COUNTS_NO_MHC)[1],
  #     '\nNO_NORM_DROP_GENES:', dim(DROP_GENES_RAW)[1],
  #     '\nSCT_DROP_GENES:', dim(DROP_GENES_SCT)[1],
  #     '\nCPM_DROP_GENES:', dim(DROP_GENES_CPM)[1])

  # Create object - saves ctd obj to folder
  message('\nCreating CTD object ... \n\n')
  dir.create(outdir, showWarnings = FALSE)
  ctd <- EWCE::generate_celltype_data(exp = matrix_drop_genes,
                                      annotLevels = annotLevels,
                                      groupName = group_name,
                                      savePath = outdir,
                                      numberOfBins = 10)

  return(ctd)


}



