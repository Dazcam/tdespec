#' Create files for SLDSR and MAGMA gene set enrichment tests from ctd object (Top 10%)
#'
#' @return Writes a set of tsv MAGMA / SLDSR files to directory of interest
#' @export
#'
#' @importFrom dplyr inner_join filter select mutate group_by group_walk
#' @importFrom tibble as_tibble
#' @importFrom tidyselect all_of
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom readr write_tsv
#' @importFrom tools file_path_sans_ext
#' @importFrom utils globalVariables
#'
#' @param ctd_path A string of file path for ctd object
#' @param study_id A unique identifier for the study
#' @param gene_coordinates A gene coordinate reference file
#' @param levels A vector of levels
#' @param outdir A directory to save the enrichment test files to.
#' @param sldsr_ext A vector of the upstream and downstream gene / peak boundary
#' extension window for SLDSR
#'
#'
#' @examples get_gene_coordinates(c('SPI1', 'MEF2C', 'GAD1', 'TREM2', 'NEUROD1'))
create_enrich_test_files <- function(

  ctd_path = NULL,
  study_id = NULL,
  gene_coordinates = NULL,
  levels = c(1, 2),
  outdir = NULL,
  sldsr_ext = c(100, 100)

  ) {

  # Housekeeping
  options(scipen = 999) # Turn of scientifc notation
  utils::globalVariables(".") # Silence R CMD check https://stackoverflow.com/questions/66816638/


  dir.create(paste0(outdir, 'MAGMA/'),  recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, 'LDSR/'), recursive = TRUE, showWarnings = FALSE)

  # Strip basename from ctd object as study id
  if(is.null(study_id)) { tools::file_path_sans_ext(study_id) }

  for (level in levels) {

    message('Creating MAGMA input files for cluster level: ', level, '...')

    # Load ctd object need function as file is saved as .rda
    loadRData <- function(fileName){ # https://stackoverflow.com/questions/5577221
      #loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    ctd_obj <- loadRData(ctd_path)

    cell_types <- colnames(ctd_obj[[level]]$specificity_quantiles)

    MAGMA <- dplyr::as_tibble(as.matrix(ctd_obj[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
      dplyr::inner_join(gene_coordinates) %>%
      tidyr::pivot_longer(tidyselect::all_of(cell_types), names_to = 'cell_type', values_to = 'quantile') %>%
      dplyr::filter(.data$quantile == 10) %>%
      dplyr::select(.data$cell_type, .data$entrez) %>%
      with(., split(.data$entrez, cell_type))

    for(i in names(MAGMA)) {

      cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
          file = paste0(outdir, 'MAGMA/', study_id, '_', level, '.txt')
          , sep = '', append = TRUE)

    }

    upstream_ext <- sldsr_ext[1] * 1000
    downstream_ext <- sldsr_ext[1] * 1000
    window <- paste0(sldsr_ext[1], 'UP_', sldsr_ext[2], 'DOWN')

    message('Creating SLDSR input files for cluster level: ', level, '...')
    ldsr <- tibble::as_tibble(as.matrix(ctd_obj[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
      dplyr::inner_join(gene_coordinates) %>%
      tidyr::pivot_longer(tidyselect::all_of(cell_types), names_to = 'cell_type', values_to = 'quantile') %>%
      dplyr::filter(.data$quantile == 10) %>%
      dplyr::mutate(start = ifelse(.data$start - upstream_ext < 0, 0,
                                   .data$start - upstream_ext),
                    end = .data$end + downstream_ext) %>%
      dplyr::select(.data$chr, .data$start, .data$end, .data$entrez, .data$cell_type) %>%
      dplyr::group_by(.data$cell_type) %>%
      dplyr::group_walk(~ readr::write_tsv(.x[,1:4], paste0(outdir, 'LDSR/',
                                                            .y$cell_type, '.',
                                                            window, '.bed'),
                                           col_names = FALSE))

    message('Output files have been saved to: ', outdir, '.')
    message('Done.')

    }

 }




