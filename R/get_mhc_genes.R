#' Get vector of hg38 genes overlapping MHC region
#'
#' @return A dataframe of genes, positions and metadata overlapping the MHC region of chr6.
#' @export
#'
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom httr GET stop_for_status content content_type
#' @importFrom dplyr %>% relocate bind_rows as_tibble
#' @importFrom purrr map
#' @importFrom stringi stri_remove_empty
#'
#' @param build A genome build to query the MHC genes from (hg19 or hg38)
#' @param method The query method to use (biomart or rest_api)
#'
#'
#' @examples get_mhc_genes(build = 'hg38', method = 'rest_api')
#' get_mhc_genes(build = 'hg19', method = 'biomart')
get_mhc_genes <- function(build = 'hg38', method = 'rest_api') {

  stopifnot(build == 'hg19' || build == 'hg38')

  if (build == 'hg19') {

    biomart_server <- "https://grch37.ensembl.org"
    rest_server <- "https://grch37.rest.ensembl.org"
    coord <- "chr6:28,477,797-33,448,354"
    start <- "28477797"
    end <- "33448354"


  } else {

    biomart_server <- "https://www.ensembl.org"
    rest_server <- "https://rest.ensembl.org"
    coord <- "chr6:28,510,120-33,480,577"
    start <- "28510120"
    end <- "33480577"

  }

  if (method == 'biomart') {

  options(timeout = 30000)
  message('\n\nConnecting to Ensembl server via BiomaRt ...')
  mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl", host = biomart_server)

  message('\nCalling MHC region genes ...\n')
  mhc_genes_df <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position",
                                             "gene_biotype", "source", "version", "ensembl_gene_id"),
                     filters = c("chromosome_name", "start", "end"),
                     values = list(chromosome = "6", start = start, end = end),
                     mart = mart)

  }

  if (method == 'rest_api') {


    message('\n\nConnecting to Ensembl server via REST API...')
    ext <- paste0("/overlap/region/human/", coord, "?feature=gene")
    r <- httr::GET(paste(rest_server, ext, sep = ""), httr::content_type("text/csv"))
    httr::stop_for_status(r)

    mhc_genes_df <- print(httr::content(r)) %>%
      purrr::map(unlist) %>%
      purrr::map(t) %>%
      purrr::map(dplyr::as_tibble) %>%
      dplyr::bind_rows() %>%
      dplyr::relocate(.data$external_name, .data$seq_region_name, .data$start,
                      .data$end, .data$assembly_name, .data$biotype, .data$description)



  }


  mhc_genes_uniq <- stringi::stri_remove_empty(unique(mhc_genes_df$external_name), na_empty = FALSE)
  message('MHC genes detected: ', length(mhc_genes_df$external_name))
  message('MHC genes with unique HGNC ID detected: ', length(mhc_genes_uniq))
  return(mhc_genes_df)

}


