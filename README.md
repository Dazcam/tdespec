
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tdespec

<!-- badges: start -->
<!-- badges: end -->

The goal of tdespec is to generate cell specificity scores for scRNAseq
dataset as described in the [Skene et
al. (2018)](https://www.nature.com/articles/s41588-018-0129-5). This
utilises the functionality of the Bioconductor package
[EWCE](https://github.com/NathanSkene/EWCE/) to generate a CTD object
containing the specificity scores for each cell type in a snRNAseq
dataset and generates files for gene set enrichment analyses packages
[MAGMA](https://ctg.cncr.nl/software/magma) and [LD Score
regression](https://github.com/bulik/ldsc).

## Installation

You can install the development version of tdespec like so:

``` r
if (!require(devtools)) { install.packages('devtools') } 
library(devtools)
install_github("Dazcam/tdespec")
```

## Basic workflow

The basic `tdespec` workflow is a 4-stage process:

1.  Create a vector of genes overlapping the MHC region of the human
    genome.
2.  Obtain a a tibble of a gene reference file with MHC overlapping
    genes removed.
3.  Create a ctd object from and snRNAseq gene expression matrix.
4.  Create input files for gene enrichment tests (MAGMA and LDSR) and
    save to file.

``` r
library(tdespec)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
data("seurat_small")
mhc_genes <- get_mhc_genes(build = 'hg19', method = 'rest_api')
#> 
#> 
#> Connecting to Ensembl server via REST API...
#> MHC genes detected: 382
#> MHC genes with unique HGNC ID detected: 374
gene_coord <- get_gene_coordinates(mhc_genes %>% pull(external_name))
#> 
#> 
#> Dowloading MAGMA hg19 reference file ...
ctd_path <- create_ctd(seurat_small@assays$RNA@data, 
                       meta_lev_1 = seurat_small$cluster_level_1,
                       meta_lev_2 = seurat_small$cluster_level_2,
                       outdir = 'test/',
                       mhc_genes = mhc_genes,
                       group_name = 'brain_study')
#> 
#> 
#> Creating CTD object ...
#> MHC genes removed: 0
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Running CPM ...
#> 
#> Dropping uninformative genes cpm norm ...
#> Check 21310Check 1452
#> 1 core(s) assigned as workers (7 reserved).
#> Converting to sparse matrix.
#> Checking for non-expressed genes.
#> -2918 / 18392 non-expressed genes dropped
#> Checking for cells with no expressed genes.
#> DGE:: Limma...
#> 8,495 / 18,392 genes dropped @ DGE adj_pval_thresh < 1e-05
#> Time difference of 4.635189 secs
#> 
#> Gene counts:
#> 
#> RAW:21310
#> RAW_NO_MHC:21310
#> CPM_DROP_GENES:9897
#> 
#> Creating CTD object ...
#> 1 core(s) assigned as workers (7 reserved).
#> Converting to sparse matrix.
#> + Calculating normalized mean expression.
#> Converting to sparse matrix.
#> Converting to sparse matrix.
#> + Calculating normalized specificity.
#> Converting to sparse matrix.
#> Converting to sparse matrix.
#> Converting to sparse matrix.
#> Converting to sparse matrix.
#> Loading required namespace: ggdendro
#> + Saving results ==>  test//ctd_brain_study.rda
create_enrich_test_files(ctd_path = ctd_path, 
                         study_id = 'brain_study',
                         gene_coordinates = gene_coord,
                         outdir = 'test/')
#> Creating MAGMA input files for cluster level 1 ...
#> Joining with `by = join_by(hgnc)`
#> Creating SLDSR input files for cluster level: 1...
#> Joining with `by = join_by(hgnc)`
#> Output files have been saved to: test/.
#> Done.
#> Creating MAGMA input files for cluster level 2 ...
#> Joining with `by = join_by(hgnc)`
#> Creating SLDSR input files for cluster level: 2...
#> Joining with `by = join_by(hgnc)`
#> Output files have been saved to: test/.
#> Done.
```
