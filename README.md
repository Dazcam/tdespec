
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
2.  Obtain a a tibble of a gene reference file with MHC overalpping
    genes removed.
3.  Create a ctd object from and snRNAseq gene expression matrix.
4.  Create input files for gene enrichment tests (MAGMA and LDSR)

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
mhc_genes <- get_mhc_genes(build = 'hg19', method = 'rest_api')
#> 
#> 
#> Connecting to Ensembl server via REST API...
#> [[1]]
#> [[1]]$external_name
#> [1] "CSNK2B-LY6G5B-1181"
#> 
#> [[1]]$canonical_transcript
#> [1] "ENST00000375880.2"
#> 
#> [[1]]$end
#> [1] 31641323
#> 
#> [[1]]$id
#> [1] "ENSG00000263020"
#> 
#> [[1]]$start
#> [1] 31633879
#> 
#> [[1]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[1]]$version
#> [1] 1
#> 
#> [[1]]$strand
#> [1] 1
#> 
#> [[1]]$gene_id
#> [1] "ENSG00000263020"
#> 
#> [[1]]$source
#> [1] "havana"
#> 
#> [[1]]$seq_region_name
#> [1] "6"
#> 
#> [[1]]$biotype
#> [1] "protein_coding"
#> 
#> [[1]]$feature_type
#> [1] "gene"
#> 
#> [[1]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[1]]$description
#> [1] "Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:H0Y6P8]"
#> 
#> 
#> [[2]]
#> [[2]]$start
#> [1] 31637944
#> 
#> [[2]]$external_name
#> [1] "LY6G5B"
#> 
#> [[2]]$canonical_transcript
#> [1] "ENST00000375864.4"
#> 
#> [[2]]$id
#> [1] "ENSG00000240053"
#> 
#> [[2]]$end
#> [1] 31641553
#> 
#> [[2]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[2]]$strand
#> [1] 1
#> 
#> [[2]]$version
#> [1] 8
#> 
#> [[2]]$gene_id
#> [1] "ENSG00000240053"
#> 
#> [[2]]$source
#> [1] "ensembl_havana"
#> 
#> [[2]]$seq_region_name
#> [1] "6"
#> 
#> [[2]]$description
#> [1] "lymphocyte antigen 6 complex, locus G5B [Source:HGNC Symbol;Acc:13931]"
#> 
#> [[2]]$biotype
#> [1] "protein_coding"
#> 
#> [[2]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[2]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[3]]
#> [[3]]$seq_region_name
#> [1] "6"
#> 
#> [[3]]$source
#> [1] "ensembl_havana"
#> 
#> [[3]]$description
#> [1] "hydroxysteroid (17-beta) dehydrogenase 8 [Source:HGNC Symbol;Acc:3554]"
#> 
#> [[3]]$feature_type
#> [1] "gene"
#> 
#> [[3]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[3]]$biotype
#> [1] "protein_coding"
#> 
#> [[3]]$start
#> [1] 33172419
#> 
#> [[3]]$id
#> [1] "ENSG00000204228"
#> 
#> [[3]]$end
#> [1] 33174608
#> 
#> [[3]]$canonical_transcript
#> [1] "ENST00000374662.3"
#> 
#> [[3]]$external_name
#> [1] "HSD17B8"
#> 
#> [[3]]$gene_id
#> [1] "ENSG00000204228"
#> 
#> [[3]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[3]]$version
#> [1] 3
#> 
#> [[3]]$strand
#> [1] 1
#> 
#> 
#> [[4]]
#> [[4]]$gene_id
#> [1] "ENSG00000204428"
#> 
#> [[4]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[4]]$version
#> [1] 8
#> 
#> [[4]]$strand
#> [1] -1
#> 
#> [[4]]$start
#> [1] 31644461
#> 
#> [[4]]$end
#> [1] 31651817
#> 
#> [[4]]$id
#> [1] "ENSG00000204428"
#> 
#> [[4]]$external_name
#> [1] "LY6G5C"
#> 
#> [[4]]$canonical_transcript
#> [1] "ENST00000383237.4"
#> 
#> [[4]]$description
#> [1] "lymphocyte antigen 6 complex, locus G5C [Source:HGNC Symbol;Acc:13932]"
#> 
#> [[4]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[4]]$feature_type
#> [1] "gene"
#> 
#> [[4]]$biotype
#> [1] "protein_coding"
#> 
#> [[4]]$seq_region_name
#> [1] "6"
#> 
#> [[4]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[5]]
#> [[5]]$seq_region_name
#> [1] "6"
#> 
#> [[5]]$source
#> [1] "ensembl_havana"
#> 
#> [[5]]$description
#> [1] "ring finger protein 1 [Source:HGNC Symbol;Acc:10018]"
#> 
#> [[5]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[5]]$feature_type
#> [1] "gene"
#> 
#> [[5]]$biotype
#> [1] "protein_coding"
#> 
#> [[5]]$start
#> [1] 33176272
#> 
#> [[5]]$id
#> [1] "ENSG00000204227"
#> 
#> [[5]]$end
#> [1] 33180499
#> 
#> [[5]]$canonical_transcript
#> [1] "ENST00000374656.4"
#> 
#> [[5]]$external_name
#> [1] "RING1"
#> 
#> [[5]]$gene_id
#> [1] "ENSG00000204227"
#> 
#> [[5]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[5]]$version
#> [1] 4
#> 
#> [[5]]$strand
#> [1] 1
#> 
#> 
#> [[6]]
#> [[6]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[6]]$feature_type
#> [1] "gene"
#> 
#> [[6]]$biotype
#> [1] "protein_coding"
#> 
#> [[6]]$description
#> [1] "abhydrolase domain containing 16A [Source:HGNC Symbol;Acc:13921]"
#> 
#> [[6]]$seq_region_name
#> [1] "6"
#> 
#> [[6]]$source
#> [1] "ensembl_havana"
#> 
#> [[6]]$gene_id
#> [1] "ENSG00000204427"
#> 
#> [[6]]$version
#> [1] 7
#> 
#> [[6]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[6]]$strand
#> [1] -1
#> 
#> [[6]]$id
#> [1] "ENSG00000204427"
#> 
#> [[6]]$end
#> [1] 31671221
#> 
#> [[6]]$external_name
#> [1] "ABHD16A"
#> 
#> [[6]]$canonical_transcript
#> [1] "ENST00000395952.3"
#> 
#> [[6]]$start
#> [1] 31654726
#> 
#> 
#> [[7]]
#> [[7]]$description
#> [1] "vacuolar protein sorting 52 homolog (S. cerevisiae) [Source:HGNC Symbol;Acc:10518]"
#> 
#> [[7]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[7]]$feature_type
#> [1] "gene"
#> 
#> [[7]]$biotype
#> [1] "protein_coding"
#> 
#> [[7]]$seq_region_name
#> [1] "6"
#> 
#> [[7]]$source
#> [1] "ensembl_havana"
#> 
#> [[7]]$gene_id
#> [1] "ENSG00000223501"
#> 
#> [[7]]$version
#> [1] 4
#> 
#> [[7]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[7]]$strand
#> [1] -1
#> 
#> [[7]]$start
#> [1] 33218049
#> 
#> [[7]]$end
#> [1] 33239824
#> 
#> [[7]]$id
#> [1] "ENSG00000223501"
#> 
#> [[7]]$external_name
#> [1] "VPS52"
#> 
#> [[7]]$canonical_transcript
#> [1] "ENST00000445902.2"
#> 
#> 
#> [[8]]
#> [[8]]$start
#> [1] 28555155
#> 
#> [[8]]$end
#> [1] 28559524
#> 
#> [[8]]$id
#> [1] "ENSG00000246350"
#> 
#> [[8]]$canonical_transcript
#> [1] "ENST00000499525.1"
#> 
#> [[8]]$external_name
#> [1] "RP5-1186N24.3"
#> 
#> [[8]]$gene_id
#> [1] "ENSG00000246350"
#> 
#> [[8]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[8]]$version
#> [1] 1
#> 
#> [[8]]$strand
#> [1] 1
#> 
#> [[8]]$seq_region_name
#> [1] "6"
#> 
#> [[8]]$source
#> [1] "havana"
#> 
#> [[8]]$description
#> NULL
#> 
#> [[8]]$feature_type
#> [1] "gene"
#> 
#> [[8]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[8]]$biotype
#> [1] "processed_transcript"
#> 
#> 
#> [[9]]
#> [[9]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[9]]$version
#> [1] 1
#> 
#> [[9]]$strand
#> [1] -1
#> 
#> [[9]]$gene_id
#> [1] "ENSG00000271440"
#> 
#> [[9]]$external_name
#> [1] "RP11-373N24.2"
#> 
#> [[9]]$canonical_transcript
#> [1] "ENST00000605542.1"
#> 
#> [[9]]$end
#> [1] 28601371
#> 
#> [[9]]$id
#> [1] "ENSG00000271440"
#> 
#> [[9]]$start
#> [1] 28601158
#> 
#> [[9]]$biotype
#> [1] "pseudogene"
#> 
#> [[9]]$feature_type
#> [1] "gene"
#> 
#> [[9]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[9]]$description
#> NULL
#> 
#> [[9]]$source
#> [1] "havana"
#> 
#> [[9]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[10]]
#> [[10]]$start
#> [1] 28615974
#> 
#> [[10]]$external_name
#> [1] "AL121932.1"
#> 
#> [[10]]$canonical_transcript
#> [1] "ENST00000606479.1"
#> 
#> [[10]]$end
#> [1] 28616042
#> 
#> [[10]]$id
#> [1] "ENSG00000272278"
#> 
#> [[10]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[10]]$strand
#> [1] -1
#> 
#> [[10]]$version
#> [1] 1
#> 
#> [[10]]$gene_id
#> [1] "ENSG00000272278"
#> 
#> [[10]]$source
#> [1] "ensembl"
#> 
#> [[10]]$seq_region_name
#> [1] "6"
#> 
#> [[10]]$description
#> NULL
#> 
#> [[10]]$biotype
#> [1] "miRNA"
#> 
#> [[10]]$feature_type
#> [1] "gene"
#> 
#> [[10]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> 
#> [[11]]
#> [[11]]$seq_region_name
#> [1] "6"
#> 
#> [[11]]$source
#> [1] "havana"
#> 
#> [[11]]$description
#> [1] "long intergenic non-protein coding RNA 533 [Source:HGNC Symbol;Acc:18690]"
#> 
#> [[11]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[11]]$feature_type
#> [1] "gene"
#> 
#> [[11]]$biotype
#> [1] "lincRNA"
#> 
#> [[11]]$start
#> [1] 28616063
#> 
#> [[11]]$end
#> [1] 28616843
#> 
#> [[11]]$id
#> [1] "ENSG00000235570"
#> 
#> [[11]]$external_name
#> [1] "LINC00533"
#> 
#> [[11]]$canonical_transcript
#> [1] "ENST00000456103.1"
#> 
#> [[11]]$gene_id
#> [1] "ENSG00000235570"
#> 
#> [[11]]$strand
#> [1] -1
#> 
#> [[11]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[11]]$version
#> [1] 1
#> 
#> 
#> [[12]]
#> [[12]]$seq_region_name
#> [1] "6"
#> 
#> [[12]]$source
#> [1] "havana"
#> 
#> [[12]]$feature_type
#> [1] "gene"
#> 
#> [[12]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[12]]$biotype
#> [1] "pseudogene"
#> 
#> [[12]]$description
#> [1] "ribosomal protein SA pseudogene 2 [Source:HGNC Symbol;Acc:18771]"
#> 
#> [[12]]$end
#> [1] 28700681
#> 
#> [[12]]$id
#> [1] "ENSG00000237425"
#> 
#> [[12]]$external_name
#> [1] "RPSAP2"
#> 
#> [[12]]$canonical_transcript
#> [1] "ENST00000416785.1"
#> 
#> [[12]]$start
#> [1] 28699794
#> 
#> [[12]]$gene_id
#> [1] "ENSG00000237425"
#> 
#> [[12]]$strand
#> [1] 1
#> 
#> [[12]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[12]]$version
#> [1] 1
#> 
#> 
#> [[13]]
#> [[13]]$description
#> NULL
#> 
#> [[13]]$feature_type
#> [1] "gene"
#> 
#> [[13]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[13]]$biotype
#> [1] "miRNA"
#> 
#> [[13]]$seq_region_name
#> [1] "6"
#> 
#> [[13]]$source
#> [1] "ensembl"
#> 
#> [[13]]$gene_id
#> [1] "ENSG00000221191"
#> 
#> [[13]]$strand
#> [1] 1
#> 
#> [[13]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[13]]$version
#> [1] 1
#> 
#> [[13]]$start
#> [1] 28743693
#> 
#> [[13]]$end
#> [1] 28743800
#> 
#> [[13]]$id
#> [1] "ENSG00000221191"
#> 
#> [[13]]$external_name
#> [1] "AL662890.1"
#> 
#> [[13]]$canonical_transcript
#> [1] "ENST00000408264.1"
#> 
#> 
#> [[14]]
#> [[14]]$source
#> [1] "havana"
#> 
#> [[14]]$seq_region_name
#> [1] "6"
#> 
#> [[14]]$description
#> NULL
#> 
#> [[14]]$biotype
#> [1] "pseudogene"
#> 
#> [[14]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[14]]$feature_type
#> [1] "gene"
#> 
#> [[14]]$start
#> [1] 28751410
#> 
#> [[14]]$external_name
#> [1] "NOL5BP"
#> 
#> [[14]]$canonical_transcript
#> [1] "ENST00000440030.1"
#> 
#> [[14]]$id
#> [1] "ENSG00000235559"
#> 
#> [[14]]$end
#> [1] 28751781
#> 
#> [[14]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[14]]$strand
#> [1] -1
#> 
#> [[14]]$version
#> [1] 1
#> 
#> [[14]]$gene_id
#> [1] "ENSG00000235559"
#> 
#> 
#> [[15]]
#> [[15]]$description
#> NULL
#> 
#> [[15]]$feature_type
#> [1] "gene"
#> 
#> [[15]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[15]]$biotype
#> [1] "miRNA"
#> 
#> [[15]]$seq_region_name
#> [1] "6"
#> 
#> [[15]]$source
#> [1] "ensembl"
#> 
#> [[15]]$gene_id
#> [1] "ENSG00000265764"
#> 
#> [[15]]$version
#> [1] 1
#> 
#> [[15]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[15]]$strand
#> [1] -1
#> 
#> [[15]]$start
#> [1] 28763749
#> 
#> [[15]]$id
#> [1] "ENSG00000265764"
#> 
#> [[15]]$end
#> [1] 28763812
#> 
#> [[15]]$external_name
#> [1] "AL662890.3"
#> 
#> [[15]]$canonical_transcript
#> [1] "ENST00000577570.1"
#> 
#> 
#> [[16]]
#> [[16]]$start
#> [1] 28805646
#> 
#> [[16]]$id
#> [1] "ENSG00000225173"
#> 
#> [[16]]$end
#> [1] 28806783
#> 
#> [[16]]$external_name
#> [1] "XXbac-BPG308K3.5"
#> 
#> [[16]]$canonical_transcript
#> [1] "ENST00000457253.1"
#> 
#> [[16]]$gene_id
#> [1] "ENSG00000225173"
#> 
#> [[16]]$strand
#> [1] -1
#> 
#> [[16]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[16]]$version
#> [1] 1
#> 
#> [[16]]$seq_region_name
#> [1] "6"
#> 
#> [[16]]$source
#> [1] "havana"
#> 
#> [[16]]$description
#> NULL
#> 
#> [[16]]$feature_type
#> [1] "gene"
#> 
#> [[16]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[16]]$biotype
#> [1] "lincRNA"
#> 
#> 
#> [[17]]
#> [[17]]$description
#> NULL
#> 
#> [[17]]$biotype
#> [1] "lincRNA"
#> 
#> [[17]]$feature_type
#> [1] "gene"
#> 
#> [[17]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[17]]$source
#> [1] "havana"
#> 
#> [[17]]$seq_region_name
#> [1] "6"
#> 
#> [[17]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[17]]$version
#> [1] 2
#> 
#> [[17]]$strand
#> [1] -1
#> 
#> [[17]]$gene_id
#> [1] "ENSG00000225595"
#> 
#> [[17]]$start
#> [1] 28827402
#> 
#> [[17]]$canonical_transcript
#> [1] "ENST00000440244.1"
#> 
#> [[17]]$external_name
#> [1] "XXbac-BPG308K3.6"
#> 
#> [[17]]$end
#> [1] 28832407
#> 
#> [[17]]$id
#> [1] "ENSG00000225595"
#> 
#> 
#> [[18]]
#> [[18]]$start
#> [1] 28829193
#> 
#> [[18]]$canonical_transcript
#> [1] "ENST00000396818.2"
#> 
#> [[18]]$external_name
#> [1] "RPL13P"
#> 
#> [[18]]$id
#> [1] "ENSG00000213916"
#> 
#> [[18]]$end
#> [1] 28829824
#> 
#> [[18]]$version
#> [1] 2
#> 
#> [[18]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[18]]$strand
#> [1] 1
#> 
#> [[18]]$gene_id
#> [1] "ENSG00000213916"
#> 
#> [[18]]$source
#> [1] "havana"
#> 
#> [[18]]$seq_region_name
#> [1] "6"
#> 
#> [[18]]$description
#> [1] "ribosomal protein L13 pseudogene [Source:HGNC Symbol;Acc:13978]"
#> 
#> [[18]]$biotype
#> [1] "pseudogene"
#> 
#> [[18]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[18]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[19]]
#> [[19]]$start
#> [1] 28856609
#> 
#> [[19]]$end
#> [1] 28857683
#> 
#> [[19]]$id
#> [1] "ENSG00000233366"
#> 
#> [[19]]$canonical_transcript
#> [1] "ENST00000419323.1"
#> 
#> [[19]]$external_name
#> [1] "ZNF90P2"
#> 
#> [[19]]$gene_id
#> [1] "ENSG00000233366"
#> 
#> [[19]]$strand
#> [1] -1
#> 
#> [[19]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[19]]$version
#> [1] 1
#> 
#> [[19]]$seq_region_name
#> [1] "6"
#> 
#> [[19]]$source
#> [1] "havana"
#> 
#> [[19]]$description
#> [1] "zinc finger protein 90 pseudogene 2 [Source:HGNC Symbol;Acc:21687]"
#> 
#> [[19]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[19]]$feature_type
#> [1] "gene"
#> 
#> [[19]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[20]]
#> [[20]]$external_name
#> [1] "HCG14"
#> 
#> [[20]]$canonical_transcript
#> [1] "ENST00000444986.1"
#> 
#> [[20]]$end
#> [1] 28865099
#> 
#> [[20]]$id
#> [1] "ENSG00000224157"
#> 
#> [[20]]$start
#> [1] 28864307
#> 
#> [[20]]$version
#> [1] 1
#> 
#> [[20]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[20]]$strand
#> [1] 1
#> 
#> [[20]]$gene_id
#> [1] "ENSG00000224157"
#> 
#> [[20]]$source
#> [1] "havana"
#> 
#> [[20]]$seq_region_name
#> [1] "6"
#> 
#> [[20]]$biotype
#> [1] "antisense"
#> 
#> [[20]]$feature_type
#> [1] "gene"
#> 
#> [[20]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[20]]$description
#> [1] "HLA complex group 14 (non-protein coding) [Source:HGNC Symbol;Acc:18323]"
#> 
#> 
#> [[21]]
#> [[21]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[21]]$strand
#> [1] 1
#> 
#> [[21]]$version
#> [1] 1
#> 
#> [[21]]$gene_id
#> [1] "ENSG00000212240"
#> 
#> [[21]]$start
#> [1] 28883422
#> 
#> [[21]]$canonical_transcript
#> [1] "ENST00000390938.1"
#> 
#> [[21]]$external_name
#> [1] "RNU6-930P"
#> 
#> [[21]]$id
#> [1] "ENSG00000212240"
#> 
#> [[21]]$end
#> [1] 28883521
#> 
#> [[21]]$description
#> [1] "RNA, U6 small nuclear 930, pseudogene [Source:HGNC Symbol;Acc:47893]"
#> 
#> [[21]]$biotype
#> [1] "snRNA"
#> 
#> [[21]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[21]]$feature_type
#> [1] "gene"
#> 
#> [[21]]$source
#> [1] "ensembl"
#> 
#> [[21]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[22]]
#> [[22]]$start
#> [1] 28936907
#> 
#> [[22]]$end
#> [1] 28938237
#> 
#> [[22]]$id
#> [1] "ENSG00000228666"
#> 
#> [[22]]$canonical_transcript
#> [1] "ENST00000452408.1"
#> 
#> [[22]]$external_name
#> [1] "KRT18P1"
#> 
#> [[22]]$gene_id
#> [1] "ENSG00000228666"
#> 
#> [[22]]$version
#> [1] 1
#> 
#> [[22]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[22]]$strand
#> [1] -1
#> 
#> [[22]]$seq_region_name
#> [1] "6"
#> 
#> [[22]]$source
#> [1] "havana"
#> 
#> [[22]]$description
#> [1] "keratin 18 pseudogene 1 [Source:HGNC Symbol;Acc:6434]"
#> 
#> [[22]]$feature_type
#> [1] "gene"
#> 
#> [[22]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[22]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[23]]
#> [[23]]$source
#> [1] "ensembl"
#> 
#> [[23]]$seq_region_name
#> [1] "6"
#> 
#> [[23]]$description
#> [1] "RNA, 7SL, cytoplasmic 471, pseudogene [Source:HGNC Symbol;Acc:46487]"
#> 
#> [[23]]$biotype
#> [1] "misc_RNA"
#> 
#> [[23]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[23]]$feature_type
#> [1] "gene"
#> 
#> [[23]]$start
#> [1] 28945252
#> 
#> [[23]]$canonical_transcript
#> [1] "ENST00000580476.1"
#> 
#> [[23]]$external_name
#> [1] "RN7SL471P"
#> 
#> [[23]]$end
#> [1] 28945547
#> 
#> [[23]]$id
#> [1] "ENSG00000263426"
#> 
#> [[23]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[23]]$version
#> [1] 1
#> 
#> [[23]]$strand
#> [1] 1
#> 
#> [[23]]$gene_id
#> [1] "ENSG00000263426"
#> 
#> 
#> [[24]]
#> [[24]]$gene_id
#> [1] "ENSG00000227214"
#> 
#> [[24]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[24]]$version
#> [1] 2
#> 
#> [[24]]$strand
#> [1] 1
#> 
#> [[24]]$id
#> [1] "ENSG00000227214"
#> 
#> [[24]]$end
#> [1] 28955261
#> 
#> [[24]]$external_name
#> [1] "HCG15"
#> 
#> [[24]]$canonical_transcript
#> [1] "ENST00000438190.1"
#> 
#> [[24]]$start
#> [1] 28953980
#> 
#> [[24]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[24]]$feature_type
#> [1] "gene"
#> 
#> [[24]]$biotype
#> [1] "antisense"
#> 
#> [[24]]$description
#> [1] "HLA complex group 15 (non-protein coding) [Source:HGNC Symbol;Acc:18361]"
#> 
#> [[24]]$seq_region_name
#> [1] "6"
#> 
#> [[24]]$source
#> [1] "havana"
#> 
#> 
#> [[25]]
#> [[25]]$source
#> [1] "havana"
#> 
#> [[25]]$seq_region_name
#> [1] "6"
#> 
#> [[25]]$description
#> [1] "HLA complex group 16 (non-protein coding) [Source:HGNC Symbol;Acc:20424]"
#> 
#> [[25]]$biotype
#> [1] "antisense"
#> 
#> [[25]]$feature_type
#> [1] "gene"
#> 
#> [[25]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[25]]$start
#> [1] 28954577
#> 
#> [[25]]$canonical_transcript
#> [1] "ENST00000445544.1"
#> 
#> [[25]]$external_name
#> [1] "HCG16"
#> 
#> [[25]]$end
#> [1] 28956313
#> 
#> [[25]]$id
#> [1] "ENSG00000244349"
#> 
#> [[25]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[25]]$strand
#> [1] 1
#> 
#> [[25]]$version
#> [1] 1
#> 
#> [[25]]$gene_id
#> [1] "ENSG00000244349"
#> 
#> 
#> [[26]]
#> [[26]]$description
#> [1] "olfactory receptor, family 2, subfamily AD, member 1 pseudogene [Source:HGNC Symbol;Acc:14749]"
#> 
#> [[26]]$biotype
#> [1] "pseudogene"
#> 
#> [[26]]$feature_type
#> [1] "gene"
#> 
#> [[26]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[26]]$source
#> [1] "havana"
#> 
#> [[26]]$seq_region_name
#> [1] "6"
#> 
#> [[26]]$version
#> [1] 1
#> 
#> [[26]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[26]]$strand
#> [1] -1
#> 
#> [[26]]$gene_id
#> [1] "ENSG00000223677"
#> 
#> [[26]]$start
#> [1] 28994457
#> 
#> [[26]]$external_name
#> [1] "OR2AD1P"
#> 
#> [[26]]$canonical_transcript
#> [1] "ENST00000441992.1"
#> 
#> [[26]]$end
#> [1] 28995384
#> 
#> [[26]]$id
#> [1] "ENSG00000223677"
#> 
#> 
#> [[27]]
#> [[27]]$source
#> [1] "havana"
#> 
#> [[27]]$seq_region_name
#> [1] "6"
#> 
#> [[27]]$biotype
#> [1] "pseudogene"
#> 
#> [[27]]$feature_type
#> [1] "gene"
#> 
#> [[27]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[27]]$description
#> [1] "olfactory receptor, family 2, subfamily P, member 1 pseudogene [Source:HGNC Symbol;Acc:8272]"
#> 
#> [[27]]$external_name
#> [1] "OR2P1P"
#> 
#> [[27]]$canonical_transcript
#> [1] "ENST00000413771.1"
#> 
#> [[27]]$end
#> [1] 29040376
#> 
#> [[27]]$id
#> [1] "ENSG00000236909"
#> 
#> [[27]]$start
#> [1] 29039601
#> 
#> [[27]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[27]]$version
#> [1] 1
#> 
#> [[27]]$strand
#> [1] -1
#> 
#> [[27]]$gene_id
#> [1] "ENSG00000236909"
#> 
#> 
#> [[28]]
#> [[28]]$strand
#> [1] -1
#> 
#> [[28]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[28]]$version
#> [1] 1
#> 
#> [[28]]$gene_id
#> [1] "ENSG00000228364"
#> 
#> [[28]]$canonical_transcript
#> [1] "ENST00000455287.1"
#> 
#> [[28]]$external_name
#> [1] "SAR1P1"
#> 
#> [[28]]$id
#> [1] "ENSG00000228364"
#> 
#> [[28]]$end
#> [1] 29044945
#> 
#> [[28]]$start
#> [1] 29044350
#> 
#> [[28]]$biotype
#> [1] "pseudogene"
#> 
#> [[28]]$feature_type
#> [1] "gene"
#> 
#> [[28]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[28]]$description
#> NULL
#> 
#> [[28]]$source
#> [1] "havana"
#> 
#> [[28]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[29]]
#> [[29]]$seq_region_name
#> [1] "6"
#> 
#> [[29]]$source
#> [1] "havana"
#> 
#> [[29]]$description
#> NULL
#> 
#> [[29]]$feature_type
#> [1] "gene"
#> 
#> [[29]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[29]]$biotype
#> [1] "lincRNA"
#> 
#> [[29]]$start
#> [1] 29091987
#> 
#> [[29]]$end
#> [1] 29096685
#> 
#> [[29]]$id
#> [1] "ENSG00000227206"
#> 
#> [[29]]$external_name
#> [1] "XXbac-BCX196D17.5"
#> 
#> [[29]]$canonical_transcript
#> [1] "ENST00000419702.1"
#> 
#> [[29]]$gene_id
#> [1] "ENSG00000227206"
#> 
#> [[29]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[29]]$version
#> [1] 1
#> 
#> [[29]]$strand
#> [1] 1
#> 
#> 
#> [[30]]
#> [[30]]$biotype
#> [1] "pseudogene"
#> 
#> [[30]]$feature_type
#> [1] "gene"
#> 
#> [[30]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[30]]$description
#> [1] "olfactory receptor, family 2, subfamily N, member 1 pseudogene [Source:HGNC Symbol;Acc:8271]"
#> 
#> [[30]]$source
#> [1] "havana"
#> 
#> [[30]]$seq_region_name
#> [1] "6"
#> 
#> [[30]]$version
#> [1] 3
#> 
#> [[30]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[30]]$strand
#> [1] -1
#> 
#> [[30]]$gene_id
#> [1] "ENSG00000203492"
#> 
#> [[30]]$canonical_transcript
#> [1] "ENST00000421799.1"
#> 
#> [[30]]$external_name
#> [1] "OR2N1P"
#> 
#> [[30]]$end
#> [1] 29106593
#> 
#> [[30]]$id
#> [1] "ENSG00000203492"
#> 
#> [[30]]$start
#> [1] 29105657
#> 
#> 
#> [[31]]
#> [[31]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[31]]$feature_type
#> [1] "gene"
#> 
#> [[31]]$biotype
#> [1] "pseudogene"
#> 
#> [[31]]$description
#> [1] "olfactory receptor, family 2, subfamily J, member 4 pseudogene [Source:HGNC Symbol;Acc:8262]"
#> 
#> [[31]]$seq_region_name
#> [1] "6"
#> 
#> [[31]]$source
#> [1] "havana"
#> 
#> [[31]]$gene_id
#> [1] "ENSG00000224233"
#> 
#> [[31]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[31]]$version
#> [1] 1
#> 
#> [[31]]$strand
#> [1] 1
#> 
#> [[31]]$end
#> [1] 29150222
#> 
#> [[31]]$id
#> [1] "ENSG00000224233"
#> 
#> [[31]]$canonical_transcript
#> [1] "ENST00000453647.1"
#> 
#> [[31]]$external_name
#> [1] "OR2J4P"
#> 
#> [[31]]$start
#> [1] 29149287
#> 
#> 
#> [[32]]
#> [[32]]$start
#> [1] 29183013
#> 
#> [[32]]$end
#> [1] 29183956
#> 
#> [[32]]$id
#> [1] "ENSG00000230598"
#> 
#> [[32]]$external_name
#> [1] "OR2H4P"
#> 
#> [[32]]$canonical_transcript
#> [1] "ENST00000415595.1"
#> 
#> [[32]]$gene_id
#> [1] "ENSG00000230598"
#> 
#> [[32]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[32]]$version
#> [1] 1
#> 
#> [[32]]$strand
#> [1] 1
#> 
#> [[32]]$seq_region_name
#> [1] "6"
#> 
#> [[32]]$source
#> [1] "havana"
#> 
#> [[32]]$description
#> [1] "olfactory receptor, family 2, subfamily H, member 4 pseudogene [Source:HGNC Symbol;Acc:8255]"
#> 
#> [[32]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[32]]$feature_type
#> [1] "gene"
#> 
#> [[32]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[33]]
#> [[33]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[33]]$version
#> [1] 2
#> 
#> [[33]]$strand
#> [1] 1
#> 
#> [[33]]$gene_id
#> [1] "ENSG00000232505"
#> 
#> [[33]]$start
#> [1] 29191750
#> 
#> [[33]]$canonical_transcript
#> [1] "ENST00000441381.1"
#> 
#> [[33]]$external_name
#> [1] "XXbac-BPG308J9.3"
#> 
#> [[33]]$id
#> [1] "ENSG00000232505"
#> 
#> [[33]]$end
#> [1] 29258217
#> 
#> [[33]]$description
#> NULL
#> 
#> [[33]]$biotype
#> [1] "lincRNA"
#> 
#> [[33]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[33]]$feature_type
#> [1] "gene"
#> 
#> [[33]]$source
#> [1] "havana"
#> 
#> [[33]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[34]]
#> [[34]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[34]]$strand
#> [1] -1
#> 
#> [[34]]$version
#> [1] 2
#> 
#> [[34]]$gene_id
#> [1] "ENSG00000213911"
#> 
#> [[34]]$start
#> [1] 29197004
#> 
#> [[34]]$external_name
#> [1] "OR2G1P"
#> 
#> [[34]]$canonical_transcript
#> [1] "ENST00000396808.2"
#> 
#> [[34]]$id
#> [1] "ENSG00000213911"
#> 
#> [[34]]$end
#> [1] 29197934
#> 
#> [[34]]$description
#> [1] "olfactory receptor, family 2, subfamily G, member 1 pseudogene [Source:HGNC Symbol;Acc:8251]"
#> 
#> [[34]]$biotype
#> [1] "pseudogene"
#> 
#> [[34]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[34]]$feature_type
#> [1] "gene"
#> 
#> [[34]]$source
#> [1] "havana"
#> 
#> [[34]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[35]]
#> [[35]]$seq_region_name
#> [1] "6"
#> 
#> [[35]]$source
#> [1] "ensembl_havana"
#> 
#> [[35]]$description
#> [1] "olfactory receptor, family 2, subfamily U, member 1 pseudogene [Source:HGNC Symbol;Acc:8278]"
#> 
#> [[35]]$feature_type
#> [1] "gene"
#> 
#> [[35]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[35]]$biotype
#> [1] "pseudogene"
#> 
#> [[35]]$start
#> [1] 29230480
#> 
#> [[35]]$end
#> [1] 29231856
#> 
#> [[35]]$id
#> [1] "ENSG00000204697"
#> 
#> [[35]]$external_name
#> [1] "OR2U1P"
#> 
#> [[35]]$canonical_transcript
#> [1] "ENST00000377164.1"
#> 
#> [[35]]$gene_id
#> [1] "ENSG00000204697"
#> 
#> [[35]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[35]]$version
#> [1] 3
#> 
#> [[35]]$strand
#> [1] -1
#> 
#> 
#> [[36]]
#> [[36]]$strand
#> [1] -1
#> 
#> [[36]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[36]]$version
#> [1] 1
#> 
#> [[36]]$gene_id
#> [1] "ENSG00000242524"
#> 
#> [[36]]$start
#> [1] 29236239
#> 
#> [[36]]$canonical_transcript
#> [1] "ENST00000465951.1"
#> 
#> [[36]]$external_name
#> [1] "OR2U2P"
#> 
#> [[36]]$id
#> [1] "ENSG00000242524"
#> 
#> [[36]]$end
#> [1] 29237198
#> 
#> [[36]]$description
#> [1] "olfactory receptor, family 2, subfamily U, member 2 pseudogene [Source:HGNC Symbol;Acc:8279]"
#> 
#> [[36]]$biotype
#> [1] "pseudogene"
#> 
#> [[36]]$feature_type
#> [1] "gene"
#> 
#> [[36]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[36]]$source
#> [1] "havana"
#> 
#> [[36]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[37]]
#> [[37]]$description
#> [1] "olfactory receptor, family 2, subfamily B, member 4 pseudogene [Source:HGNC Symbol;Acc:8239]"
#> 
#> [[37]]$biotype
#> [1] "pseudogene"
#> 
#> [[37]]$feature_type
#> [1] "gene"
#> 
#> [[37]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[37]]$source
#> [1] "havana"
#> 
#> [[37]]$seq_region_name
#> [1] "6"
#> 
#> [[37]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[37]]$version
#> [1] 5
#> 
#> [[37]]$strand
#> [1] 1
#> 
#> [[37]]$gene_id
#> [1] "ENSG00000197171"
#> 
#> [[37]]$start
#> [1] 29258488
#> 
#> [[37]]$external_name
#> [1] "OR2B4P"
#> 
#> [[37]]$canonical_transcript
#> [1] "ENST00000356293.5"
#> 
#> [[37]]$end
#> [1] 29259430
#> 
#> [[37]]$id
#> [1] "ENSG00000197171"
#> 
#> 
#> [[38]]
#> [[38]]$source
#> [1] "havana"
#> 
#> [[38]]$seq_region_name
#> [1] "6"
#> 
#> [[38]]$biotype
#> [1] "pseudogene"
#> 
#> [[38]]$feature_type
#> [1] "gene"
#> 
#> [[38]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[38]]$description
#> [1] "DEAD (Asp-Glu-Ala-Asp) box helicase 6 pseudogene 1 [Source:HGNC Symbol;Acc:13948]"
#> 
#> [[38]]$canonical_transcript
#> [1] "ENST00000441966.1"
#> 
#> [[38]]$external_name
#> [1] "DDX6P1"
#> 
#> [[38]]$end
#> [1] 29298838
#> 
#> [[38]]$id
#> [1] "ENSG00000230056"
#> 
#> [[38]]$start
#> [1] 29297403
#> 
#> [[38]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[38]]$version
#> [1] 1
#> 
#> [[38]]$strand
#> [1] -1
#> 
#> [[38]]$gene_id
#> [1] "ENSG00000230056"
#> 
#> 
#> [[39]]
#> [[39]]$canonical_transcript
#> [1] "ENST00000457888.2"
#> 
#> [[39]]$external_name
#> [1] "UBDP1"
#> 
#> [[39]]$id
#> [1] "ENSG00000225797"
#> 
#> [[39]]$end
#> [1] 29437583
#> 
#> [[39]]$start
#> [1] 29432373
#> 
#> [[39]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[39]]$strand
#> [1] -1
#> 
#> [[39]]$version
#> [1] 2
#> 
#> [[39]]$gene_id
#> [1] "ENSG00000225797"
#> 
#> [[39]]$source
#> [1] "havana"
#> 
#> [[39]]$seq_region_name
#> [1] "6"
#> 
#> [[39]]$biotype
#> [1] "pseudogene"
#> 
#> [[39]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[39]]$feature_type
#> [1] "gene"
#> 
#> [[39]]$description
#> [1] "ubiquitin D pseudogene 1 [Source:HGNC Symbol;Acc:18796]"
#> 
#> 
#> [[40]]
#> [[40]]$seq_region_name
#> [1] "6"
#> 
#> [[40]]$source
#> [1] "havana"
#> 
#> [[40]]$feature_type
#> [1] "gene"
#> 
#> [[40]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[40]]$biotype
#> [1] "pseudogene"
#> 
#> [[40]]$description
#> [1] "MAS1 oncogene-like pseudogene 1 [Source:HGNC Symbol;Acc:33448]"
#> 
#> [[40]]$id
#> [1] "ENSG00000230164"
#> 
#> [[40]]$end
#> [1] 29443848
#> 
#> [[40]]$external_name
#> [1] "MAS1LP1"
#> 
#> [[40]]$canonical_transcript
#> [1] "ENST00000419416.1"
#> 
#> [[40]]$start
#> [1] 29442802
#> 
#> [[40]]$gene_id
#> [1] "ENSG00000230164"
#> 
#> [[40]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[40]]$strand
#> [1] -1
#> 
#> [[40]]$version
#> [1] 1
#> 
#> 
#> [[41]]
#> [[41]]$source
#> [1] "havana"
#> 
#> [[41]]$seq_region_name
#> [1] "6"
#> 
#> [[41]]$biotype
#> [1] "pseudogene"
#> 
#> [[41]]$feature_type
#> [1] "gene"
#> 
#> [[41]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[41]]$description
#> [1] "ribosomal protein S17 pseudogene 1 [Source:HGNC Symbol;Acc:13982]"
#> 
#> [[41]]$external_name
#> [1] "RPS17P1"
#> 
#> [[41]]$canonical_transcript
#> [1] "ENST00000396783.2"
#> 
#> [[41]]$id
#> [1] "ENSG00000213900"
#> 
#> [[41]]$end
#> [1] 29457453
#> 
#> [[41]]$start
#> [1] 29457048
#> 
#> [[41]]$version
#> [1] 2
#> 
#> [[41]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[41]]$strand
#> [1] -1
#> 
#> [[41]]$gene_id
#> [1] "ENSG00000213900"
#> 
#> 
#> [[42]]
#> [[42]]$gene_id
#> [1] "ENSG00000229274"
#> 
#> [[42]]$version
#> [1] 1
#> 
#> [[42]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[42]]$strand
#> [1] 1
#> 
#> [[42]]$start
#> [1] 29465286
#> 
#> [[42]]$end
#> [1] 29478333
#> 
#> [[42]]$id
#> [1] "ENSG00000229274"
#> 
#> [[42]]$canonical_transcript
#> [1] "ENST00000436804.1"
#> 
#> [[42]]$external_name
#> [1] "XXbac-BPG13B8.10"
#> 
#> [[42]]$description
#> NULL
#> 
#> [[42]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[42]]$feature_type
#> [1] "gene"
#> 
#> [[42]]$biotype
#> [1] "lincRNA"
#> 
#> [[42]]$seq_region_name
#> [1] "6"
#> 
#> [[42]]$source
#> [1] "havana"
#> 
#> 
#> [[43]]
#> [[43]]$seq_region_name
#> [1] "6"
#> 
#> [[43]]$source
#> [1] "havana"
#> 
#> [[43]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[43]]$feature_type
#> [1] "gene"
#> 
#> [[43]]$biotype
#> [1] "antisense"
#> 
#> [[43]]$description
#> [1] "long intergenic non-protein coding RNA 1015 [Source:HGNC Symbol;Acc:48988]"
#> 
#> [[43]]$id
#> [1] "ENSG00000224582"
#> 
#> [[43]]$end
#> [1] 29501345
#> 
#> [[43]]$external_name
#> [1] "LINC01015"
#> 
#> [[43]]$canonical_transcript
#> [1] "ENST00000434657.1"
#> 
#> [[43]]$start
#> [1] 29497197
#> 
#> [[43]]$gene_id
#> [1] "ENSG00000224582"
#> 
#> [[43]]$strand
#> [1] 1
#> 
#> [[43]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[43]]$version
#> [1] 1
#> 
#> 
#> [[44]]
#> [[44]]$feature_type
#> [1] "gene"
#> 
#> [[44]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[44]]$biotype
#> [1] "pseudogene"
#> 
#> [[44]]$description
#> [1] "G protein-coupled receptor 53, pseudogene [Source:HGNC Symbol;Acc:4509]"
#> 
#> [[44]]$seq_region_name
#> [1] "6"
#> 
#> [[44]]$source
#> [1] "havana"
#> 
#> [[44]]$gene_id
#> [1] "ENSG00000229281"
#> 
#> [[44]]$version
#> [1] 1
#> 
#> [[44]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[44]]$strand
#> [1] -1
#> 
#> [[44]]$id
#> [1] "ENSG00000229281"
#> 
#> [[44]]$end
#> [1] 29506564
#> 
#> [[44]]$canonical_transcript
#> [1] "ENST00000457298.1"
#> 
#> [[44]]$external_name
#> [1] "GPR53P"
#> 
#> [[44]]$start
#> [1] 29505481
#> 
#> 
#> [[45]]
#> [[45]]$seq_region_name
#> [1] "6"
#> 
#> [[45]]$source
#> [1] "ensembl_havana"
#> 
#> [[45]]$description
#> [1] "olfactory receptor, family 2, subfamily I, member 1 pseudogene [Source:HGNC Symbol;Acc:8258]"
#> 
#> [[45]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[45]]$feature_type
#> [1] "gene"
#> 
#> [[45]]$biotype
#> [1] "pseudogene"
#> 
#> [[45]]$start
#> [1] 29520996
#> 
#> [[45]]$id
#> [1] "ENSG00000237988"
#> 
#> [[45]]$end
#> [1] 29521943
#> 
#> [[45]]$external_name
#> [1] "OR2I1P"
#> 
#> [[45]]$canonical_transcript
#> [1] "ENST00000453522.1"
#> 
#> [[45]]$gene_id
#> [1] "ENSG00000237988"
#> 
#> [[45]]$version
#> [1] 2
#> 
#> [[45]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[45]]$strand
#> [1] 1
#> 
#> 
#> [[46]]
#> [[46]]$seq_region_name
#> [1] "6"
#> 
#> [[46]]$source
#> [1] "havana"
#> 
#> [[46]]$description
#> [1] "olfactory receptor, family 2, subfamily H, member 5 pseudogene [Source:HGNC Symbol;Acc:8256]"
#> 
#> [[46]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[46]]$feature_type
#> [1] "gene"
#> 
#> [[46]]$biotype
#> [1] "pseudogene"
#> 
#> [[46]]$start
#> [1] 29541686
#> 
#> [[46]]$end
#> [1] 29542613
#> 
#> [[46]]$id
#> [1] "ENSG00000232173"
#> 
#> [[46]]$external_name
#> [1] "OR2H5P"
#> 
#> [[46]]$canonical_transcript
#> [1] "ENST00000430565.1"
#> 
#> [[46]]$gene_id
#> [1] "ENSG00000232173"
#> 
#> [[46]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[46]]$strand
#> [1] 1
#> 
#> [[46]]$version
#> [1] 2
#> 
#> 
#> [[47]]
#> [[47]]$gene_id
#> [1] "ENSG00000227609"
#> 
#> [[47]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[47]]$version
#> [1] 2
#> 
#> [[47]]$strand
#> [1] 1
#> 
#> [[47]]$start
#> [1] 29545236
#> 
#> [[47]]$end
#> [1] 29545525
#> 
#> [[47]]$id
#> [1] "ENSG00000227609"
#> 
#> [[47]]$external_name
#> [1] "TMEM183AP1"
#> 
#> [[47]]$canonical_transcript
#> [1] "ENST00000432739.1"
#> 
#> [[47]]$description
#> [1] "transmembrane protein 183A pseudogene 1 [Source:HGNC Symbol;Acc:33458]"
#> 
#> [[47]]$feature_type
#> [1] "gene"
#> 
#> [[47]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[47]]$biotype
#> [1] "pseudogene"
#> 
#> [[47]]$seq_region_name
#> [1] "6"
#> 
#> [[47]]$source
#> [1] "havana"
#> 
#> 
#> [[48]]
#> [[48]]$gene_id
#> [1] "ENSG00000201330"
#> 
#> [[48]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[48]]$version
#> [1] 1
#> 
#> [[48]]$strand
#> [1] 1
#> 
#> [[48]]$id
#> [1] "ENSG00000201330"
#> 
#> [[48]]$end
#> [1] 29550109
#> 
#> [[48]]$external_name
#> [1] "SNORD32B"
#> 
#> [[48]]$canonical_transcript
#> [1] "ENST00000364460.1"
#> 
#> [[48]]$start
#> [1] 29550026
#> 
#> [[48]]$feature_type
#> [1] "gene"
#> 
#> [[48]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[48]]$biotype
#> [1] "snoRNA"
#> 
#> [[48]]$description
#> [1] "small nucleolar RNA, C/D box 32B [Source:HGNC Symbol;Acc:32719]"
#> 
#> [[48]]$seq_region_name
#> [1] "6"
#> 
#> [[48]]$source
#> [1] "ensembl"
#> 
#> 
#> [[49]]
#> [[49]]$gene_id
#> [1] "ENSG00000231301"
#> 
#> [[49]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[49]]$version
#> [1] 1
#> 
#> [[49]]$strand
#> [1] 1
#> 
#> [[49]]$start
#> [1] 29550285
#> 
#> [[49]]$id
#> [1] "ENSG00000231301"
#> 
#> [[49]]$end
#> [1] 29550802
#> 
#> [[49]]$external_name
#> [1] "RPL13AP"
#> 
#> [[49]]$canonical_transcript
#> [1] "ENST00000441979.1"
#> 
#> [[49]]$description
#> [1] "ribosomal protein L13a pseudogene [Source:HGNC Symbol;Acc:13977]"
#> 
#> [[49]]$feature_type
#> [1] "gene"
#> 
#> [[49]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[49]]$biotype
#> [1] "pseudogene"
#> 
#> [[49]]$seq_region_name
#> [1] "6"
#> 
#> [[49]]$source
#> [1] "havana"
#> 
#> 
#> [[50]]
#> [[50]]$strand
#> [1] -1
#> 
#> [[50]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[50]]$version
#> [1] 1
#> 
#> [[50]]$gene_id
#> [1] "ENSG00000235238"
#> 
#> [[50]]$start
#> [1] 29603837
#> 
#> [[50]]$external_name
#> [1] "SUMO2P1"
#> 
#> [[50]]$canonical_transcript
#> [1] "ENST00000445436.1"
#> 
#> [[50]]$end
#> [1] 29604120
#> 
#> [[50]]$id
#> [1] "ENSG00000235238"
#> 
#> [[50]]$description
#> [1] "SUMO2 pseudogene 1 [Source:HGNC Symbol;Acc:13985]"
#> 
#> [[50]]$biotype
#> [1] "pseudogene"
#> 
#> [[50]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[50]]$feature_type
#> [1] "gene"
#> 
#> [[50]]$source
#> [1] "havana"
#> 
#> [[50]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[51]]
#> [[51]]$version
#> [1] 1
#> 
#> [[51]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[51]]$strand
#> [1] -1
#> 
#> [[51]]$gene_id
#> [1] "ENSG00000233916"
#> 
#> [[51]]$external_name
#> [1] "ZDHHC20P1"
#> 
#> [[51]]$canonical_transcript
#> [1] "ENST00000435787.1"
#> 
#> [[51]]$end
#> [1] 29676324
#> 
#> [[51]]$id
#> [1] "ENSG00000233916"
#> 
#> [[51]]$start
#> [1] 29675902
#> 
#> [[51]]$biotype
#> [1] "pseudogene"
#> 
#> [[51]]$feature_type
#> [1] "gene"
#> 
#> [[51]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[51]]$description
#> [1] "zinc finger, DHHC-type containing 20 pseudogene 1 [Source:HGNC Symbol;Acc:33456]"
#> 
#> [[51]]$source
#> [1] "havana"
#> 
#> [[51]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[52]]
#> [[52]]$gene_id
#> [1] "ENSG00000225864"
#> 
#> [[52]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[52]]$strand
#> [1] -1
#> 
#> [[52]]$version
#> [1] 1
#> 
#> [[52]]$end
#> [1] 29691748
#> 
#> [[52]]$id
#> [1] "ENSG00000225864"
#> 
#> [[52]]$external_name
#> [1] "HCG4P11"
#> 
#> [[52]]$canonical_transcript
#> [1] "ENST00000427340.1"
#> 
#> [[52]]$start
#> [1] 29690758
#> 
#> [[52]]$feature_type
#> [1] "gene"
#> 
#> [[52]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[52]]$biotype
#> [1] "pseudogene"
#> 
#> [[52]]$description
#> [1] "HLA complex group 4 pseudogene 11 [Source:HGNC Symbol;Acc:22930]"
#> 
#> [[52]]$seq_region_name
#> [1] "6"
#> 
#> [[52]]$source
#> [1] "havana"
#> 
#> 
#> [[53]]
#> [[53]]$start
#> [1] 29694378
#> 
#> [[53]]$canonical_transcript
#> [1] "ENST00000458236.1"
#> 
#> [[53]]$external_name
#> [1] "HLA-F-AS1"
#> 
#> [[53]]$end
#> [1] 29716826
#> 
#> [[53]]$id
#> [1] "ENSG00000214922"
#> 
#> [[53]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[53]]$strand
#> [1] -1
#> 
#> [[53]]$version
#> [1] 5
#> 
#> [[53]]$gene_id
#> [1] "ENSG00000214922"
#> 
#> [[53]]$source
#> [1] "havana"
#> 
#> [[53]]$seq_region_name
#> [1] "6"
#> 
#> [[53]]$description
#> [1] "HLA-F antisense RNA 1 [Source:HGNC Symbol;Acc:26645]"
#> 
#> [[53]]$biotype
#> [1] "processed_transcript"
#> 
#> [[53]]$feature_type
#> [1] "gene"
#> 
#> [[53]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> 
#> [[54]]
#> [[54]]$seq_region_name
#> [1] "6"
#> 
#> [[54]]$source
#> [1] "havana"
#> 
#> [[54]]$feature_type
#> [1] "gene"
#> 
#> [[54]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[54]]$biotype
#> [1] "pseudogene"
#> 
#> [[54]]$description
#> [1] "ribosomal protein L23a pseudogene 1 [Source:HGNC Symbol;Acc:10318]"
#> 
#> [[54]]$id
#> [1] "ENSG00000239257"
#> 
#> [[54]]$end
#> [1] 29694916
#> 
#> [[54]]$canonical_transcript
#> [1] "ENST00000428990.1"
#> 
#> [[54]]$external_name
#> [1] "RPL23AP1"
#> 
#> [[54]]$start
#> [1] 29694446
#> 
#> [[54]]$gene_id
#> [1] "ENSG00000239257"
#> 
#> [[54]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[54]]$version
#> [1] 1
#> 
#> [[54]]$strand
#> [1] -1
#> 
#> 
#> [[55]]
#> [[55]]$start
#> [1] 29709508
#> 
#> [[55]]$end
#> [1] 29716746
#> 
#> [[55]]$id
#> [1] "ENSG00000273340"
#> 
#> [[55]]$external_name
#> [1] "MICE"
#> 
#> [[55]]$canonical_transcript
#> [1] "ENST00000510438.1"
#> 
#> [[55]]$gene_id
#> [1] "ENSG00000273340"
#> 
#> [[55]]$strand
#> [1] -1
#> 
#> [[55]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[55]]$version
#> [1] 1
#> 
#> [[55]]$seq_region_name
#> [1] "6"
#> 
#> [[55]]$source
#> [1] "havana"
#> 
#> [[55]]$description
#> [1] "MHC class I polypeptide-related sequence E (pseudogene) [Source:HGNC Symbol;Acc:7094]"
#> 
#> [[55]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[55]]$feature_type
#> [1] "gene"
#> 
#> [[55]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[56]]
#> [[56]]$start
#> [1] 29716066
#> 
#> [[56]]$external_name
#> [1] "HCG9P5"
#> 
#> [[56]]$canonical_transcript
#> [1] "ENST00000435408.1"
#> 
#> [[56]]$id
#> [1] "ENSG00000227758"
#> 
#> [[56]]$end
#> [1] 29716290
#> 
#> [[56]]$strand
#> [1] 1
#> 
#> [[56]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[56]]$version
#> [1] 1
#> 
#> [[56]]$gene_id
#> [1] "ENSG00000227758"
#> 
#> [[56]]$source
#> [1] "havana"
#> 
#> [[56]]$seq_region_name
#> [1] "6"
#> 
#> [[56]]$description
#> [1] "HLA complex group 9 pseudogene 5 [Source:HGNC Symbol;Acc:30983]"
#> 
#> [[56]]$biotype
#> [1] "pseudogene"
#> 
#> [[56]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[56]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[57]]
#> [[57]]$start
#> [1] 29718148
#> 
#> [[57]]$end
#> [1] 29718249
#> 
#> [[57]]$id
#> [1] "ENSG00000199290"
#> 
#> [[57]]$canonical_transcript
#> [1] "ENST00000362420.1"
#> 
#> [[57]]$external_name
#> [1] "Y_RNA"
#> 
#> [[57]]$gene_id
#> [1] "ENSG00000199290"
#> 
#> [[57]]$version
#> [1] 1
#> 
#> [[57]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[57]]$strand
#> [1] 1
#> 
#> [[57]]$seq_region_name
#> [1] "6"
#> 
#> [[57]]$source
#> [1] "ensembl"
#> 
#> [[57]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[57]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[57]]$feature_type
#> [1] "gene"
#> 
#> [[57]]$biotype
#> [1] "misc_RNA"
#> 
#> 
#> [[58]]
#> [[58]]$description
#> [1] "interferon induced transmembrane protein 4 pseudogene [Source:HGNC Symbol;Acc:21669]"
#> 
#> [[58]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[58]]$feature_type
#> [1] "gene"
#> 
#> [[58]]$biotype
#> [1] "pseudogene"
#> 
#> [[58]]$seq_region_name
#> [1] "6"
#> 
#> [[58]]$source
#> [1] "havana"
#> 
#> [[58]]$gene_id
#> [1] "ENSG00000235821"
#> 
#> [[58]]$version
#> [1] 1
#> 
#> [[58]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[58]]$strand
#> [1] -1
#> 
#> [[58]]$start
#> [1] 29718506
#> 
#> [[58]]$end
#> [1] 29718925
#> 
#> [[58]]$id
#> [1] "ENSG00000235821"
#> 
#> [[58]]$external_name
#> [1] "IFITM4P"
#> 
#> [[58]]$canonical_transcript
#> [1] "ENST00000441380.1"
#> 
#> 
#> [[59]]
#> [[59]]$source
#> [1] "havana"
#> 
#> [[59]]$seq_region_name
#> [1] "6"
#> 
#> [[59]]$biotype
#> [1] "lincRNA"
#> 
#> [[59]]$feature_type
#> [1] "gene"
#> 
#> [[59]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[59]]$description
#> NULL
#> 
#> [[59]]$external_name
#> [1] "XXbac-BPG170G13.32"
#> 
#> [[59]]$canonical_transcript
#> [1] "ENST00000606834.1"
#> 
#> [[59]]$id
#> [1] "ENSG00000272236"
#> 
#> [[59]]$end
#> [1] 29719984
#> 
#> [[59]]$start
#> [1] 29719742
#> 
#> [[59]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[59]]$strand
#> [1] -1
#> 
#> [[59]]$version
#> [1] 1
#> 
#> [[59]]$gene_id
#> [1] "ENSG00000272236"
#> 
#> 
#> [[60]]
#> [[60]]$seq_region_name
#> [1] "6"
#> 
#> [[60]]$source
#> [1] "havana"
#> 
#> [[60]]$description
#> NULL
#> 
#> [[60]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[60]]$feature_type
#> [1] "gene"
#> 
#> [[60]]$biotype
#> [1] "pseudogene"
#> 
#> [[60]]$start
#> [1] 29731035
#> 
#> [[60]]$id
#> [1] "ENSG00000270896"
#> 
#> [[60]]$end
#> [1] 29731169
#> 
#> [[60]]$external_name
#> [1] "XXbac-BPG170G13.31"
#> 
#> [[60]]$canonical_transcript
#> [1] "ENST00000604495.1"
#> 
#> [[60]]$gene_id
#> [1] "ENSG00000270896"
#> 
#> [[60]]$version
#> [1] 1
#> 
#> [[60]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[60]]$strand
#> [1] -1
#> 
#> 
#> [[61]]
#> [[61]]$source
#> [1] "ensembl_havana"
#> 
#> [[61]]$seq_region_name
#> [1] "6"
#> 
#> [[61]]$biotype
#> [1] "pseudogene"
#> 
#> [[61]]$feature_type
#> [1] "gene"
#> 
#> [[61]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[61]]$description
#> [1] "major histocompatibility complex, class I, V (pseudogene) [Source:HGNC Symbol;Acc:23482]"
#> 
#> [[61]]$external_name
#> [1] "HLA-V"
#> 
#> [[61]]$canonical_transcript
#> [1] "ENST00000476601.1"
#> 
#> [[61]]$id
#> [1] "ENSG00000181126"
#> 
#> [[61]]$end
#> [1] 29765588
#> 
#> [[61]]$start
#> [1] 29758731
#> 
#> [[61]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[61]]$strand
#> [1] 1
#> 
#> [[61]]$version
#> [1] 9
#> 
#> [[61]]$gene_id
#> [1] "ENSG00000181126"
#> 
#> 
#> [[62]]
#> [[62]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[62]]$version
#> [1] 3
#> 
#> [[62]]$strand
#> [1] -1
#> 
#> [[62]]$gene_id
#> [1] "ENSG00000176998"
#> 
#> [[62]]$start
#> [1] 29758816
#> 
#> [[62]]$canonical_transcript
#> [1] "ENST00000320533.2"
#> 
#> [[62]]$external_name
#> [1] "HCG4"
#> 
#> [[62]]$id
#> [1] "ENSG00000176998"
#> 
#> [[62]]$end
#> [1] 29760850
#> 
#> [[62]]$description
#> [1] "HLA complex group 4 (non-protein coding) [Source:HGNC Symbol;Acc:21241]"
#> 
#> [[62]]$biotype
#> [1] "pseudogene"
#> 
#> [[62]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[62]]$feature_type
#> [1] "gene"
#> 
#> [[62]]$source
#> [1] "ensembl_havana"
#> 
#> [[62]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[63]]
#> [[63]]$source
#> [1] "havana"
#> 
#> [[63]]$seq_region_name
#> [1] "6"
#> 
#> [[63]]$biotype
#> [1] "pseudogene"
#> 
#> [[63]]$feature_type
#> [1] "gene"
#> 
#> [[63]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[63]]$description
#> [1] "major histocompatibility complex, class I, P (pseudogene) [Source:HGNC Symbol;Acc:21196]"
#> 
#> [[63]]$external_name
#> [1] "HLA-P"
#> 
#> [[63]]$canonical_transcript
#> [1] "ENST00000568192.1"
#> 
#> [[63]]$id
#> [1] "ENSG00000261548"
#> 
#> [[63]]$end
#> [1] 29770202
#> 
#> [[63]]$start
#> [1] 29768192
#> 
#> [[63]]$version
#> [1] 1
#> 
#> [[63]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[63]]$strand
#> [1] 1
#> 
#> [[63]]$gene_id
#> [1] "ENSG00000261548"
#> 
#> 
#> [[64]]
#> [[64]]$description
#> [1] "ribosomal protein L7a pseudogene 7 [Source:HGNC Symbol;Acc:21395]"
#> 
#> [[64]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[64]]$feature_type
#> [1] "gene"
#> 
#> [[64]]$biotype
#> [1] "pseudogene"
#> 
#> [[64]]$seq_region_name
#> [1] "6"
#> 
#> [[64]]$source
#> [1] "havana"
#> 
#> [[64]]$gene_id
#> [1] "ENSG00000213880"
#> 
#> [[64]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[64]]$strand
#> [1] -1
#> 
#> [[64]]$version
#> [1] 3
#> 
#> [[64]]$start
#> [1] 29770972
#> 
#> [[64]]$end
#> [1] 29771768
#> 
#> [[64]]$id
#> [1] "ENSG00000213880"
#> 
#> [[64]]$canonical_transcript
#> [1] "ENST00000453429.1"
#> 
#> [[64]]$external_name
#> [1] "RPL7AP7"
#> 
#> 
#> [[65]]
#> [[65]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[65]]$strand
#> [1] -1
#> 
#> [[65]]$version
#> [1] 1
#> 
#> [[65]]$gene_id
#> [1] "ENSG00000237042"
#> 
#> [[65]]$start
#> [1] 29780342
#> 
#> [[65]]$canonical_transcript
#> [1] "ENST00000416822.1"
#> 
#> [[65]]$external_name
#> [1] "MICG"
#> 
#> [[65]]$end
#> [1] 29780473
#> 
#> [[65]]$id
#> [1] "ENSG00000237042"
#> 
#> [[65]]$description
#> [1] "MHC class I polypeptide-related sequence G (pseudogene) [Source:HGNC Symbol;Acc:16802]"
#> 
#> [[65]]$biotype
#> [1] "pseudogene"
#> 
#> [[65]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[65]]$feature_type
#> [1] "gene"
#> 
#> [[65]]$source
#> [1] "havana"
#> 
#> [[65]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[66]]
#> [[66]]$seq_region_name
#> [1] "6"
#> 
#> [[66]]$source
#> [1] "havana"
#> 
#> [[66]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[66]]$feature_type
#> [1] "gene"
#> 
#> [[66]]$biotype
#> [1] "pseudogene"
#> 
#> [[66]]$description
#> [1] "HLA complex group 4 pseudogene 8 [Source:HGNC Symbol;Acc:22927]"
#> 
#> [[66]]$end
#> [1] 29796141
#> 
#> [[66]]$id
#> [1] "ENSG00000229142"
#> 
#> [[66]]$external_name
#> [1] "HCG4P8"
#> 
#> [[66]]$canonical_transcript
#> [1] "ENST00000443049.1"
#> 
#> [[66]]$start
#> [1] 29795162
#> 
#> [[66]]$gene_id
#> [1] "ENSG00000229142"
#> 
#> [[66]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[66]]$version
#> [1] 1
#> 
#> [[66]]$strand
#> [1] -1
#> 
#> 
#> [[67]]
#> [[67]]$source
#> [1] "havana"
#> 
#> [[67]]$seq_region_name
#> [1] "6"
#> 
#> [[67]]$description
#> [1] "MHC class I polypeptide-related sequence F (pseudogene) [Source:HGNC Symbol;Acc:16801]"
#> 
#> [[67]]$biotype
#> [1] "pseudogene"
#> 
#> [[67]]$feature_type
#> [1] "gene"
#> 
#> [[67]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[67]]$start
#> [1] 29820140
#> 
#> [[67]]$canonical_transcript
#> [1] "ENST00000432679.1"
#> 
#> [[67]]$external_name
#> [1] "MICF"
#> 
#> [[67]]$id
#> [1] "ENSG00000233265"
#> 
#> [[67]]$end
#> [1] 29820268
#> 
#> [[67]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[67]]$version
#> [1] 1
#> 
#> [[67]]$strand
#> [1] -1
#> 
#> [[67]]$gene_id
#> [1] "ENSG00000233265"
#> 
#> 
#> [[68]]
#> [[68]]$gene_id
#> [1] "ENSG00000230521"
#> 
#> [[68]]$strand
#> [1] -1
#> 
#> [[68]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[68]]$version
#> [1] 1
#> 
#> [[68]]$id
#> [1] "ENSG00000230521"
#> 
#> [[68]]$end
#> [1] 29856045
#> 
#> [[68]]$external_name
#> [1] "HCG4P7"
#> 
#> [[68]]$canonical_transcript
#> [1] "ENST00000420084.1"
#> 
#> [[68]]$start
#> [1] 29855071
#> 
#> [[68]]$feature_type
#> [1] "gene"
#> 
#> [[68]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[68]]$biotype
#> [1] "pseudogene"
#> 
#> [[68]]$description
#> [1] "HLA complex group 4 pseudogene 7 [Source:HGNC Symbol;Acc:22926]"
#> 
#> [[68]]$seq_region_name
#> [1] "6"
#> 
#> [[68]]$source
#> [1] "havana"
#> 
#> 
#> [[69]]
#> [[69]]$source
#> [1] "ensembl_havana"
#> 
#> [[69]]$seq_region_name
#> [1] "6"
#> 
#> [[69]]$biotype
#> [1] "pseudogene"
#> 
#> [[69]]$feature_type
#> [1] "gene"
#> 
#> [[69]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[69]]$description
#> [1] "major histocompatibility complex, class I, H (pseudogene) [Source:HGNC Symbol;Acc:4965]"
#> 
#> [[69]]$external_name
#> [1] "HLA-H"
#> 
#> [[69]]$canonical_transcript
#> [1] "ENST00000383326.4"
#> 
#> [[69]]$id
#> [1] "ENSG00000206341"
#> 
#> [[69]]$end
#> [1] 29858259
#> 
#> [[69]]$start
#> [1] 29855350
#> 
#> [[69]]$version
#> [1] 6
#> 
#> [[69]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[69]]$strand
#> [1] 1
#> 
#> [[69]]$gene_id
#> [1] "ENSG00000206341"
#> 
#> 
#> [[70]]
#> [[70]]$feature_type
#> [1] "gene"
#> 
#> [[70]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[70]]$biotype
#> [1] "pseudogene"
#> 
#> [[70]]$description
#> [1] "major histocompatibility complex, class I, T (pseudogene) [Source:HGNC Symbol;Acc:23478]"
#> 
#> [[70]]$seq_region_name
#> [1] "6"
#> 
#> [[70]]$source
#> [1] "havana"
#> 
#> [[70]]$gene_id
#> [1] "ENSG00000231130"
#> 
#> [[70]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[70]]$version
#> [1] 1
#> 
#> [[70]]$strand
#> [1] 1
#> 
#> [[70]]$end
#> [1] 29865563
#> 
#> [[70]]$id
#> [1] "ENSG00000231130"
#> 
#> [[70]]$external_name
#> [1] "HLA-T"
#> 
#> [[70]]$canonical_transcript
#> [1] "ENST00000429813.1"
#> 
#> [[70]]$start
#> [1] 29864431
#> 
#> 
#> [[71]]
#> [[71]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[71]]$strand
#> [1] 1
#> 
#> [[71]]$version
#> [1] 1
#> 
#> [[71]]$gene_id
#> [1] "ENSG00000233677"
#> 
#> [[71]]$start
#> [1] 29874320
#> 
#> [[71]]$external_name
#> [1] "DDX39BP1"
#> 
#> [[71]]$canonical_transcript
#> [1] "ENST00000422507.1"
#> 
#> [[71]]$end
#> [1] 29874686
#> 
#> [[71]]$id
#> [1] "ENSG00000233677"
#> 
#> [[71]]$description
#> [1] "DEAD (Asp-Glu-Ala-Asp) box polypeptide 39B pseudogene 1 [Source:HGNC Symbol;Acc:33450]"
#> 
#> [[71]]$biotype
#> [1] "pseudogene"
#> 
#> [[71]]$feature_type
#> [1] "gene"
#> 
#> [[71]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[71]]$source
#> [1] "havana"
#> 
#> [[71]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[72]]
#> [[72]]$description
#> [1] "mitochondrial coiled-coil domain 1 pseudogene 1 [Source:HGNC Symbol;Acc:33449]"
#> 
#> [[72]]$biotype
#> [1] "pseudogene"
#> 
#> [[72]]$feature_type
#> [1] "gene"
#> 
#> [[72]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[72]]$source
#> [1] "havana"
#> 
#> [[72]]$seq_region_name
#> [1] "6"
#> 
#> [[72]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[72]]$version
#> [1] 1
#> 
#> [[72]]$strand
#> [1] -1
#> 
#> [[72]]$gene_id
#> [1] "ENSG00000235963"
#> 
#> [[72]]$start
#> [1] 29875560
#> 
#> [[72]]$external_name
#> [1] "MCCD1P1"
#> 
#> [[72]]$canonical_transcript
#> [1] "ENST00000445732.1"
#> 
#> [[72]]$end
#> [1] 29876422
#> 
#> [[72]]$id
#> [1] "ENSG00000235963"
#> 
#> 
#> [[73]]
#> [[73]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[73]]$feature_type
#> [1] "gene"
#> 
#> [[73]]$biotype
#> [1] "pseudogene"
#> 
#> [[73]]$description
#> [1] "HLA complex group 4B (non-protein coding) [Source:HGNC Symbol;Acc:22919]"
#> 
#> [[73]]$seq_region_name
#> [1] "6"
#> 
#> [[73]]$source
#> [1] "havana"
#> 
#> [[73]]$gene_id
#> [1] "ENSG00000227262"
#> 
#> [[73]]$strand
#> [1] -1
#> 
#> [[73]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[73]]$version
#> [1] 3
#> 
#> [[73]]$id
#> [1] "ENSG00000227262"
#> 
#> [[73]]$end
#> [1] 29894750
#> 
#> [[73]]$external_name
#> [1] "HCG4B"
#> 
#> [[73]]$canonical_transcript
#> [1] "ENST00000450128.1"
#> 
#> [[73]]$start
#> [1] 29893760
#> 
#> 
#> [[74]]
#> [[74]]$strand
#> [1] 1
#> 
#> [[74]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[74]]$version
#> [1] 2
#> 
#> [[74]]$gene_id
#> [1] "ENSG00000230795"
#> 
#> [[74]]$start
#> [1] 29894236
#> 
#> [[74]]$external_name
#> [1] "HLA-K"
#> 
#> [[74]]$canonical_transcript
#> [1] "ENST00000430151.1"
#> 
#> [[74]]$end
#> [1] 29897009
#> 
#> [[74]]$id
#> [1] "ENSG00000230795"
#> 
#> [[74]]$description
#> [1] "major histocompatibility complex, class I, K (pseudogene) [Source:HGNC Symbol;Acc:4969]"
#> 
#> [[74]]$biotype
#> [1] "pseudogene"
#> 
#> [[74]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[74]]$feature_type
#> [1] "gene"
#> 
#> [[74]]$source
#> [1] "havana"
#> 
#> [[74]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[75]]
#> [[75]]$biotype
#> [1] "pseudogene"
#> 
#> [[75]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[75]]$feature_type
#> [1] "gene"
#> 
#> [[75]]$description
#> [1] "major histocompatibility complex, class I, U (pseudogene) [Source:HGNC Symbol;Acc:23477]"
#> 
#> [[75]]$source
#> [1] "havana"
#> 
#> [[75]]$seq_region_name
#> [1] "6"
#> 
#> [[75]]$strand
#> [1] 1
#> 
#> [[75]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[75]]$version
#> [1] 1
#> 
#> [[75]]$gene_id
#> [1] "ENSG00000228078"
#> 
#> [[75]]$canonical_transcript
#> [1] "ENST00000418981.1"
#> 
#> [[75]]$external_name
#> [1] "HLA-U"
#> 
#> [[75]]$id
#> [1] "ENSG00000228078"
#> 
#> [[75]]$end
#> [1] 29902063
#> 
#> [[75]]$start
#> [1] 29901878
#> 
#> 
#> [[76]]
#> [[76]]$external_name
#> [1] "HCG4P5"
#> 
#> [[76]]$canonical_transcript
#> [1] "ENST00000429656.1"
#> 
#> [[76]]$id
#> [1] "ENSG00000227766"
#> 
#> [[76]]$end
#> [1] 29910844
#> 
#> [[76]]$start
#> [1] 29909852
#> 
#> [[76]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[76]]$version
#> [1] 1
#> 
#> [[76]]$strand
#> [1] -1
#> 
#> [[76]]$gene_id
#> [1] "ENSG00000227766"
#> 
#> [[76]]$source
#> [1] "havana"
#> 
#> [[76]]$seq_region_name
#> [1] "6"
#> 
#> [[76]]$biotype
#> [1] "pseudogene"
#> 
#> [[76]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[76]]$feature_type
#> [1] "gene"
#> 
#> [[76]]$description
#> [1] "HLA complex group 4 pseudogene 5 [Source:HGNC Symbol;Acc:22925]"
#> 
#> 
#> [[77]]
#> [[77]]$gene_id
#> [1] "ENSG00000235290"
#> 
#> [[77]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[77]]$strand
#> [1] 1
#> 
#> [[77]]$version
#> [1] 1
#> 
#> [[77]]$start
#> [1] 29924373
#> 
#> [[77]]$end
#> [1] 29926347
#> 
#> [[77]]$id
#> [1] "ENSG00000235290"
#> 
#> [[77]]$canonical_transcript
#> [1] "ENST00000439514.1"
#> 
#> [[77]]$external_name
#> [1] "HLA-W"
#> 
#> [[77]]$description
#> [1] "major histocompatibility complex, class I, W (pseudogene) [Source:HGNC Symbol;Acc:23425]"
#> 
#> [[77]]$feature_type
#> [1] "gene"
#> 
#> [[77]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[77]]$biotype
#> [1] "pseudogene"
#> 
#> [[77]]$seq_region_name
#> [1] "6"
#> 
#> [[77]]$source
#> [1] "havana"
#> 
#> 
#> [[78]]
#> [[78]]$source
#> [1] "havana"
#> 
#> [[78]]$seq_region_name
#> [1] "6"
#> 
#> [[78]]$description
#> [1] "MHC class I polypeptide-related sequence D (pseudogene) [Source:HGNC Symbol;Acc:7093]"
#> 
#> [[78]]$biotype
#> [1] "pseudogene"
#> 
#> [[78]]$feature_type
#> [1] "gene"
#> 
#> [[78]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[78]]$start
#> [1] 29938578
#> 
#> [[78]]$external_name
#> [1] "MICD"
#> 
#> [[78]]$canonical_transcript
#> [1] "ENST00000413248.1"
#> 
#> [[78]]$id
#> [1] "ENSG00000229390"
#> 
#> [[78]]$end
#> [1] 29940241
#> 
#> [[78]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[78]]$version
#> [1] 1
#> 
#> [[78]]$strand
#> [1] -1
#> 
#> [[78]]$gene_id
#> [1] "ENSG00000229390"
#> 
#> 
#> [[79]]
#> [[79]]$canonical_transcript
#> [1] "ENST00000376800.3"
#> 
#> [[79]]$external_name
#> [1] "HCG9"
#> 
#> [[79]]$id
#> [1] "ENSG00000204625"
#> 
#> [[79]]$end
#> [1] 29946183
#> 
#> [[79]]$start
#> [1] 29942889
#> 
#> [[79]]$version
#> [1] 6
#> 
#> [[79]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[79]]$strand
#> [1] 1
#> 
#> [[79]]$gene_id
#> [1] "ENSG00000204625"
#> 
#> [[79]]$source
#> [1] "havana"
#> 
#> [[79]]$seq_region_name
#> [1] "6"
#> 
#> [[79]]$biotype
#> [1] "lincRNA"
#> 
#> [[79]]$feature_type
#> [1] "gene"
#> 
#> [[79]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[79]]$description
#> [1] "HLA complex group 9 (non-protein coding) [Source:HGNC Symbol;Acc:21243]"
#> 
#> 
#> [[80]]
#> [[80]]$gene_id
#> [1] "ENSG00000238024"
#> 
#> [[80]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[80]]$strand
#> [1] 1
#> 
#> [[80]]$version
#> [1] 1
#> 
#> [[80]]$id
#> [1] "ENSG00000238024"
#> 
#> [[80]]$end
#> [1] 29961382
#> 
#> [[80]]$external_name
#> [1] "DDX39BP2"
#> 
#> [[80]]$canonical_transcript
#> [1] "ENST00000453923.1"
#> 
#> [[80]]$start
#> [1] 29960986
#> 
#> [[80]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[80]]$feature_type
#> [1] "gene"
#> 
#> [[80]]$biotype
#> [1] "pseudogene"
#> 
#> [[80]]$description
#> [1] "DEAD (Asp-Glu-Ala-Asp) box polypeptide 39B pseudogene 2 [Source:HGNC Symbol;Acc:33461]"
#> 
#> [[80]]$seq_region_name
#> [1] "6"
#> 
#> [[80]]$source
#> [1] "havana"
#> 
#> 
#> [[81]]
#> [[81]]$feature_type
#> [1] "gene"
#> 
#> [[81]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[81]]$biotype
#> [1] "pseudogene"
#> 
#> [[81]]$description
#> [1] "mitochondrial coiled-coil domain 1 pseudogene 2 [Source:HGNC Symbol;Acc:33462]"
#> 
#> [[81]]$seq_region_name
#> [1] "6"
#> 
#> [[81]]$source
#> [1] "havana"
#> 
#> [[81]]$gene_id
#> [1] "ENSG00000224312"
#> 
#> [[81]]$strand
#> [1] -1
#> 
#> [[81]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[81]]$version
#> [1] 1
#> 
#> [[81]]$id
#> [1] "ENSG00000224312"
#> 
#> [[81]]$end
#> [1] 29963091
#> 
#> [[81]]$canonical_transcript
#> [1] "ENST00000423228.1"
#> 
#> [[81]]$external_name
#> [1] "MCCD1P2"
#> 
#> [[81]]$start
#> [1] 29962214
#> 
#> 
#> [[82]]
#> [[82]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[82]]$feature_type
#> [1] "gene"
#> 
#> [[82]]$biotype
#> [1] "antisense"
#> 
#> [[82]]$description
#> [1] "ZNRD1 antisense RNA 1 [Source:HGNC Symbol;Acc:13924]"
#> 
#> [[82]]$seq_region_name
#> [1] "6"
#> 
#> [[82]]$source
#> [1] "ensembl_havana"
#> 
#> [[82]]$gene_id
#> [1] "ENSG00000204623"
#> 
#> [[82]]$version
#> [1] 4
#> 
#> [[82]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[82]]$strand
#> [1] -1
#> 
#> [[82]]$end
#> [1] 30029417
#> 
#> [[82]]$id
#> [1] "ENSG00000204623"
#> 
#> [[82]]$canonical_transcript
#> [1] "ENST00000376797.3"
#> 
#> [[82]]$external_name
#> [1] "ZNRD1-AS1"
#> 
#> [[82]]$start
#> [1] 29968788
#> 
#> 
#> [[83]]
#> [[83]]$gene_id
#> [1] "ENSG00000237669"
#> 
#> [[83]]$version
#> [1] 1
#> 
#> [[83]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[83]]$strand
#> [1] -1
#> 
#> [[83]]$start
#> [1] 29973898
#> 
#> [[83]]$end
#> [1] 29974893
#> 
#> [[83]]$id
#> [1] "ENSG00000237669"
#> 
#> [[83]]$canonical_transcript
#> [1] "ENST00000458060.1"
#> 
#> [[83]]$external_name
#> [1] "HCG4P3"
#> 
#> [[83]]$description
#> [1] "HLA complex group 4 pseudogene 3 [Source:HGNC Symbol;Acc:22922]"
#> 
#> [[83]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[83]]$feature_type
#> [1] "gene"
#> 
#> [[83]]$biotype
#> [1] "pseudogene"
#> 
#> [[83]]$seq_region_name
#> [1] "6"
#> 
#> [[83]]$source
#> [1] "havana"
#> 
#> 
#> [[84]]
#> [[84]]$strand
#> [1] 1
#> 
#> [[84]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[84]]$version
#> [1] 6
#> 
#> [[84]]$gene_id
#> [1] "ENSG00000204622"
#> 
#> [[84]]$start
#> [1] 29974360
#> 
#> [[84]]$canonical_transcript
#> [1] "ENST00000462773.1"
#> 
#> [[84]]$external_name
#> [1] "HLA-J"
#> 
#> [[84]]$end
#> [1] 29977733
#> 
#> [[84]]$id
#> [1] "ENSG00000204622"
#> 
#> [[84]]$description
#> [1] "major histocompatibility complex, class I, J (pseudogene) [Source:HGNC Symbol;Acc:4967]"
#> 
#> [[84]]$biotype
#> [1] "pseudogene"
#> 
#> [[84]]$feature_type
#> [1] "gene"
#> 
#> [[84]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[84]]$source
#> [1] "havana"
#> 
#> [[84]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[85]]
#> [[85]]$source
#> [1] "havana"
#> 
#> [[85]]$seq_region_name
#> [1] "6"
#> 
#> [[85]]$description
#> [1] "eukaryotic translation termination factor 1 pseudogene 1 [Source:HGNC Symbol;Acc:3478]"
#> 
#> [[85]]$biotype
#> [1] "pseudogene"
#> 
#> [[85]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[85]]$feature_type
#> [1] "gene"
#> 
#> [[85]]$start
#> [1] 29999490
#> 
#> [[85]]$external_name
#> [1] "ETF1P1"
#> 
#> [[85]]$canonical_transcript
#> [1] "ENST00000427296.1"
#> 
#> [[85]]$end
#> [1] 30000780
#> 
#> [[85]]$id
#> [1] "ENSG00000232757"
#> 
#> [[85]]$version
#> [1] 1
#> 
#> [[85]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[85]]$strand
#> [1] 1
#> 
#> [[85]]$gene_id
#> [1] "ENSG00000232757"
#> 
#> 
#> [[86]]
#> [[86]]$start
#> [1] 30056598
#> 
#> [[86]]$id
#> [1] "ENSG00000252901"
#> 
#> [[86]]$end
#> [1] 30056730
#> 
#> [[86]]$canonical_transcript
#> [1] "ENST00000517092.1"
#> 
#> [[86]]$external_name
#> [1] "AL669914.1"
#> 
#> [[86]]$gene_id
#> [1] "ENSG00000252901"
#> 
#> [[86]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[86]]$strand
#> [1] -1
#> 
#> [[86]]$version
#> [1] 1
#> 
#> [[86]]$seq_region_name
#> [1] "6"
#> 
#> [[86]]$source
#> [1] "ensembl"
#> 
#> [[86]]$description
#> NULL
#> 
#> [[86]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[86]]$feature_type
#> [1] "gene"
#> 
#> [[86]]$biotype
#> [1] "miRNA"
#> 
#> 
#> [[87]]
#> [[87]]$seq_region_name
#> [1] "6"
#> 
#> [[87]]$source
#> [1] "havana"
#> 
#> [[87]]$description
#> [1] "TRIM31 antisense RNA 1 [Source:HGNC Symbol;Acc:39761]"
#> 
#> [[87]]$feature_type
#> [1] "gene"
#> 
#> [[87]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[87]]$biotype
#> [1] "antisense"
#> 
#> [[87]]$start
#> [1] 30073017
#> 
#> [[87]]$id
#> [1] "ENSG00000231226"
#> 
#> [[87]]$end
#> [1] 30082501
#> 
#> [[87]]$external_name
#> [1] "TRIM31-AS1"
#> 
#> [[87]]$canonical_transcript
#> [1] "ENST00000440874.1"
#> 
#> [[87]]$gene_id
#> [1] "ENSG00000231226"
#> 
#> [[87]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[87]]$strand
#> [1] 1
#> 
#> [[87]]$version
#> [1] 1
#> 
#> 
#> [[88]]
#> [[88]]$biotype
#> [1] "snoRNA"
#> 
#> [[88]]$feature_type
#> [1] "gene"
#> 
#> [[88]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[88]]$description
#> [1] "Small nucleolar RNA SNORA48 [Source:RFAM;Acc:RF00554]"
#> 
#> [[88]]$source
#> [1] "ensembl"
#> 
#> [[88]]$seq_region_name
#> [1] "6"
#> 
#> [[88]]$version
#> [1] 1
#> 
#> [[88]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[88]]$strand
#> [1] -1
#> 
#> [[88]]$gene_id
#> [1] "ENSG00000252228"
#> 
#> [[88]]$external_name
#> [1] "SNORA48"
#> 
#> [[88]]$canonical_transcript
#> [1] "ENST00000516419.1"
#> 
#> [[88]]$id
#> [1] "ENSG00000252228"
#> 
#> [[88]]$end
#> [1] 30100743
#> 
#> [[88]]$start
#> [1] 30100582
#> 
#> 
#> [[89]]
#> [[89]]$gene_id
#> [1] "ENSG00000233892"
#> 
#> [[89]]$strand
#> [1] -1
#> 
#> [[89]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[89]]$version
#> [1] 1
#> 
#> [[89]]$start
#> [1] 30154575
#> 
#> [[89]]$end
#> [1] 30156391
#> 
#> [[89]]$id
#> [1] "ENSG00000233892"
#> 
#> [[89]]$canonical_transcript
#> [1] "ENST00000446875.1"
#> 
#> [[89]]$external_name
#> [1] "PAIP1P1"
#> 
#> [[89]]$description
#> [1] "poly(A) binding protein interacting protein 1 pseudogene 1 [Source:HGNC Symbol;Acc:18240]"
#> 
#> [[89]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[89]]$feature_type
#> [1] "gene"
#> 
#> [[89]]$biotype
#> [1] "pseudogene"
#> 
#> [[89]]$seq_region_name
#> [1] "6"
#> 
#> [[89]]$source
#> [1] "havana"
#> 
#> 
#> [[90]]
#> [[90]]$feature_type
#> [1] "gene"
#> 
#> [[90]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[90]]$biotype
#> [1] "lincRNA"
#> 
#> [[90]]$description
#> [1] "HLA complex group 17 (non-protein coding) [Source:HGNC Symbol;Acc:31339]"
#> 
#> [[90]]$seq_region_name
#> [1] "6"
#> 
#> [[90]]$source
#> [1] "havana"
#> 
#> [[90]]$gene_id
#> [1] "ENSG00000270604"
#> 
#> [[90]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[90]]$strand
#> [1] -1
#> 
#> [[90]]$version
#> [1] 1
#> 
#> [[90]]$id
#> [1] "ENSG00000270604"
#> 
#> [[90]]$end
#> [1] 30293911
#> 
#> [[90]]$external_name
#> [1] "HCG17"
#> 
#> [[90]]$canonical_transcript
#> [1] "ENST00000453558.1"
#> 
#> [[90]]$start
#> [1] 30201816
#> 
#> 
#> [[91]]
#> [[91]]$feature_type
#> [1] "gene"
#> 
#> [[91]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[91]]$biotype
#> [1] "pseudogene"
#> 
#> [[91]]$description
#> [1] "tripartite motif containing 26B, pseudogene [Source:HGNC Symbol;Acc:31338]"
#> 
#> [[91]]$seq_region_name
#> [1] "6"
#> 
#> [[91]]$source
#> [1] "havana"
#> 
#> [[91]]$gene_id
#> [1] "ENSG00000236475"
#> 
#> [[91]]$version
#> [1] 1
#> 
#> [[91]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[91]]$strand
#> [1] 1
#> 
#> [[91]]$id
#> [1] "ENSG00000236475"
#> 
#> [[91]]$end
#> [1] 30210056
#> 
#> [[91]]$external_name
#> [1] "TRIM26BP"
#> 
#> [[91]]$canonical_transcript
#> [1] "ENST00000427723.1"
#> 
#> [[91]]$start
#> [1] 30206078
#> 
#> 
#> [[92]]
#> [[92]]$source
#> [1] "ensembl_havana"
#> 
#> [[92]]$seq_region_name
#> [1] "6"
#> 
#> [[92]]$description
#> [1] "major histocompatibility complex, class I, L (pseudogene) [Source:HGNC Symbol;Acc:4970]"
#> 
#> [[92]]$biotype
#> [1] "pseudogene"
#> 
#> [[92]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[92]]$feature_type
#> [1] "gene"
#> 
#> [[92]]$start
#> [1] 30227361
#> 
#> [[92]]$canonical_transcript
#> [1] "ENST00000463348.1"
#> 
#> [[92]]$external_name
#> [1] "HLA-L"
#> 
#> [[92]]$id
#> [1] "ENSG00000243753"
#> 
#> [[92]]$end
#> [1] 30260791
#> 
#> [[92]]$version
#> [1] 1
#> 
#> [[92]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[92]]$strand
#> [1] 1
#> 
#> [[92]]$gene_id
#> [1] "ENSG00000243753"
#> 
#> 
#> [[93]]
#> [[93]]$gene_id
#> [1] "ENSG00000231074"
#> 
#> [[93]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[93]]$version
#> [1] 4
#> 
#> [[93]]$strand
#> [1] -1
#> 
#> [[93]]$start
#> [1] 30255174
#> 
#> [[93]]$id
#> [1] "ENSG00000231074"
#> 
#> [[93]]$end
#> [1] 30294927
#> 
#> [[93]]$external_name
#> [1] "HCG18"
#> 
#> [[93]]$canonical_transcript
#> [1] "ENST00000426882.1"
#> 
#> [[93]]$description
#> [1] "HLA complex group 18 (non-protein coding) [Source:HGNC Symbol;Acc:31337]"
#> 
#> [[93]]$feature_type
#> [1] "gene"
#> 
#> [[93]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[93]]$biotype
#> [1] "antisense"
#> 
#> [[93]]$seq_region_name
#> [1] "6"
#> 
#> [[93]]$source
#> [1] "havana"
#> 
#> 
#> [[94]]
#> [[94]]$version
#> [1] 1
#> 
#> [[94]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[94]]$strand
#> [1] 1
#> 
#> [[94]]$gene_id
#> [1] "ENSG00000224372"
#> 
#> [[94]]$start
#> [1] 30319193
#> 
#> [[94]]$canonical_transcript
#> [1] "ENST00000437516.1"
#> 
#> [[94]]$external_name
#> [1] "HLA-N"
#> 
#> [[94]]$end
#> [1] 30319327
#> 
#> [[94]]$id
#> [1] "ENSG00000224372"
#> 
#> [[94]]$description
#> [1] "major histocompatibility complex, class I, N (pseudogene) [Source:HGNC Symbol;Acc:19406]"
#> 
#> [[94]]$biotype
#> [1] "pseudogene"
#> 
#> [[94]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[94]]$feature_type
#> [1] "gene"
#> 
#> [[94]]$source
#> [1] "havana"
#> 
#> [[94]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[95]]
#> [[95]]$feature_type
#> [1] "gene"
#> 
#> [[95]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[95]]$biotype
#> [1] "pseudogene"
#> 
#> [[95]]$description
#> [1] "HLA complex group 19 pseudogene [Source:HGNC Symbol;Acc:31336]"
#> 
#> [[95]]$seq_region_name
#> [1] "6"
#> 
#> [[95]]$source
#> [1] "havana"
#> 
#> [[95]]$gene_id
#> [1] "ENSG00000224486"
#> 
#> [[95]]$version
#> [1] 1
#> 
#> [[95]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[95]]$strand
#> [1] -1
#> 
#> [[95]]$id
#> [1] "ENSG00000224486"
#> 
#> [[95]]$end
#> [1] 30327688
#> 
#> [[95]]$external_name
#> [1] "HCG19P"
#> 
#> [[95]]$canonical_transcript
#> [1] "ENST00000448756.1"
#> 
#> [[95]]$start
#> [1] 30327055
#> 
#> 
#> [[96]]
#> [[96]]$start
#> [1] 30330889
#> 
#> [[96]]$end
#> [1] 30332057
#> 
#> [[96]]$id
#> [1] "ENSG00000236405"
#> 
#> [[96]]$canonical_transcript
#> [1] "ENST00000441056.1"
#> 
#> [[96]]$external_name
#> [1] "UBQLN1P1"
#> 
#> [[96]]$gene_id
#> [1] "ENSG00000236405"
#> 
#> [[96]]$version
#> [1] 1
#> 
#> [[96]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[96]]$strand
#> [1] -1
#> 
#> [[96]]$seq_region_name
#> [1] "6"
#> 
#> [[96]]$source
#> [1] "havana"
#> 
#> [[96]]$description
#> [1] "ubiquilin 1 pseudogene 1 [Source:HGNC Symbol;Acc:21633]"
#> 
#> [[96]]$feature_type
#> [1] "gene"
#> 
#> [[96]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[96]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[97]]
#> [[97]]$seq_region_name
#> [1] "6"
#> 
#> [[97]]$source
#> [1] "havana"
#> 
#> [[97]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[97]]$feature_type
#> [1] "gene"
#> 
#> [[97]]$biotype
#> [1] "pseudogene"
#> 
#> [[97]]$description
#> [1] "MHC class I polypeptide-related sequence C (pseudogene) [Source:HGNC Symbol;Acc:7092]"
#> 
#> [[97]]$id
#> [1] "ENSG00000226577"
#> 
#> [[97]]$end
#> [1] 30387096
#> 
#> [[97]]$external_name
#> [1] "MICC"
#> 
#> [[97]]$canonical_transcript
#> [1] "ENST00000445710.1"
#> 
#> [[97]]$start
#> [1] 30382492
#> 
#> [[97]]$gene_id
#> [1] "ENSG00000226577"
#> 
#> [[97]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[97]]$version
#> [1] 1
#> 
#> [[97]]$strand
#> [1] 1
#> 
#> 
#> [[98]]
#> [[98]]$seq_region_name
#> [1] "6"
#> 
#> [[98]]$source
#> [1] "havana"
#> 
#> [[98]]$description
#> [1] "thymopoietin pseudogene 1 [Source:HGNC Symbol;Acc:21364]"
#> 
#> [[98]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[98]]$feature_type
#> [1] "gene"
#> 
#> [[98]]$biotype
#> [1] "pseudogene"
#> 
#> [[98]]$start
#> [1] 30434229
#> 
#> [[98]]$id
#> [1] "ENSG00000229068"
#> 
#> [[98]]$end
#> [1] 30435771
#> 
#> [[98]]$external_name
#> [1] "TMPOP1"
#> 
#> [[98]]$canonical_transcript
#> [1] "ENST00000424642.1"
#> 
#> [[98]]$gene_id
#> [1] "ENSG00000229068"
#> 
#> [[98]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[98]]$version
#> [1] 1
#> 
#> [[98]]$strand
#> [1] 1
#> 
#> 
#> [[99]]
#> [[99]]$source
#> [1] "havana"
#> 
#> [[99]]$seq_region_name
#> [1] "6"
#> 
#> [[99]]$description
#> [1] "succinate-CoA ligase, ADP-forming, beta subunit pseudogene 1 [Source:HGNC Symbol;Acc:21632]"
#> 
#> [[99]]$biotype
#> [1] "pseudogene"
#> 
#> [[99]]$feature_type
#> [1] "gene"
#> 
#> [[99]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[99]]$start
#> [1] 30436659
#> 
#> [[99]]$canonical_transcript
#> [1] "ENST00000425839.1"
#> 
#> [[99]]$external_name
#> [1] "SUCLA2P1"
#> 
#> [[99]]$id
#> [1] "ENSG00000224936"
#> 
#> [[99]]$end
#> [1] 30438028
#> 
#> [[99]]$version
#> [1] 1
#> 
#> [[99]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[99]]$strand
#> [1] 1
#> 
#> [[99]]$gene_id
#> [1] "ENSG00000224936"
#> 
#> 
#> [[100]]
#> [[100]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[100]]$feature_type
#> [1] "gene"
#> 
#> [[100]]$biotype
#> [1] "pseudogene"
#> 
#> [[100]]$description
#> [1] "RAN, member RAS oncogene family pseudogene 1 [Source:HGNC Symbol;Acc:21631]"
#> 
#> [[100]]$seq_region_name
#> [1] "6"
#> 
#> [[100]]$source
#> [1] "ensembl_havana"
#> 
#> [[100]]$gene_id
#> [1] "ENSG00000236603"
#> 
#> [[100]]$strand
#> [1] 1
#> 
#> [[100]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[100]]$version
#> [1] 1
#> 
#> [[100]]$end
#> [1] 30454367
#> 
#> [[100]]$id
#> [1] "ENSG00000236603"
#> 
#> [[100]]$canonical_transcript
#> [1] "ENST00000455094.1"
#> 
#> [[100]]$external_name
#> [1] "RANP1"
#> 
#> [[100]]$start
#> [1] 30453717
#> 
#> 
#> [[101]]
#> [[101]]$start
#> [1] 30484043
#> 
#> [[101]]$end
#> [1] 30486994
#> 
#> [[101]]$id
#> [1] "ENSG00000235781"
#> 
#> [[101]]$external_name
#> [1] "XXbac-BPG249D20.9"
#> 
#> [[101]]$canonical_transcript
#> [1] "ENST00000415195.1"
#> 
#> [[101]]$gene_id
#> [1] "ENSG00000235781"
#> 
#> [[101]]$version
#> [1] 1
#> 
#> [[101]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[101]]$strand
#> [1] 1
#> 
#> [[101]]$seq_region_name
#> [1] "6"
#> 
#> [[101]]$source
#> [1] "havana"
#> 
#> [[101]]$description
#> NULL
#> 
#> [[101]]$feature_type
#> [1] "gene"
#> 
#> [[101]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[101]]$biotype
#> [1] "antisense"
#> 
#> 
#> [[102]]
#> [[102]]$start
#> [1] 30552109
#> 
#> [[102]]$end
#> [1] 30552194
#> 
#> [[102]]$id
#> [1] "ENSG00000216101"
#> 
#> [[102]]$external_name
#> [1] "MIR877"
#> 
#> [[102]]$canonical_transcript
#> [1] "ENST00000401282.1"
#> 
#> [[102]]$gene_id
#> [1] "ENSG00000216101"
#> 
#> [[102]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[102]]$version
#> [1] 1
#> 
#> [[102]]$strand
#> [1] 1
#> 
#> [[102]]$seq_region_name
#> [1] "6"
#> 
#> [[102]]$source
#> [1] "ensembl"
#> 
#> [[102]]$description
#> [1] "microRNA 877 [Source:HGNC Symbol;Acc:33660]"
#> 
#> [[102]]$feature_type
#> [1] "gene"
#> 
#> [[102]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[102]]$biotype
#> [1] "miRNA"
#> 
#> 
#> [[103]]
#> [[103]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[103]]$strand
#> [1] -1
#> 
#> [[103]]$version
#> [1] 1
#> 
#> [[103]]$gene_id
#> [1] "ENSG00000222894"
#> 
#> [[103]]$start
#> [1] 30584006
#> 
#> [[103]]$canonical_transcript
#> [1] "ENST00000410962.1"
#> 
#> [[103]]$external_name
#> [1] "AL662800.1"
#> 
#> [[103]]$end
#> [1] 30584082
#> 
#> [[103]]$id
#> [1] "ENSG00000222894"
#> 
#> [[103]]$description
#> NULL
#> 
#> [[103]]$biotype
#> [1] "miRNA"
#> 
#> [[103]]$feature_type
#> [1] "gene"
#> 
#> [[103]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[103]]$source
#> [1] "ensembl"
#> 
#> [[103]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[104]]
#> [[104]]$biotype
#> [1] "pseudogene"
#> 
#> [[104]]$feature_type
#> [1] "gene"
#> 
#> [[104]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[104]]$description
#> [1] "prothymosin, alpha pseudogene 1 (gene sequence 26) [Source:HGNC Symbol;Acc:9624]"
#> 
#> [[104]]$source
#> [1] "havana"
#> 
#> [[104]]$seq_region_name
#> [1] "6"
#> 
#> [[104]]$version
#> [1] 2
#> 
#> [[104]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[104]]$strand
#> [1] 1
#> 
#> [[104]]$gene_id
#> [1] "ENSG00000228415"
#> 
#> [[104]]$canonical_transcript
#> [1] "ENST00000455552.1"
#> 
#> [[104]]$external_name
#> [1] "PTMAP1"
#> 
#> [[104]]$id
#> [1] "ENSG00000228415"
#> 
#> [[104]]$end
#> [1] 30601675
#> 
#> [[104]]$start
#> [1] 30601409
#> 
#> 
#> [[105]]
#> [[105]]$start
#> [1] 30616443
#> 
#> [[105]]$canonical_transcript
#> [1] "ENST00000583820.1"
#> 
#> [[105]]$external_name
#> [1] "AL662800.2"
#> 
#> [[105]]$id
#> [1] "ENSG00000266183"
#> 
#> [[105]]$end
#> [1] 30616536
#> 
#> [[105]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[105]]$strand
#> [1] 1
#> 
#> [[105]]$version
#> [1] 1
#> 
#> [[105]]$gene_id
#> [1] "ENSG00000266183"
#> 
#> [[105]]$source
#> [1] "ensembl"
#> 
#> [[105]]$seq_region_name
#> [1] "6"
#> 
#> [[105]]$description
#> NULL
#> 
#> [[105]]$biotype
#> [1] "miRNA"
#> 
#> [[105]]$feature_type
#> [1] "gene"
#> 
#> [[105]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> 
#> [[106]]
#> [[106]]$seq_region_name
#> [1] "6"
#> 
#> [[106]]$source
#> [1] "havana"
#> 
#> [[106]]$description
#> [1] "ribosomal protein L7 pseudogene 4 [Source:HGNC Symbol;Acc:21634]"
#> 
#> [[106]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[106]]$feature_type
#> [1] "gene"
#> 
#> [[106]]$biotype
#> [1] "pseudogene"
#> 
#> [[106]]$start
#> [1] 30664583
#> 
#> [[106]]$end
#> [1] 30665312
#> 
#> [[106]]$id
#> [1] "ENSG00000230449"
#> 
#> [[106]]$canonical_transcript
#> [1] "ENST00000430239.1"
#> 
#> [[106]]$external_name
#> [1] "RPL7P4"
#> 
#> [[106]]$gene_id
#> [1] "ENSG00000230449"
#> 
#> [[106]]$version
#> [1] 1
#> 
#> [[106]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[106]]$strand
#> [1] -1
#> 
#> 
#> [[107]]
#> [[107]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[107]]$strand
#> [1] 1
#> 
#> [[107]]$version
#> [1] 1
#> 
#> [[107]]$gene_id
#> [1] "ENSG00000224328"
#> 
#> [[107]]$canonical_transcript
#> [1] "ENST00000442150.1"
#> 
#> [[107]]$external_name
#> [1] "MDC1-AS1"
#> 
#> [[107]]$id
#> [1] "ENSG00000224328"
#> 
#> [[107]]$end
#> [1] 30680961
#> 
#> [[107]]$start
#> [1] 30670844
#> 
#> [[107]]$biotype
#> [1] "antisense"
#> 
#> [[107]]$feature_type
#> [1] "gene"
#> 
#> [[107]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[107]]$description
#> [1] "MDC1 antisense RNA 1 [Source:HGNC Symbol;Acc:39764]"
#> 
#> [[107]]$source
#> [1] "havana"
#> 
#> [[107]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[108]]
#> [[108]]$seq_region_name
#> [1] "6"
#> 
#> [[108]]$source
#> [1] "havana"
#> 
#> [[108]]$description
#> NULL
#> 
#> [[108]]$feature_type
#> [1] "gene"
#> 
#> [[108]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[108]]$biotype
#> [1] "antisense"
#> 
#> [[108]]$start
#> [1] 30690882
#> 
#> [[108]]$end
#> [1] 30691654
#> 
#> [[108]]$id
#> [1] "ENSG00000272540"
#> 
#> [[108]]$external_name
#> [1] "XXbac-BPG252P9.9"
#> 
#> [[108]]$canonical_transcript
#> [1] "ENST00000607476.1"
#> 
#> [[108]]$gene_id
#> [1] "ENSG00000272540"
#> 
#> [[108]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[108]]$strand
#> [1] -1
#> 
#> [[108]]$version
#> [1] 1
#> 
#> 
#> [[109]]
#> [[109]]$start
#> [1] 30704081
#> 
#> [[109]]$id
#> [1] "ENSG00000201988"
#> 
#> [[109]]$end
#> [1] 30704191
#> 
#> [[109]]$canonical_transcript
#> [1] "ENST00000365118.2"
#> 
#> [[109]]$external_name
#> [1] "Y_RNA"
#> 
#> [[109]]$gene_id
#> [1] "ENSG00000201988"
#> 
#> [[109]]$strand
#> [1] -1
#> 
#> [[109]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[109]]$version
#> [1] 2
#> 
#> [[109]]$seq_region_name
#> [1] "6"
#> 
#> [[109]]$source
#> [1] "ensembl"
#> 
#> [[109]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[109]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[109]]$feature_type
#> [1] "gene"
#> 
#> [[109]]$biotype
#> [1] "misc_RNA"
#> 
#> 
#> [[110]]
#> [[110]]$seq_region_name
#> [1] "6"
#> 
#> [[110]]$source
#> [1] "havana"
#> 
#> [[110]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[110]]$feature_type
#> [1] "gene"
#> 
#> [[110]]$biotype
#> [1] "antisense"
#> 
#> [[110]]$description
#> NULL
#> 
#> [[110]]$end
#> [1] 30711369
#> 
#> [[110]]$id
#> [1] "ENSG00000272273"
#> 
#> [[110]]$external_name
#> [1] "XXbac-BPG252P9.10"
#> 
#> [[110]]$canonical_transcript
#> [1] "ENST00000607333.1"
#> 
#> [[110]]$start
#> [1] 30710706
#> 
#> [[110]]$gene_id
#> [1] "ENSG00000272273"
#> 
#> [[110]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[110]]$strand
#> [1] 1
#> 
#> [[110]]$version
#> [1] 1
#> 
#> 
#> [[111]]
#> [[111]]$description
#> [1] "RNA, 7SL, cytoplasmic 353, pseudogene [Source:HGNC Symbol;Acc:46369]"
#> 
#> [[111]]$biotype
#> [1] "misc_RNA"
#> 
#> [[111]]$feature_type
#> [1] "gene"
#> 
#> [[111]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[111]]$source
#> [1] "ensembl"
#> 
#> [[111]]$seq_region_name
#> [1] "6"
#> 
#> [[111]]$strand
#> [1] 1
#> 
#> [[111]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[111]]$version
#> [1] 1
#> 
#> [[111]]$gene_id
#> [1] "ENSG00000263608"
#> 
#> [[111]]$start
#> [1] 30718815
#> 
#> [[111]]$external_name
#> [1] "RN7SL353P"
#> 
#> [[111]]$canonical_transcript
#> [1] "ENST00000579902.1"
#> 
#> [[111]]$end
#> [1] 30719062
#> 
#> [[111]]$id
#> [1] "ENSG00000263608"
#> 
#> 
#> [[112]]
#> [[112]]$description
#> [1] "HLA complex group 20 (non-protein coding) [Source:HGNC Symbol;Acc:31334]"
#> 
#> [[112]]$feature_type
#> [1] "gene"
#> 
#> [[112]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[112]]$biotype
#> [1] "lincRNA"
#> 
#> [[112]]$seq_region_name
#> [1] "6"
#> 
#> [[112]]$source
#> [1] "havana"
#> 
#> [[112]]$gene_id
#> [1] "ENSG00000228022"
#> 
#> [[112]]$version
#> [1] 1
#> 
#> [[112]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[112]]$strand
#> [1] 1
#> 
#> [[112]]$start
#> [1] 30734602
#> 
#> [[112]]$end
#> [1] 30760027
#> 
#> [[112]]$id
#> [1] "ENSG00000228022"
#> 
#> [[112]]$canonical_transcript
#> [1] "ENST00000439406.1"
#> 
#> [[112]]$external_name
#> [1] "HCG20"
#> 
#> 
#> [[113]]
#> [[113]]$description
#> [1] "long intergenic non-protein coding RNA 243 [Source:HGNC Symbol;Acc:30956]"
#> 
#> [[113]]$feature_type
#> [1] "gene"
#> 
#> [[113]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[113]]$biotype
#> [1] "processed_transcript"
#> 
#> [[113]]$seq_region_name
#> [1] "6"
#> 
#> [[113]]$source
#> [1] "havana"
#> 
#> [[113]]$gene_id
#> [1] "ENSG00000214894"
#> 
#> [[113]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[113]]$version
#> [1] 2
#> 
#> [[113]]$strand
#> [1] -1
#> 
#> [[113]]$start
#> [1] 30766431
#> 
#> [[113]]$id
#> [1] "ENSG00000214894"
#> 
#> [[113]]$end
#> [1] 30798436
#> 
#> [[113]]$external_name
#> [1] "LINC00243"
#> 
#> [[113]]$canonical_transcript
#> [1] "ENST00000399196.1"
#> 
#> 
#> [[114]]
#> [[114]]$seq_region_name
#> [1] "6"
#> 
#> [[114]]$source
#> [1] "havana"
#> 
#> [[114]]$description
#> NULL
#> 
#> [[114]]$feature_type
#> [1] "gene"
#> 
#> [[114]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[114]]$biotype
#> [1] "lincRNA"
#> 
#> [[114]]$start
#> [1] 30807303
#> 
#> [[114]]$end
#> [1] 30815936
#> 
#> [[114]]$id
#> [1] "ENSG00000237923"
#> 
#> [[114]]$canonical_transcript
#> [1] "ENST00000442852.1"
#> 
#> [[114]]$external_name
#> [1] "XXbac-BPG27H4.8"
#> 
#> [[114]]$gene_id
#> [1] "ENSG00000237923"
#> 
#> [[114]]$strand
#> [1] -1
#> 
#> [[114]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[114]]$version
#> [1] 1
#> 
#> 
#> [[115]]
#> [[115]]$external_name
#> [1] "RN7SKP186"
#> 
#> [[115]]$canonical_transcript
#> [1] "ENST00000365371.1"
#> 
#> [[115]]$end
#> [1] 30832329
#> 
#> [[115]]$id
#> [1] "ENSG00000202241"
#> 
#> [[115]]$start
#> [1] 30832027
#> 
#> [[115]]$version
#> [1] 1
#> 
#> [[115]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[115]]$strand
#> [1] 1
#> 
#> [[115]]$gene_id
#> [1] "ENSG00000202241"
#> 
#> [[115]]$source
#> [1] "ensembl"
#> 
#> [[115]]$seq_region_name
#> [1] "6"
#> 
#> [[115]]$biotype
#> [1] "misc_RNA"
#> 
#> [[115]]$feature_type
#> [1] "gene"
#> 
#> [[115]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[115]]$description
#> [1] "RNA, 7SK small nuclear pseudogene 186 [Source:HGNC Symbol;Acc:45910]"
#> 
#> 
#> [[116]]
#> [[116]]$version
#> [1] 1
#> 
#> [[116]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[116]]$strand
#> [1] -1
#> 
#> [[116]]$gene_id
#> [1] "ENSG00000237775"
#> 
#> [[116]]$external_name
#> [1] "DDR1-AS1"
#> 
#> [[116]]$canonical_transcript
#> [1] "ENST00000458361.1"
#> 
#> [[116]]$end
#> [1] 30843695
#> 
#> [[116]]$id
#> [1] "ENSG00000237775"
#> 
#> [[116]]$start
#> [1] 30834759
#> 
#> [[116]]$biotype
#> [1] "antisense"
#> 
#> [[116]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[116]]$feature_type
#> [1] "gene"
#> 
#> [[116]]$description
#> [1] "DDR1 antisense RNA 1 (head to head) [Source:HGNC Symbol;Acc:28694]"
#> 
#> [[116]]$source
#> [1] "havana"
#> 
#> [[116]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[117]]
#> [[117]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[117]]$feature_type
#> [1] "gene"
#> 
#> [[117]]$biotype
#> [1] "miRNA"
#> 
#> [[117]]$description
#> [1] "microRNA 4640 [Source:HGNC Symbol;Acc:41561]"
#> 
#> [[117]]$seq_region_name
#> [1] "6"
#> 
#> [[117]]$source
#> [1] "ensembl"
#> 
#> [[117]]$gene_id
#> [1] "ENSG00000264594"
#> 
#> [[117]]$version
#> [1] 1
#> 
#> [[117]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[117]]$strand
#> [1] 1
#> 
#> [[117]]$end
#> [1] 30858749
#> 
#> [[117]]$id
#> [1] "ENSG00000264594"
#> 
#> [[117]]$external_name
#> [1] "MIR4640"
#> 
#> [[117]]$canonical_transcript
#> [1] "ENST00000581824.1"
#> 
#> [[117]]$start
#> [1] 30858660
#> 
#> 
#> [[118]]
#> [[118]]$version
#> [1] 1
#> 
#> [[118]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[118]]$strand
#> [1] 1
#> 
#> [[118]]$gene_id
#> [1] "ENSG00000264731"
#> 
#> [[118]]$start
#> [1] 30874655
#> 
#> [[118]]$canonical_transcript
#> [1] "ENST00000580375.1"
#> 
#> [[118]]$external_name
#> [1] "RN7SL175P"
#> 
#> [[118]]$end
#> [1] 30874913
#> 
#> [[118]]$id
#> [1] "ENSG00000264731"
#> 
#> [[118]]$description
#> [1] "RNA, 7SL, cytoplasmic 175, pseudogene [Source:HGNC Symbol;Acc:46191]"
#> 
#> [[118]]$biotype
#> [1] "misc_RNA"
#> 
#> [[118]]$feature_type
#> [1] "gene"
#> 
#> [[118]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[118]]$source
#> [1] "ensembl"
#> 
#> [[118]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[119]]
#> [[119]]$gene_id
#> [1] "ENSG00000252761"
#> 
#> [[119]]$version
#> [1] 1
#> 
#> [[119]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[119]]$strand
#> [1] -1
#> 
#> [[119]]$start
#> [1] 30900395
#> 
#> [[119]]$end
#> [1] 30900486
#> 
#> [[119]]$id
#> [1] "ENSG00000252761"
#> 
#> [[119]]$external_name
#> [1] "Y_RNA"
#> 
#> [[119]]$canonical_transcript
#> [1] "ENST00000516952.1"
#> 
#> [[119]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[119]]$feature_type
#> [1] "gene"
#> 
#> [[119]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[119]]$biotype
#> [1] "misc_RNA"
#> 
#> [[119]]$seq_region_name
#> [1] "6"
#> 
#> [[119]]$source
#> [1] "ensembl"
#> 
#> 
#> [[120]]
#> [[120]]$id
#> [1] "ENSG00000233529"
#> 
#> [[120]]$end
#> [1] 30922639
#> 
#> [[120]]$external_name
#> [1] "HCG21"
#> 
#> [[120]]$canonical_transcript
#> [1] "ENST00000419481.1"
#> 
#> [[120]]$start
#> [1] 30913756
#> 
#> [[120]]$gene_id
#> [1] "ENSG00000233529"
#> 
#> [[120]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[120]]$strand
#> [1] -1
#> 
#> [[120]]$version
#> [1] 1
#> 
#> [[120]]$seq_region_name
#> [1] "6"
#> 
#> [[120]]$source
#> [1] "havana"
#> 
#> [[120]]$feature_type
#> [1] "gene"
#> 
#> [[120]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[120]]$biotype
#> [1] "antisense"
#> 
#> [[120]]$description
#> [1] "HLA complex group 21 (non-protein coding) [Source:HGNC Symbol;Acc:31335]"
#> 
#> 
#> [[121]]
#> [[121]]$start
#> [1] 31021227
#> 
#> [[121]]$canonical_transcript
#> [1] "ENST00000426185.1"
#> 
#> [[121]]$external_name
#> [1] "HCG22"
#> 
#> [[121]]$id
#> [1] "ENSG00000228789"
#> 
#> [[121]]$end
#> [1] 31027667
#> 
#> [[121]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[121]]$strand
#> [1] 1
#> 
#> [[121]]$version
#> [1] 2
#> 
#> [[121]]$gene_id
#> [1] "ENSG00000228789"
#> 
#> [[121]]$source
#> [1] "havana"
#> 
#> [[121]]$seq_region_name
#> [1] "6"
#> 
#> [[121]]$description
#> [1] "HLA complex group 22 [Source:HGNC Symbol;Acc:27780]"
#> 
#> [[121]]$biotype
#> [1] "lincRNA"
#> 
#> [[121]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[121]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[122]]
#> [[122]]$strand
#> [1] 1
#> 
#> [[122]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[122]]$version
#> [1] 1
#> 
#> [[122]]$gene_id
#> [1] "ENSG00000222895"
#> 
#> [[122]]$start
#> [1] 31050787
#> 
#> [[122]]$external_name
#> [1] "RNU6-1133P"
#> 
#> [[122]]$canonical_transcript
#> [1] "ENST00000410963.1"
#> 
#> [[122]]$end
#> [1] 31050886
#> 
#> [[122]]$id
#> [1] "ENSG00000222895"
#> 
#> [[122]]$description
#> [1] "RNA, U6 small nuclear 1133, pseudogene [Source:HGNC Symbol;Acc:48096]"
#> 
#> [[122]]$biotype
#> [1] "snRNA"
#> 
#> [[122]]$feature_type
#> [1] "gene"
#> 
#> [[122]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[122]]$source
#> [1] "ensembl"
#> 
#> [[122]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[123]]
#> [[123]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[123]]$strand
#> [1] 1
#> 
#> [[123]]$version
#> [1] 1
#> 
#> [[123]]$gene_id
#> [1] "ENSG00000238211"
#> 
#> [[123]]$external_name
#> [1] "POLR2LP"
#> 
#> [[123]]$canonical_transcript
#> [1] "ENST00000444785.1"
#> 
#> [[123]]$end
#> [1] 31108690
#> 
#> [[123]]$id
#> [1] "ENSG00000238211"
#> 
#> [[123]]$start
#> [1] 31108504
#> 
#> [[123]]$biotype
#> [1] "pseudogene"
#> 
#> [[123]]$feature_type
#> [1] "gene"
#> 
#> [[123]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[123]]$description
#> [1] "polymerase (RNA) II (DNA directed) polypeptide L pseudogene [Source:HGNC Symbol;Acc:31340]"
#> 
#> [[123]]$source
#> [1] "havana"
#> 
#> [[123]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[124]]
#> [[124]]$description
#> [1] "psoriasis susceptibility 1 candidate 3 (non-protein coding) [Source:HGNC Symbol;Acc:17203]"
#> 
#> [[124]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[124]]$feature_type
#> [1] "gene"
#> 
#> [[124]]$biotype
#> [1] "sense_intronic"
#> 
#> [[124]]$seq_region_name
#> [1] "6"
#> 
#> [[124]]$source
#> [1] "ensembl_havana"
#> 
#> [[124]]$gene_id
#> [1] "ENSG00000204528"
#> 
#> [[124]]$strand
#> [1] -1
#> 
#> [[124]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[124]]$version
#> [1] 3
#> 
#> [[124]]$start
#> [1] 31141512
#> 
#> [[124]]$end
#> [1] 31145676
#> 
#> [[124]]$id
#> [1] "ENSG00000204528"
#> 
#> [[124]]$external_name
#> [1] "PSORS1C3"
#> 
#> [[124]]$canonical_transcript
#> [1] "ENST00000412143.1"
#> 
#> 
#> [[125]]
#> [[125]]$canonical_transcript
#> [1] "ENST00000606367.1"
#> 
#> [[125]]$external_name
#> [1] "XXbac-BPG299F13.17"
#> 
#> [[125]]$id
#> [1] "ENSG00000272501"
#> 
#> [[125]]$end
#> [1] 31165814
#> 
#> [[125]]$start
#> [1] 31162977
#> 
#> [[125]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[125]]$version
#> [1] 1
#> 
#> [[125]]$strand
#> [1] -1
#> 
#> [[125]]$gene_id
#> [1] "ENSG00000272501"
#> 
#> [[125]]$source
#> [1] "havana"
#> 
#> [[125]]$seq_region_name
#> [1] "6"
#> 
#> [[125]]$biotype
#> [1] "antisense"
#> 
#> [[125]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[125]]$feature_type
#> [1] "gene"
#> 
#> [[125]]$description
#> NULL
#> 
#> 
#> [[126]]
#> [[126]]$description
#> NULL
#> 
#> [[126]]$biotype
#> [1] "lincRNA"
#> 
#> [[126]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[126]]$feature_type
#> [1] "gene"
#> 
#> [[126]]$source
#> [1] "havana"
#> 
#> [[126]]$seq_region_name
#> [1] "6"
#> 
#> [[126]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[126]]$strand
#> [1] -1
#> 
#> [[126]]$version
#> [1] 1
#> 
#> [[126]]$gene_id
#> [1] "ENSG00000271821"
#> 
#> [[126]]$start
#> [1] 31167942
#> 
#> [[126]]$external_name
#> [1] "XXbac-BPG299F13.14"
#> 
#> [[126]]$canonical_transcript
#> [1] "ENST00000606909.1"
#> 
#> [[126]]$id
#> [1] "ENSG00000271821"
#> 
#> [[126]]$end
#> [1] 31169695
#> 
#> 
#> [[127]]
#> [[127]]$strand
#> [1] 1
#> 
#> [[127]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[127]]$version
#> [1] 1
#> 
#> [[127]]$gene_id
#> [1] "ENSG00000255726"
#> 
#> [[127]]$start
#> [1] 31190690
#> 
#> [[127]]$canonical_transcript
#> [1] "ENST00000546340.1"
#> 
#> [[127]]$external_name
#> [1] "XXbac-BPG299F13.15"
#> 
#> [[127]]$end
#> [1] 31190870
#> 
#> [[127]]$id
#> [1] "ENSG00000255726"
#> 
#> [[127]]$description
#> NULL
#> 
#> [[127]]$biotype
#> [1] "pseudogene"
#> 
#> [[127]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[127]]$feature_type
#> [1] "gene"
#> 
#> [[127]]$source
#> [1] "havana"
#> 
#> [[127]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[128]]
#> [[128]]$strand
#> [1] -1
#> 
#> [[128]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[128]]$version
#> [1] 1
#> 
#> [[128]]$gene_id
#> [1] "ENSG00000255899"
#> 
#> [[128]]$start
#> [1] 31192119
#> 
#> [[128]]$canonical_transcript
#> [1] "ENST00000538455.1"
#> 
#> [[128]]$external_name
#> [1] "XXbac-BPG299F13.16"
#> 
#> [[128]]$id
#> [1] "ENSG00000255899"
#> 
#> [[128]]$end
#> [1] 31192835
#> 
#> [[128]]$description
#> NULL
#> 
#> [[128]]$biotype
#> [1] "pseudogene"
#> 
#> [[128]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[128]]$feature_type
#> [1] "gene"
#> 
#> [[128]]$source
#> [1] "havana"
#> 
#> [[128]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[129]]
#> [[129]]$seq_region_name
#> [1] "6"
#> 
#> [[129]]$source
#> [1] "havana"
#> 
#> [[129]]$description
#> [1] "ubiquitin specific peptidase 8 pseudogene 1 [Source:HGNC Symbol;Acc:13987]"
#> 
#> [[129]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[129]]$feature_type
#> [1] "gene"
#> 
#> [[129]]$biotype
#> [1] "pseudogene"
#> 
#> [[129]]$start
#> [1] 31243349
#> 
#> [[129]]$end
#> [1] 31246531
#> 
#> [[129]]$id
#> [1] "ENSG00000214892"
#> 
#> [[129]]$external_name
#> [1] "USP8P1"
#> 
#> [[129]]$canonical_transcript
#> [1] "ENST00000494673.1"
#> 
#> [[129]]$gene_id
#> [1] "ENSG00000214892"
#> 
#> [[129]]$strand
#> [1] 1
#> 
#> [[129]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[129]]$version
#> [1] 4
#> 
#> 
#> [[130]]
#> [[130]]$id
#> [1] "ENSG00000227939"
#> 
#> [[130]]$end
#> [1] 31249296
#> 
#> [[130]]$external_name
#> [1] "RPL3P2"
#> 
#> [[130]]$canonical_transcript
#> [1] "ENST00000413027.1"
#> 
#> [[130]]$start
#> [1] 31248094
#> 
#> [[130]]$gene_id
#> [1] "ENSG00000227939"
#> 
#> [[130]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[130]]$version
#> [1] 1
#> 
#> [[130]]$strand
#> [1] 1
#> 
#> [[130]]$seq_region_name
#> [1] "6"
#> 
#> [[130]]$source
#> [1] "havana"
#> 
#> [[130]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[130]]$feature_type
#> [1] "gene"
#> 
#> [[130]]$biotype
#> [1] "pseudogene"
#> 
#> [[130]]$description
#> [1] "ribosomal protein L3 pseudogene 2 [Source:HGNC Symbol;Acc:17834]"
#> 
#> 
#> [[131]]
#> [[131]]$source
#> [1] "havana"
#> 
#> [[131]]$seq_region_name
#> [1] "6"
#> 
#> [[131]]$description
#> [1] "WAS protein family, member 5, pseudogene [Source:HGNC Symbol;Acc:21665]"
#> 
#> [[131]]$biotype
#> [1] "pseudogene"
#> 
#> [[131]]$feature_type
#> [1] "gene"
#> 
#> [[131]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[131]]$start
#> [1] 31255287
#> 
#> [[131]]$canonical_transcript
#> [1] "ENST00000428639.1"
#> 
#> [[131]]$external_name
#> [1] "WASF5P"
#> 
#> [[131]]$id
#> [1] "ENSG00000231402"
#> 
#> [[131]]$end
#> [1] 31256741
#> 
#> [[131]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[131]]$version
#> [1] 1
#> 
#> [[131]]$strand
#> [1] -1
#> 
#> [[131]]$gene_id
#> [1] "ENSG00000231402"
#> 
#> 
#> [[132]]
#> [[132]]$start
#> [1] 31261685
#> 
#> [[132]]$id
#> [1] "ENSG00000256166"
#> 
#> [[132]]$end
#> [1] 31269419
#> 
#> [[132]]$external_name
#> [1] "XXbac-BPG248L24.13"
#> 
#> [[132]]$canonical_transcript
#> [1] "ENST00000539514.1"
#> 
#> [[132]]$gene_id
#> [1] "ENSG00000256166"
#> 
#> [[132]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[132]]$strand
#> [1] -1
#> 
#> [[132]]$version
#> [1] 1
#> 
#> [[132]]$seq_region_name
#> [1] "6"
#> 
#> [[132]]$source
#> [1] "havana"
#> 
#> [[132]]$description
#> NULL
#> 
#> [[132]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[132]]$feature_type
#> [1] "gene"
#> 
#> [[132]]$biotype
#> [1] "lincRNA"
#> 
#> 
#> [[133]]
#> [[133]]$source
#> [1] "havana"
#> 
#> [[133]]$seq_region_name
#> [1] "6"
#> 
#> [[133]]$description
#> NULL
#> 
#> [[133]]$biotype
#> [1] "pseudogene"
#> 
#> [[133]]$feature_type
#> [1] "gene"
#> 
#> [[133]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[133]]$start
#> [1] 31275592
#> 
#> [[133]]$canonical_transcript
#> [1] "ENST00000421191.1"
#> 
#> [[133]]$external_name
#> [1] "XXbac-BPG248L24.10"
#> 
#> [[133]]$end
#> [1] 31276326
#> 
#> [[133]]$id
#> [1] "ENSG00000229836"
#> 
#> [[133]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[133]]$strand
#> [1] -1
#> 
#> [[133]]$version
#> [1] 1
#> 
#> [[133]]$gene_id
#> [1] "ENSG00000229836"
#> 
#> 
#> [[134]]
#> [[134]]$external_name
#> [1] "XXbac-BPG248L24.12"
#> 
#> [[134]]$canonical_transcript
#> [1] "ENST00000603274.1"
#> 
#> [[134]]$end
#> [1] 31325414
#> 
#> [[134]]$id
#> [1] "ENSG00000271581"
#> 
#> [[134]]$start
#> [1] 31324424
#> 
#> [[134]]$version
#> [1] 1
#> 
#> [[134]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[134]]$strand
#> [1] 1
#> 
#> [[134]]$gene_id
#> [1] "ENSG00000271581"
#> 
#> [[134]]$source
#> [1] "havana"
#> 
#> [[134]]$seq_region_name
#> [1] "6"
#> 
#> [[134]]$biotype
#> [1] "pseudogene"
#> 
#> [[134]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[134]]$feature_type
#> [1] "gene"
#> 
#> [[134]]$description
#> NULL
#> 
#> 
#> [[135]]
#> [[135]]$gene_id
#> [1] "ENSG00000228432"
#> 
#> [[135]]$version
#> [1] 1
#> 
#> [[135]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[135]]$strand
#> [1] -1
#> 
#> [[135]]$start
#> [1] 31334129
#> 
#> [[135]]$end
#> [1] 31334675
#> 
#> [[135]]$id
#> [1] "ENSG00000228432"
#> 
#> [[135]]$external_name
#> [1] "DHFRP2"
#> 
#> [[135]]$canonical_transcript
#> [1] "ENST00000414224.1"
#> 
#> [[135]]$description
#> [1] "dihydrofolate reductase pseudogene 2 [Source:HGNC Symbol;Acc:2863]"
#> 
#> [[135]]$feature_type
#> [1] "gene"
#> 
#> [[135]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[135]]$biotype
#> [1] "pseudogene"
#> 
#> [[135]]$seq_region_name
#> [1] "6"
#> 
#> [[135]]$source
#> [1] "havana"
#> 
#> 
#> [[136]]
#> [[136]]$source
#> [1] "ensembl"
#> 
#> [[136]]$seq_region_name
#> [1] "6"
#> 
#> [[136]]$biotype
#> [1] "snRNA"
#> 
#> [[136]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[136]]$feature_type
#> [1] "gene"
#> 
#> [[136]]$description
#> [1] "RNA, U6 small nuclear 283, pseudogene [Source:HGNC Symbol;Acc:47246]"
#> 
#> [[136]]$canonical_transcript
#> [1] "ENST00000364788.1"
#> 
#> [[136]]$external_name
#> [1] "RNU6-283P"
#> 
#> [[136]]$id
#> [1] "ENSG00000201658"
#> 
#> [[136]]$end
#> [1] 31338017
#> 
#> [[136]]$start
#> [1] 31337911
#> 
#> [[136]]$version
#> [1] 1
#> 
#> [[136]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[136]]$strand
#> [1] 1
#> 
#> [[136]]$gene_id
#> [1] "ENSG00000201658"
#> 
#> 
#> [[137]]
#> [[137]]$seq_region_name
#> [1] "6"
#> 
#> [[137]]$source
#> [1] "ensembl"
#> 
#> [[137]]$description
#> NULL
#> 
#> [[137]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[137]]$feature_type
#> [1] "gene"
#> 
#> [[137]]$biotype
#> [1] "miRNA"
#> 
#> [[137]]$start
#> [1] 31342095
#> 
#> [[137]]$id
#> [1] "ENSG00000265294"
#> 
#> [[137]]$end
#> [1] 31342180
#> 
#> [[137]]$canonical_transcript
#> [1] "ENST00000578513.1"
#> 
#> [[137]]$external_name
#> [1] "AL671883.1"
#> 
#> [[137]]$gene_id
#> [1] "ENSG00000265294"
#> 
#> [[137]]$strand
#> [1] -1
#> 
#> [[137]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[137]]$version
#> [1] 1
#> 
#> 
#> [[138]]
#> [[138]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[138]]$strand
#> [1] 1
#> 
#> [[138]]$version
#> [1] 1
#> 
#> [[138]]$gene_id
#> [1] "ENSG00000230994"
#> 
#> [[138]]$start
#> [1] 31345196
#> 
#> [[138]]$external_name
#> [1] "FGFR3P1"
#> 
#> [[138]]$canonical_transcript
#> [1] "ENST00000449999.1"
#> 
#> [[138]]$end
#> [1] 31345796
#> 
#> [[138]]$id
#> [1] "ENSG00000230994"
#> 
#> [[138]]$description
#> [1] "fibroblast growth factor receptor 3 pseudogene 1 [Source:HGNC Symbol;Acc:21664]"
#> 
#> [[138]]$biotype
#> [1] "pseudogene"
#> 
#> [[138]]$feature_type
#> [1] "gene"
#> 
#> [[138]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[138]]$source
#> [1] "havana"
#> 
#> [[138]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[139]]
#> [[139]]$canonical_transcript
#> [1] "ENST00000424108.1"
#> 
#> [[139]]$external_name
#> [1] "ZDHHC20P2"
#> 
#> [[139]]$end
#> [1] 31348616
#> 
#> [[139]]$id
#> [1] "ENSG00000223702"
#> 
#> [[139]]$start
#> [1] 31348188
#> 
#> [[139]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[139]]$version
#> [1] 1
#> 
#> [[139]]$strand
#> [1] 1
#> 
#> [[139]]$gene_id
#> [1] "ENSG00000223702"
#> 
#> [[139]]$source
#> [1] "havana"
#> 
#> [[139]]$seq_region_name
#> [1] "6"
#> 
#> [[139]]$biotype
#> [1] "pseudogene"
#> 
#> [[139]]$feature_type
#> [1] "gene"
#> 
#> [[139]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[139]]$description
#> [1] "zinc finger, DHHC-type containing 20 pseudogene 2 [Source:HGNC Symbol;Acc:33457]"
#> 
#> 
#> [[140]]
#> [[140]]$gene_id
#> [1] "ENSG00000225851"
#> 
#> [[140]]$strand
#> [1] -1
#> 
#> [[140]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[140]]$version
#> [1] 1
#> 
#> [[140]]$end
#> [1] 31350065
#> 
#> [[140]]$id
#> [1] "ENSG00000225851"
#> 
#> [[140]]$external_name
#> [1] "HLA-S"
#> 
#> [[140]]$canonical_transcript
#> [1] "ENST00000425174.1"
#> 
#> [[140]]$start
#> [1] 31349851
#> 
#> [[140]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[140]]$feature_type
#> [1] "gene"
#> 
#> [[140]]$biotype
#> [1] "pseudogene"
#> 
#> [[140]]$description
#> [1] "major histocompatibility complex, class I, S (pseudogene) [Source:HGNC Symbol;Acc:19395]"
#> 
#> [[140]]$seq_region_name
#> [1] "6"
#> 
#> [[140]]$source
#> [1] "havana"
#> 
#> 
#> [[141]]
#> [[141]]$source
#> [1] "havana"
#> 
#> [[141]]$seq_region_name
#> [1] "6"
#> 
#> [[141]]$biotype
#> [1] "lincRNA"
#> 
#> [[141]]$feature_type
#> [1] "gene"
#> 
#> [[141]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[141]]$description
#> NULL
#> 
#> [[141]]$external_name
#> [1] "XXbac-BPG181B23.7"
#> 
#> [[141]]$canonical_transcript
#> [1] "ENST00000606743.1"
#> 
#> [[141]]$id
#> [1] "ENSG00000272221"
#> 
#> [[141]]$end
#> [1] 31363272
#> 
#> [[141]]$start
#> [1] 31362066
#> 
#> [[141]]$strand
#> [1] -1
#> 
#> [[141]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[141]]$version
#> [1] 1
#> 
#> [[141]]$gene_id
#> [1] "ENSG00000272221"
#> 
#> 
#> [[142]]
#> [[142]]$source
#> [1] "ensembl_havana"
#> 
#> [[142]]$seq_region_name
#> [1] "6"
#> 
#> [[142]]$biotype
#> [1] "sense_overlapping"
#> 
#> [[142]]$feature_type
#> [1] "gene"
#> 
#> [[142]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[142]]$description
#> [1] "HLA complex P5 (non-protein coding) [Source:HGNC Symbol;Acc:21659]"
#> 
#> [[142]]$external_name
#> [1] "HCP5"
#> 
#> [[142]]$canonical_transcript
#> [1] "ENST00000414046.2"
#> 
#> [[142]]$id
#> [1] "ENSG00000206337"
#> 
#> [[142]]$end
#> [1] 31445283
#> 
#> [[142]]$start
#> [1] 31368479
#> 
#> [[142]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[142]]$strand
#> [1] 1
#> 
#> [[142]]$version
#> [1] 6
#> 
#> [[142]]$gene_id
#> [1] "ENSG00000206337"
#> 
#> 
#> [[143]]
#> [[143]]$gene_id
#> [1] "ENSG00000199332"
#> 
#> [[143]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[143]]$strand
#> [1] -1
#> 
#> [[143]]$version
#> [1] 1
#> 
#> [[143]]$start
#> [1] 31369929
#> 
#> [[143]]$end
#> [1] 31370027
#> 
#> [[143]]$id
#> [1] "ENSG00000199332"
#> 
#> [[143]]$external_name
#> [1] "Y_RNA"
#> 
#> [[143]]$canonical_transcript
#> [1] "ENST00000362462.1"
#> 
#> [[143]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[143]]$feature_type
#> [1] "gene"
#> 
#> [[143]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[143]]$biotype
#> [1] "misc_RNA"
#> 
#> [[143]]$seq_region_name
#> [1] "6"
#> 
#> [[143]]$source
#> [1] "ensembl"
#> 
#> 
#> [[144]]
#> [[144]]$id
#> [1] "ENSG00000230174"
#> 
#> [[144]]$end
#> [1] 31414750
#> 
#> [[144]]$canonical_transcript
#> [1] "ENST00000430364.1"
#> 
#> [[144]]$external_name
#> [1] "LINC01149"
#> 
#> [[144]]$start
#> [1] 31409444
#> 
#> [[144]]$gene_id
#> [1] "ENSG00000230174"
#> 
#> [[144]]$version
#> [1] 1
#> 
#> [[144]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[144]]$strand
#> [1] 1
#> 
#> [[144]]$seq_region_name
#> [1] "6"
#> 
#> [[144]]$source
#> [1] "havana"
#> 
#> [[144]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[144]]$feature_type
#> [1] "gene"
#> 
#> [[144]]$biotype
#> [1] "lincRNA"
#> 
#> [[144]]$description
#> [1] "long intergenic non-protein coding RNA 1149 [Source:HGNC Symbol;Acc:39757]"
#> 
#> 
#> [[145]]
#> [[145]]$strand
#> [1] -1
#> 
#> [[145]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[145]]$version
#> [1] 1
#> 
#> [[145]]$gene_id
#> [1] "ENSG00000233902"
#> 
#> [[145]]$start
#> [1] 31430505
#> 
#> [[145]]$canonical_transcript
#> [1] "ENST00000440087.1"
#> 
#> [[145]]$external_name
#> [1] "XXbac-BPG181B23.6"
#> 
#> [[145]]$end
#> [1] 31431113
#> 
#> [[145]]$id
#> [1] "ENSG00000233902"
#> 
#> [[145]]$description
#> NULL
#> 
#> [[145]]$biotype
#> [1] "pseudogene"
#> 
#> [[145]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[145]]$feature_type
#> [1] "gene"
#> 
#> [[145]]$source
#> [1] "havana"
#> 
#> [[145]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[146]]
#> [[146]]$start
#> [1] 31464466
#> 
#> [[146]]$end
#> [1] 31464567
#> 
#> [[146]]$id
#> [1] "ENSG00000201680"
#> 
#> [[146]]$external_name
#> [1] "Y_RNA"
#> 
#> [[146]]$canonical_transcript
#> [1] "ENST00000383850.1"
#> 
#> [[146]]$gene_id
#> [1] "ENSG00000201680"
#> 
#> [[146]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[146]]$version
#> [1] 1
#> 
#> [[146]]$strand
#> [1] -1
#> 
#> [[146]]$seq_region_name
#> [1] "6"
#> 
#> [[146]]$source
#> [1] "ensembl"
#> 
#> [[146]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[146]]$feature_type
#> [1] "gene"
#> 
#> [[146]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[146]]$biotype
#> [1] "misc_RNA"
#> 
#> 
#> [[147]]
#> [[147]]$start
#> [1] 31483756
#> 
#> [[147]]$end
#> [1] 31483988
#> 
#> [[147]]$id
#> [1] "ENSG00000256851"
#> 
#> [[147]]$external_name
#> [1] "XXbac-BPG16N22.5"
#> 
#> [[147]]$canonical_transcript
#> [1] "ENST00000538358.1"
#> 
#> [[147]]$gene_id
#> [1] "ENSG00000256851"
#> 
#> [[147]]$version
#> [1] 1
#> 
#> [[147]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[147]]$strand
#> [1] 1
#> 
#> [[147]]$seq_region_name
#> [1] "6"
#> 
#> [[147]]$source
#> [1] "havana"
#> 
#> [[147]]$description
#> NULL
#> 
#> [[147]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[147]]$feature_type
#> [1] "gene"
#> 
#> [[147]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[148]]
#> [[148]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[148]]$version
#> [1] 2
#> 
#> [[148]]$strand
#> [1] -1
#> 
#> [[148]]$gene_id
#> [1] "ENSG00000219797"
#> 
#> [[148]]$start
#> [1] 31487257
#> 
#> [[148]]$canonical_transcript
#> [1] "ENST00000403866.2"
#> 
#> [[148]]$external_name
#> [1] "PPIAP9"
#> 
#> [[148]]$id
#> [1] "ENSG00000219797"
#> 
#> [[148]]$end
#> [1] 31488068
#> 
#> [[148]]$description
#> [1] "peptidylprolyl isomerase A (cyclophilin A) pseudogene 9 [Source:HGNC Symbol;Acc:9272]"
#> 
#> [[148]]$biotype
#> [1] "pseudogene"
#> 
#> [[148]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[148]]$feature_type
#> [1] "gene"
#> 
#> [[148]]$source
#> [1] "havana"
#> 
#> [[148]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[149]]
#> [[149]]$biotype
#> [1] "pseudogene"
#> 
#> [[149]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[149]]$feature_type
#> [1] "gene"
#> 
#> [[149]]$description
#> [1] "ribosomal protein L15 pseudogene 4 [Source:HGNC Symbol;Acc:21663]"
#> 
#> [[149]]$source
#> [1] "havana"
#> 
#> [[149]]$seq_region_name
#> [1] "6"
#> 
#> [[149]]$strand
#> [1] 1
#> 
#> [[149]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[149]]$version
#> [1] 1
#> 
#> [[149]]$gene_id
#> [1] "ENSG00000225499"
#> 
#> [[149]]$canonical_transcript
#> [1] "ENST00000416625.1"
#> 
#> [[149]]$external_name
#> [1] "RPL15P4"
#> 
#> [[149]]$id
#> [1] "ENSG00000225499"
#> 
#> [[149]]$end
#> [1] 31496470
#> 
#> [[149]]$start
#> [1] 31495891
#> 
#> 
#> [[150]]
#> [[150]]$gene_id
#> [1] "ENSG00000201785"
#> 
#> [[150]]$strand
#> [1] -1
#> 
#> [[150]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[150]]$version
#> [1] 1
#> 
#> [[150]]$start
#> [1] 31504151
#> 
#> [[150]]$end
#> [1] 31504226
#> 
#> [[150]]$id
#> [1] "ENSG00000201785"
#> 
#> [[150]]$canonical_transcript
#> [1] "ENST00000364915.1"
#> 
#> [[150]]$external_name
#> [1] "SNORD117"
#> 
#> [[150]]$description
#> [1] "small nucleolar RNA, C/D box 117 [Source:HGNC Symbol;Acc:32742]"
#> 
#> [[150]]$feature_type
#> [1] "gene"
#> 
#> [[150]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[150]]$biotype
#> [1] "snoRNA"
#> 
#> [[150]]$seq_region_name
#> [1] "6"
#> 
#> [[150]]$source
#> [1] "ensembl"
#> 
#> 
#> [[151]]
#> [[151]]$seq_region_name
#> [1] "6"
#> 
#> [[151]]$source
#> [1] "ensembl"
#> 
#> [[151]]$description
#> [1] "small nucleolar RNA, C/D box 84 [Source:HGNC Symbol;Acc:32743]"
#> 
#> [[151]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[151]]$feature_type
#> [1] "gene"
#> 
#> [[151]]$biotype
#> [1] "snoRNA"
#> 
#> [[151]]$start
#> [1] 31508878
#> 
#> [[151]]$id
#> [1] "ENSG00000265236"
#> 
#> [[151]]$end
#> [1] 31508955
#> 
#> [[151]]$external_name
#> [1] "SNORD84"
#> 
#> [[151]]$canonical_transcript
#> [1] "ENST00000584275.1"
#> 
#> [[151]]$gene_id
#> [1] "ENSG00000265236"
#> 
#> [[151]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[151]]$strand
#> [1] -1
#> 
#> [[151]]$version
#> [1] 1
#> 
#> 
#> [[152]]
#> [[152]]$start
#> [1] 31510081
#> 
#> [[152]]$id
#> [1] "ENSG00000234006"
#> 
#> [[152]]$end
#> [1] 31510915
#> 
#> [[152]]$canonical_transcript
#> [1] "ENST00000420520.1"
#> 
#> [[152]]$external_name
#> [1] "DDX39B-AS1"
#> 
#> [[152]]$gene_id
#> [1] "ENSG00000234006"
#> 
#> [[152]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[152]]$strand
#> [1] 1
#> 
#> [[152]]$version
#> [1] 1
#> 
#> [[152]]$seq_region_name
#> [1] "6"
#> 
#> [[152]]$source
#> [1] "havana"
#> 
#> [[152]]$description
#> [1] "DDX39B antisense RNA 1 [Source:HGNC Symbol;Acc:39771]"
#> 
#> [[152]]$feature_type
#> [1] "gene"
#> 
#> [[152]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[152]]$biotype
#> [1] "antisense"
#> 
#> 
#> [[153]]
#> [[153]]$seq_region_name
#> [1] "6"
#> 
#> [[153]]$source
#> [1] "havana"
#> 
#> [[153]]$feature_type
#> [1] "gene"
#> 
#> [[153]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[153]]$biotype
#> [1] "pseudogene"
#> 
#> [[153]]$description
#> [1] "ubiquinol-cytochrome c reductase hinge protein pseudogene 1 [Source:HGNC Symbol;Acc:31341]"
#> 
#> [[153]]$end
#> [1] 31579133
#> 
#> [[153]]$id
#> [1] "ENSG00000230622"
#> 
#> [[153]]$canonical_transcript
#> [1] "ENST00000424609.1"
#> 
#> [[153]]$external_name
#> [1] "UQCRHP1"
#> 
#> [[153]]$start
#> [1] 31578860
#> 
#> [[153]]$gene_id
#> [1] "ENSG00000230622"
#> 
#> [[153]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[153]]$version
#> [1] 1
#> 
#> [[153]]$strand
#> [1] -1
#> 
#> 
#> [[154]]
#> [[154]]$gene_id
#> [1] "ENSG00000200816"
#> 
#> [[154]]$version
#> [1] 1
#> 
#> [[154]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[154]]$strand
#> [1] 1
#> 
#> [[154]]$id
#> [1] "ENSG00000200816"
#> 
#> [[154]]$end
#> [1] 31590987
#> 
#> [[154]]$canonical_transcript
#> [1] "ENST00000363946.1"
#> 
#> [[154]]$external_name
#> [1] "SNORA38"
#> 
#> [[154]]$start
#> [1] 31590856
#> 
#> [[154]]$feature_type
#> [1] "gene"
#> 
#> [[154]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[154]]$biotype
#> [1] "snoRNA"
#> 
#> [[154]]$description
#> [1] "small nucleolar RNA, H/ACA box 38 [Source:HGNC Symbol;Acc:32631]"
#> 
#> [[154]]$seq_region_name
#> [1] "6"
#> 
#> [[154]]$source
#> [1] "ensembl"
#> 
#> 
#> [[155]]
#> [[155]]$end
#> [1] 31628498
#> 
#> [[155]]$id
#> [1] "ENSG00000227198"
#> 
#> [[155]]$canonical_transcript
#> [1] "ENST00000422049.1"
#> 
#> [[155]]$external_name
#> [1] "C6orf47-AS1"
#> 
#> [[155]]$start
#> [1] 31626106
#> 
#> [[155]]$gene_id
#> [1] "ENSG00000227198"
#> 
#> [[155]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[155]]$version
#> [1] 1
#> 
#> [[155]]$strand
#> [1] 1
#> 
#> [[155]]$seq_region_name
#> [1] "6"
#> 
#> [[155]]$source
#> [1] "havana"
#> 
#> [[155]]$feature_type
#> [1] "gene"
#> 
#> [[155]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[155]]$biotype
#> [1] "antisense"
#> 
#> [[155]]$description
#> [1] "C6orf47 antisense RNA 1 [Source:HGNC Symbol;Acc:39767]"
#> 
#> 
#> [[156]]
#> [[156]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[156]]$biotype
#> [1] "misc_RNA"
#> 
#> [[156]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[156]]$feature_type
#> [1] "gene"
#> 
#> [[156]]$source
#> [1] "ensembl"
#> 
#> [[156]]$seq_region_name
#> [1] "6"
#> 
#> [[156]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[156]]$strand
#> [1] 1
#> 
#> [[156]]$version
#> [1] 1
#> 
#> [[156]]$gene_id
#> [1] "ENSG00000201207"
#> 
#> [[156]]$start
#> [1] 31631065
#> 
#> [[156]]$external_name
#> [1] "Y_RNA"
#> 
#> [[156]]$canonical_transcript
#> [1] "ENST00000364337.1"
#> 
#> [[156]]$end
#> [1] 31631178
#> 
#> [[156]]$id
#> [1] "ENSG00000201207"
#> 
#> 
#> [[157]]
#> [[157]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[157]]$version
#> [1] 1
#> 
#> [[157]]$strand
#> [1] -1
#> 
#> [[157]]$gene_id
#> [1] "ENSG00000266776"
#> 
#> [[157]]$canonical_transcript
#> [1] "ENST00000580775.1"
#> 
#> [[157]]$external_name
#> [1] "MIR4646"
#> 
#> [[157]]$id
#> [1] "ENSG00000266776"
#> 
#> [[157]]$end
#> [1] 31668868
#> 
#> [[157]]$start
#> [1] 31668806
#> 
#> [[157]]$biotype
#> [1] "miRNA"
#> 
#> [[157]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[157]]$feature_type
#> [1] "gene"
#> 
#> [[157]]$description
#> [1] "microRNA 4646 [Source:HGNC Symbol;Acc:41543]"
#> 
#> [[157]]$source
#> [1] "ensembl"
#> 
#> [[157]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[158]]
#> [[158]]$seq_region_name
#> [1] "6"
#> 
#> [[158]]$source
#> [1] "ensembl"
#> 
#> [[158]]$feature_type
#> [1] "gene"
#> 
#> [[158]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[158]]$biotype
#> [1] "snRNA"
#> 
#> [[158]]$description
#> [1] "RNA, U6 small nuclear 850, pseudogene [Source:HGNC Symbol;Acc:47813]"
#> 
#> [[158]]$id
#> [1] "ENSG00000252743"
#> 
#> [[158]]$end
#> [1] 31724830
#> 
#> [[158]]$canonical_transcript
#> [1] "ENST00000516934.1"
#> 
#> [[158]]$external_name
#> [1] "RNU6-850P"
#> 
#> [[158]]$start
#> [1] 31724728
#> 
#> [[158]]$gene_id
#> [1] "ENSG00000252743"
#> 
#> [[158]]$strand
#> [1] -1
#> 
#> [[158]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[158]]$version
#> [1] 1
#> 
#> 
#> [[159]]
#> [[159]]$seq_region_name
#> [1] "6"
#> 
#> [[159]]$source
#> [1] "havana"
#> 
#> [[159]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[159]]$feature_type
#> [1] "gene"
#> 
#> [[159]]$biotype
#> [1] "antisense"
#> 
#> [[159]]$description
#> [1] "SAPCD1 antisense RNA 1 [Source:HGNC Symbol;Acc:39824]"
#> 
#> [[159]]$id
#> [1] "ENSG00000235663"
#> 
#> [[159]]$end
#> [1] 31733365
#> 
#> [[159]]$canonical_transcript
#> [1] "ENST00000419679.1"
#> 
#> [[159]]$external_name
#> [1] "SAPCD1-AS1"
#> 
#> [[159]]$start
#> [1] 31732087
#> 
#> [[159]]$gene_id
#> [1] "ENSG00000235663"
#> 
#> [[159]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[159]]$strand
#> [1] -1
#> 
#> [[159]]$version
#> [1] 1
#> 
#> 
#> [[160]]
#> [[160]]$biotype
#> [1] "misc_RNA"
#> 
#> [[160]]$feature_type
#> [1] "gene"
#> 
#> [[160]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[160]]$description
#> [1] "Y RNA [Source:RFAM;Acc:RF00019]"
#> 
#> [[160]]$source
#> [1] "ensembl"
#> 
#> [[160]]$seq_region_name
#> [1] "6"
#> 
#> [[160]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[160]]$strand
#> [1] -1
#> 
#> [[160]]$version
#> [1] 1
#> 
#> [[160]]$gene_id
#> [1] "ENSG00000201555"
#> 
#> [[160]]$canonical_transcript
#> [1] "ENST00000364685.1"
#> 
#> [[160]]$external_name
#> [1] "Y_RNA"
#> 
#> [[160]]$end
#> [1] 31746682
#> 
#> [[160]]$id
#> [1] "ENSG00000201555"
#> 
#> [[160]]$start
#> [1] 31746594
#> 
#> 
#> [[161]]
#> [[161]]$feature_type
#> [1] "gene"
#> 
#> [[161]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[161]]$biotype
#> [1] "snoRNA"
#> 
#> [[161]]$description
#> [1] "small nucleolar RNA, C/D box 48 [Source:HGNC Symbol;Acc:10188]"
#> 
#> [[161]]$seq_region_name
#> [1] "6"
#> 
#> [[161]]$source
#> [1] "ensembl"
#> 
#> [[161]]$gene_id
#> [1] "ENSG00000201823"
#> 
#> [[161]]$version
#> [1] 1
#> 
#> [[161]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[161]]$strand
#> [1] 1
#> 
#> [[161]]$end
#> [1] 31803103
#> 
#> [[161]]$id
#> [1] "ENSG00000201823"
#> 
#> [[161]]$external_name
#> [1] "SNORD48"
#> 
#> [[161]]$canonical_transcript
#> [1] "ENST00000364953.1"
#> 
#> [[161]]$start
#> [1] 31803040
#> 
#> 
#> [[162]]
#> [[162]]$canonical_transcript
#> [1] "ENST00000364884.1"
#> 
#> [[162]]$external_name
#> [1] "SNORD52"
#> 
#> [[162]]$end
#> [1] 31804919
#> 
#> [[162]]$id
#> [1] "ENSG00000201754"
#> 
#> [[162]]$start
#> [1] 31804853
#> 
#> [[162]]$strand
#> [1] 1
#> 
#> [[162]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[162]]$version
#> [1] 1
#> 
#> [[162]]$gene_id
#> [1] "ENSG00000201754"
#> 
#> [[162]]$source
#> [1] "ensembl"
#> 
#> [[162]]$seq_region_name
#> [1] "6"
#> 
#> [[162]]$biotype
#> [1] "snoRNA"
#> 
#> [[162]]$feature_type
#> [1] "gene"
#> 
#> [[162]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[162]]$description
#> [1] "small nucleolar RNA, C/D box 52 [Source:HGNC Symbol;Acc:10202]"
#> 
#> 
#> [[163]]
#> [[163]]$source
#> [1] "havana"
#> 
#> [[163]]$seq_region_name
#> [1] "6"
#> 
#> [[163]]$description
#> [1] "EHMT2 antisense RNA 1 [Source:HGNC Symbol;Acc:39751]"
#> 
#> [[163]]$biotype
#> [1] "antisense"
#> 
#> [[163]]$feature_type
#> [1] "gene"
#> 
#> [[163]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[163]]$start
#> [1] 31851538
#> 
#> [[163]]$canonical_transcript
#> [1] "ENST00000434689.1"
#> 
#> [[163]]$external_name
#> [1] "EHMT2-AS1"
#> 
#> [[163]]$end
#> [1] 31851831
#> 
#> [[163]]$id
#> [1] "ENSG00000237080"
#> 
#> [[163]]$strand
#> [1] 1
#> 
#> [[163]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[163]]$version
#> [1] 1
#> 
#> [[163]]$gene_id
#> [1] "ENSG00000237080"
#> 
#> 
#> [[164]]
#> [[164]]$source
#> [1] "ensembl"
#> 
#> [[164]]$seq_region_name
#> [1] "6"
#> 
#> [[164]]$biotype
#> [1] "miRNA"
#> 
#> [[164]]$feature_type
#> [1] "gene"
#> 
#> [[164]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[164]]$description
#> [1] "microRNA 1236 [Source:HGNC Symbol;Acc:33925]"
#> 
#> [[164]]$external_name
#> [1] "MIR1236"
#> 
#> [[164]]$canonical_transcript
#> [1] "ENST00000408340.1"
#> 
#> [[164]]$end
#> [1] 31924717
#> 
#> [[164]]$id
#> [1] "ENSG00000221267"
#> 
#> [[164]]$start
#> [1] 31924616
#> 
#> [[164]]$strand
#> [1] -1
#> 
#> [[164]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[164]]$version
#> [1] 1
#> 
#> [[164]]$gene_id
#> [1] "ENSG00000221267"
#> 
#> 
#> [[165]]
#> [[165]]$seq_region_name
#> [1] "6"
#> 
#> [[165]]$source
#> [1] "havana"
#> 
#> [[165]]$description
#> [1] "C4A antisense RNA 1 [Source:HGNC Symbol;Acc:39753]"
#> 
#> [[165]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[165]]$feature_type
#> [1] "gene"
#> 
#> [[165]]$biotype
#> [1] "antisense"
#> 
#> [[165]]$start
#> [1] 31967753
#> 
#> [[165]]$id
#> [1] "ENSG00000233627"
#> 
#> [[165]]$end
#> [1] 31971298
#> 
#> [[165]]$canonical_transcript
#> [1] "ENST00000458633.1"
#> 
#> [[165]]$external_name
#> [1] "C4A-AS1"
#> 
#> [[165]]$gene_id
#> [1] "ENSG00000233627"
#> 
#> [[165]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[165]]$strand
#> [1] -1
#> 
#> [[165]]$version
#> [1] 2
#> 
#> 
#> [[166]]
#> [[166]]$feature_type
#> [1] "gene"
#> 
#> [[166]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[166]]$biotype
#> [1] "pseudogene"
#> 
#> [[166]]$description
#> [1] "cytochrome P450, family 21, subfamily A, polypeptide 1 pseudogene [Source:HGNC Symbol;Acc:2599]"
#> 
#> [[166]]$seq_region_name
#> [1] "6"
#> 
#> [[166]]$source
#> [1] "havana"
#> 
#> [[166]]$gene_id
#> [1] "ENSG00000204338"
#> 
#> [[166]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[166]]$strand
#> [1] 1
#> 
#> [[166]]$version
#> [1] 4
#> 
#> [[166]]$id
#> [1] "ENSG00000204338"
#> 
#> [[166]]$end
#> [1] 31976228
#> 
#> [[166]]$external_name
#> [1] "CYP21A1P"
#> 
#> [[166]]$canonical_transcript
#> [1] "ENST00000342991.6"
#> 
#> [[166]]$start
#> [1] 31973413
#> 
#> 
#> [[167]]
#> [[167]]$gene_id
#> [1] "ENSG00000248290"
#> 
#> [[167]]$version
#> [1] 1
#> 
#> [[167]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[167]]$strand
#> [1] -1
#> 
#> [[167]]$start
#> [1] 31976391
#> 
#> [[167]]$id
#> [1] "ENSG00000248290"
#> 
#> [[167]]$end
#> [1] 31980249
#> 
#> [[167]]$canonical_transcript
#> [1] "ENST00000507684.1"
#> 
#> [[167]]$external_name
#> [1] "TNXA"
#> 
#> [[167]]$description
#> [1] "tenascin XA (pseudogene) [Source:HGNC Symbol;Acc:11975]"
#> 
#> [[167]]$feature_type
#> [1] "gene"
#> 
#> [[167]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[167]]$biotype
#> [1] "pseudogene"
#> 
#> [[167]]$seq_region_name
#> [1] "6"
#> 
#> [[167]]$source
#> [1] "havana"
#> 
#> 
#> [[168]]
#> [[168]]$description
#> [1] "serine/threonine kinase 19 pseudogene [Source:HGNC Symbol;Acc:21668]"
#> 
#> [[168]]$feature_type
#> [1] "gene"
#> 
#> [[168]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[168]]$biotype
#> [1] "pseudogene"
#> 
#> [[168]]$seq_region_name
#> [1] "6"
#> 
#> [[168]]$source
#> [1] "havana"
#> 
#> [[168]]$gene_id
#> [1] "ENSG00000250535"
#> 
#> [[168]]$strand
#> [1] 1
#> 
#> [[168]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[168]]$version
#> [1] 1
#> 
#> [[168]]$start
#> [1] 31981047
#> 
#> [[168]]$id
#> [1] "ENSG00000250535"
#> 
#> [[168]]$end
#> [1] 31981564
#> 
#> [[168]]$canonical_transcript
#> [1] "ENST00000512515.1"
#> 
#> [[168]]$external_name
#> [1] "STK19P"
#> 
#> 
#> [[169]]
#> [[169]]$description
#> [1] "C4B antisense RNA 1 [Source:HGNC Symbol;Acc:39752]"
#> 
#> [[169]]$biotype
#> [1] "antisense"
#> 
#> [[169]]$feature_type
#> [1] "gene"
#> 
#> [[169]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[169]]$source
#> [1] "havana"
#> 
#> [[169]]$seq_region_name
#> [1] "6"
#> 
#> [[169]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[169]]$version
#> [1] 1
#> 
#> [[169]]$strand
#> [1] -1
#> 
#> [[169]]$gene_id
#> [1] "ENSG00000229776"
#> 
#> [[169]]$start
#> [1] 32000490
#> 
#> [[169]]$external_name
#> [1] "C4B-AS1"
#> 
#> [[169]]$canonical_transcript
#> [1] "ENST00000415626.1"
#> 
#> [[169]]$end
#> [1] 32004035
#> 
#> [[169]]$id
#> [1] "ENSG00000229776"
#> 
#> 
#> [[170]]
#> [[170]]$seq_region_name
#> [1] "6"
#> 
#> [[170]]$source
#> [1] "ensembl"
#> 
#> [[170]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[170]]$feature_type
#> [1] "gene"
#> 
#> [[170]]$biotype
#> [1] "rRNA"
#> 
#> [[170]]$description
#> [1] "RNA, 5S ribosomal pseudogene 206 [Source:HGNC Symbol;Acc:43106]"
#> 
#> [[170]]$end
#> [1] 32046405
#> 
#> [[170]]$id
#> [1] "ENSG00000252512"
#> 
#> [[170]]$external_name
#> [1] "RNA5SP206"
#> 
#> [[170]]$canonical_transcript
#> [1] "ENST00000516703.1"
#> 
#> [[170]]$start
#> [1] 32046285
#> 
#> [[170]]$gene_id
#> [1] "ENSG00000252512"
#> 
#> [[170]]$version
#> [1] 1
#> 
#> [[170]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[170]]$strand
#> [1] -1
#> 
#> 
#> [[171]]
#> [[171]]$description
#> NULL
#> 
#> [[171]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[171]]$feature_type
#> [1] "gene"
#> 
#> [[171]]$biotype
#> [1] "3prime_overlapping_ncRNA"
#> 
#> [[171]]$seq_region_name
#> [1] "6"
#> 
#> [[171]]$source
#> [1] "havana"
#> 
#> [[171]]$gene_id
#> [1] "ENSG00000273333"
#> 
#> [[171]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[171]]$strand
#> [1] -1
#> 
#> [[171]]$version
#> [1] 1
#> 
#> [[171]]$start
#> [1] 32152516
#> 
#> [[171]]$id
#> [1] "ENSG00000273333"
#> 
#> [[171]]$end
#> [1] 32153659
#> 
#> [[171]]$external_name
#> [1] "XXbac-BPG300A18.13"
#> 
#> [[171]]$canonical_transcript
#> [1] "ENST00000559458.1"
#> 
#> 
#> [[172]]
#> [[172]]$external_name
#> [1] "XXbac-BPG154L12.4"
#> 
#> [[172]]$canonical_transcript
#> [1] "ENST00000425033.1"
#> 
#> [[172]]$id
#> [1] "ENSG00000225914"
#> 
#> [[172]]$end
#> [1] 32233615
#> 
#> [[172]]$start
#> [1] 32223488
#> 
#> [[172]]$version
#> [1] 1
#> 
#> [[172]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[172]]$strand
#> [1] 1
#> 
#> [[172]]$gene_id
#> [1] "ENSG00000225914"
#> 
#> [[172]]$source
#> [1] "havana"
#> 
#> [[172]]$seq_region_name
#> [1] "6"
#> 
#> [[172]]$biotype
#> [1] "antisense"
#> 
#> [[172]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[172]]$feature_type
#> [1] "gene"
#> 
#> [[172]]$description
#> NULL
#> 
#> 
#> [[173]]
#> [[173]]$version
#> [1] 1
#> 
#> [[173]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[173]]$strand
#> [1] 1
#> 
#> [[173]]$gene_id
#> [1] "ENSG00000237285"
#> 
#> [[173]]$start
#> [1] 32292996
#> 
#> [[173]]$canonical_transcript
#> [1] "ENST00000432375.1"
#> 
#> [[173]]$external_name
#> [1] "HNRNPA1P2"
#> 
#> [[173]]$id
#> [1] "ENSG00000237285"
#> 
#> [[173]]$end
#> [1] 32293955
#> 
#> [[173]]$description
#> [1] "heterogeneous nuclear ribonucleoprotein A1 pseudogene 2 [Source:HGNC Symbol;Acc:13958]"
#> 
#> [[173]]$biotype
#> [1] "pseudogene"
#> 
#> [[173]]$feature_type
#> [1] "gene"
#> 
#> [[173]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[173]]$source
#> [1] "havana"
#> 
#> [[173]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[174]]
#> [[174]]$start
#> [1] 32320654
#> 
#> [[174]]$canonical_transcript
#> [1] "ENST00000411403.1"
#> 
#> [[174]]$external_name
#> [1] "RNU6-603P"
#> 
#> [[174]]$id
#> [1] "ENSG00000223335"
#> 
#> [[174]]$end
#> [1] 32320760
#> 
#> [[174]]$strand
#> [1] 1
#> 
#> [[174]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[174]]$version
#> [1] 1
#> 
#> [[174]]$gene_id
#> [1] "ENSG00000223335"
#> 
#> [[174]]$source
#> [1] "ensembl"
#> 
#> [[174]]$seq_region_name
#> [1] "6"
#> 
#> [[174]]$description
#> [1] "RNA, U6 small nuclear 603, pseudogene [Source:HGNC Symbol;Acc:47566]"
#> 
#> [[174]]$biotype
#> [1] "snRNA"
#> 
#> [[174]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[174]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[175]]
#> [[175]]$source
#> [1] "havana"
#> 
#> [[175]]$seq_region_name
#> [1] "6"
#> 
#> [[175]]$description
#> [1] "HLA complex group 23 (non-protein coding) [Source:HGNC Symbol;Acc:19713]"
#> 
#> [[175]]$biotype
#> [1] "antisense"
#> 
#> [[175]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[175]]$feature_type
#> [1] "gene"
#> 
#> [[175]]$start
#> [1] 32358287
#> 
#> [[175]]$external_name
#> [1] "HCG23"
#> 
#> [[175]]$canonical_transcript
#> [1] "ENST00000426643.1"
#> 
#> [[175]]$id
#> [1] "ENSG00000228962"
#> 
#> [[175]]$end
#> [1] 32361463
#> 
#> [[175]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[175]]$version
#> [1] 1
#> 
#> [[175]]$strand
#> [1] 1
#> 
#> [[175]]$gene_id
#> [1] "ENSG00000228962"
#> 
#> 
#> [[176]]
#> [[176]]$source
#> [1] "havana"
#> 
#> [[176]]$seq_region_name
#> [1] "6"
#> 
#> [[176]]$biotype
#> [1] "pseudogene"
#> 
#> [[176]]$feature_type
#> [1] "gene"
#> 
#> [[176]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[176]]$description
#> [1] "major histocompatibility complex, class II, DR beta 9 (pseudogene) [Source:HGNC Symbol;Acc:4957]"
#> 
#> [[176]]$canonical_transcript
#> [1] "ENST00000449413.1"
#> 
#> [[176]]$external_name
#> [1] "HLA-DRB9"
#> 
#> [[176]]$id
#> [1] "ENSG00000196301"
#> 
#> [[176]]$end
#> [1] 32441277
#> 
#> [[176]]$start
#> [1] 32427598
#> 
#> [[176]]$version
#> [1] 3
#> 
#> [[176]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[176]]$strand
#> [1] -1
#> 
#> [[176]]$gene_id
#> [1] "ENSG00000196301"
#> 
#> 
#> [[177]]
#> [[177]]$source
#> [1] "ensembl"
#> 
#> [[177]]$seq_region_name
#> [1] "6"
#> 
#> [[177]]$biotype
#> [1] "snRNA"
#> 
#> [[177]]$feature_type
#> [1] "gene"
#> 
#> [[177]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[177]]$description
#> [1] "RNA, U1 small nuclear 61, pseudogene [Source:HGNC Symbol;Acc:48403]"
#> 
#> [[177]]$external_name
#> [1] "RNU1-61P"
#> 
#> [[177]]$canonical_transcript
#> [1] "ENST00000516107.1"
#> 
#> [[177]]$id
#> [1] "ENSG00000251916"
#> 
#> [[177]]$end
#> [1] 32517867
#> 
#> [[177]]$start
#> [1] 32517717
#> 
#> [[177]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[177]]$strand
#> [1] -1
#> 
#> [[177]]$version
#> [1] 1
#> 
#> [[177]]$gene_id
#> [1] "ENSG00000251916"
#> 
#> 
#> [[178]]
#> [[178]]$id
#> [1] "ENSG00000229391"
#> 
#> [[178]]$end
#> [1] 32527799
#> 
#> [[178]]$canonical_transcript
#> [1] "ENST00000411500.1"
#> 
#> [[178]]$external_name
#> [1] "HLA-DRB6"
#> 
#> [[178]]$start
#> [1] 32520490
#> 
#> [[178]]$gene_id
#> [1] "ENSG00000229391"
#> 
#> [[178]]$strand
#> [1] -1
#> 
#> [[178]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[178]]$version
#> [1] 3
#> 
#> [[178]]$seq_region_name
#> [1] "6"
#> 
#> [[178]]$source
#> [1] "ensembl_havana"
#> 
#> [[178]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[178]]$feature_type
#> [1] "gene"
#> 
#> [[178]]$biotype
#> [1] "pseudogene"
#> 
#> [[178]]$description
#> [1] "major histocompatibility complex, class II, DR beta 6 (pseudogene) [Source:HGNC Symbol;Acc:4954]"
#> 
#> 
#> [[179]]
#> [[179]]$biotype
#> [1] "antisense"
#> 
#> [[179]]$feature_type
#> [1] "gene"
#> 
#> [[179]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[179]]$description
#> [1] "HLA-DQB1 antisense RNA 1 [Source:HGNC Symbol;Acc:39762]"
#> 
#> [[179]]$source
#> [1] "havana"
#> 
#> [[179]]$seq_region_name
#> [1] "6"
#> 
#> [[179]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[179]]$strand
#> [1] 1
#> 
#> [[179]]$version
#> [1] 1
#> 
#> [[179]]$gene_id
#> [1] "ENSG00000223534"
#> 
#> [[179]]$canonical_transcript
#> [1] "ENST00000419852.1"
#> 
#> [[179]]$external_name
#> [1] "HLA-DQB1-AS1"
#> 
#> [[179]]$id
#> [1] "ENSG00000223534"
#> 
#> [[179]]$end
#> [1] 32628506
#> 
#> [[179]]$start
#> [1] 32627657
#> 
#> 
#> [[180]]
#> [[180]]$gene_id
#> [1] "ENSG00000241287"
#> 
#> [[180]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[180]]$strand
#> [1] -1
#> 
#> [[180]]$version
#> [1] 1
#> 
#> [[180]]$start
#> [1] 32627663
#> 
#> [[180]]$end
#> [1] 32630227
#> 
#> [[180]]$id
#> [1] "ENSG00000241287"
#> 
#> [[180]]$external_name
#> [1] "XXbac-BPG254F23.6"
#> 
#> [[180]]$canonical_transcript
#> [1] "ENST00000443574.1"
#> 
#> [[180]]$description
#> NULL
#> 
#> [[180]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[180]]$feature_type
#> [1] "gene"
#> 
#> [[180]]$biotype
#> [1] "processed_transcript"
#> 
#> [[180]]$seq_region_name
#> [1] "6"
#> 
#> [[180]]$source
#> [1] "havana"
#> 
#> 
#> [[181]]
#> [[181]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[181]]$version
#> [1] 1
#> 
#> [[181]]$strand
#> [1] -1
#> 
#> [[181]]$gene_id
#> [1] "ENSG00000235040"
#> 
#> [[181]]$start
#> [1] 32673901
#> 
#> [[181]]$canonical_transcript
#> [1] "ENST00000422088.1"
#> 
#> [[181]]$external_name
#> [1] "MTCO3P1"
#> 
#> [[181]]$id
#> [1] "ENSG00000235040"
#> 
#> [[181]]$end
#> [1] 32674732
#> 
#> [[181]]$description
#> [1] "MT-CO3 pseudogene 1 [Source:HGNC Symbol;Acc:31342]"
#> 
#> [[181]]$biotype
#> [1] "pseudogene"
#> 
#> [[181]]$feature_type
#> [1] "gene"
#> 
#> [[181]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[181]]$source
#> [1] "havana"
#> 
#> [[181]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[182]]
#> [[182]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[182]]$version
#> [1] 3
#> 
#> [[182]]$strand
#> [1] 1
#> 
#> [[182]]$gene_id
#> [1] "ENSG00000232080"
#> 
#> [[182]]$start
#> [1] 32685782
#> 
#> [[182]]$external_name
#> [1] "XXbac-BPG254F23.7"
#> 
#> [[182]]$canonical_transcript
#> [1] "ENST00000585372.1"
#> 
#> [[182]]$id
#> [1] "ENSG00000232080"
#> 
#> [[182]]$end
#> [1] 32686612
#> 
#> [[182]]$description
#> NULL
#> 
#> [[182]]$biotype
#> [1] "lincRNA"
#> 
#> [[182]]$feature_type
#> [1] "gene"
#> 
#> [[182]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[182]]$source
#> [1] "havana"
#> 
#> [[182]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[183]]
#> [[183]]$start
#> [1] 32698535
#> 
#> [[183]]$canonical_transcript
#> [1] "ENST00000425319.1"
#> 
#> [[183]]$external_name
#> [1] "HLA-DQB3"
#> 
#> [[183]]$id
#> [1] "ENSG00000226030"
#> 
#> [[183]]$end
#> [1] 32699472
#> 
#> [[183]]$strand
#> [1] -1
#> 
#> [[183]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[183]]$version
#> [1] 1
#> 
#> [[183]]$gene_id
#> [1] "ENSG00000226030"
#> 
#> [[183]]$source
#> [1] "havana"
#> 
#> [[183]]$seq_region_name
#> [1] "6"
#> 
#> [[183]]$description
#> [1] "major histocompatibility complex, class II, DQ beta 3 [Source:HGNC Symbol;Acc:4946]"
#> 
#> [[183]]$biotype
#> [1] "pseudogene"
#> 
#> [[183]]$feature_type
#> [1] "gene"
#> 
#> [[183]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> 
#> [[184]]
#> [[184]]$source
#> [1] "ensembl"
#> 
#> [[184]]$seq_region_name
#> [1] "6"
#> 
#> [[184]]$biotype
#> [1] "miRNA"
#> 
#> [[184]]$feature_type
#> [1] "gene"
#> 
#> [[184]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[184]]$description
#> [1] "microRNA 3135b [Source:HGNC Symbol;Acc:41783]"
#> 
#> [[184]]$canonical_transcript
#> [1] "ENST00000581098.1"
#> 
#> [[184]]$external_name
#> [1] "MIR3135B"
#> 
#> [[184]]$end
#> [1] 32717756
#> 
#> [[184]]$id
#> [1] "ENSG00000263649"
#> 
#> [[184]]$start
#> [1] 32717689
#> 
#> [[184]]$version
#> [1] 1
#> 
#> [[184]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[184]]$strand
#> [1] -1
#> 
#> [[184]]$gene_id
#> [1] "ENSG00000263649"
#> 
#> 
#> [[185]]
#> [[185]]$gene_id
#> [1] "ENSG00000204261"
#> 
#> [[185]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[185]]$version
#> [1] 4
#> 
#> [[185]]$strand
#> [1] 1
#> 
#> [[185]]$id
#> [1] "ENSG00000204261"
#> 
#> [[185]]$end
#> [1] 32814272
#> 
#> [[185]]$canonical_transcript
#> [1] "ENST00000453426.1"
#> 
#> [[185]]$external_name
#> [1] "TAPSAR1"
#> 
#> [[185]]$start
#> [1] 32811863
#> 
#> [[185]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[185]]$feature_type
#> [1] "gene"
#> 
#> [[185]]$biotype
#> [1] "lincRNA"
#> 
#> [[185]]$description
#> [1] "TAP1 and PSMB8 antisense RNA 1 [Source:HGNC Symbol;Acc:39758]"
#> 
#> [[185]]$seq_region_name
#> [1] "6"
#> 
#> [[185]]$source
#> [1] "havana"
#> 
#> 
#> [[186]]
#> [[186]]$gene_id
#> [1] "ENSG00000234515"
#> 
#> [[186]]$version
#> [1] 1
#> 
#> [[186]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[186]]$strand
#> [1] -1
#> 
#> [[186]]$end
#> [1] 32847625
#> 
#> [[186]]$id
#> [1] "ENSG00000234515"
#> 
#> [[186]]$external_name
#> [1] "PPP1R2P1"
#> 
#> [[186]]$canonical_transcript
#> [1] "ENST00000429032.1"
#> 
#> [[186]]$start
#> [1] 32846948
#> 
#> [[186]]$feature_type
#> [1] "gene"
#> 
#> [[186]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[186]]$biotype
#> [1] "pseudogene"
#> 
#> [[186]]$description
#> [1] "protein phosphatase 1, regulatory (inhibitor) subunit 2 pseudogene 1 [Source:HGNC Symbol;Acc:9289]"
#> 
#> [[186]]$seq_region_name
#> [1] "6"
#> 
#> [[186]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[187]]
#> [[187]]$feature_type
#> [1] "gene"
#> 
#> [[187]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[187]]$biotype
#> [1] "pseudogene"
#> 
#> [[187]]$description
#> [1] "major histocompatibility complex, class I, Z (pseudogene) [Source:HGNC Symbol;Acc:19394]"
#> 
#> [[187]]$seq_region_name
#> [1] "6"
#> 
#> [[187]]$source
#> [1] "havana"
#> 
#> [[187]]$gene_id
#> [1] "ENSG00000235301"
#> 
#> [[187]]$version
#> [1] 1
#> 
#> [[187]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[187]]$strand
#> [1] 1
#> 
#> [[187]]$end
#> [1] 32864267
#> 
#> [[187]]$id
#> [1] "ENSG00000235301"
#> 
#> [[187]]$canonical_transcript
#> [1] "ENST00000434920.1"
#> 
#> [[187]]$external_name
#> [1] "HLA-Z"
#> 
#> [[187]]$start
#> [1] 32864193
#> 
#> 
#> [[188]]
#> [[188]]$gene_id
#> [1] "ENSG00000212066"
#> 
#> [[188]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[188]]$version
#> [1] 1
#> 
#> [[188]]$strand
#> [1] -1
#> 
#> [[188]]$start
#> [1] 32904693
#> 
#> [[188]]$id
#> [1] "ENSG00000212066"
#> 
#> [[188]]$end
#> [1] 32904787
#> 
#> [[188]]$external_name
#> [1] "AL645941.1"
#> 
#> [[188]]$canonical_transcript
#> [1] "ENST00000390777.1"
#> 
#> [[188]]$description
#> NULL
#> 
#> [[188]]$feature_type
#> [1] "gene"
#> 
#> [[188]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[188]]$biotype
#> [1] "miRNA"
#> 
#> [[188]]$seq_region_name
#> [1] "6"
#> 
#> [[188]]$source
#> [1] "ensembl"
#> 
#> 
#> [[189]]
#> [[189]]$start
#> [1] 32938009
#> 
#> [[189]]$id
#> [1] "ENSG00000223837"
#> 
#> [[189]]$end
#> [1] 32938663
#> 
#> [[189]]$external_name
#> [1] "BRD2-IT1"
#> 
#> [[189]]$canonical_transcript
#> [1] "ENST00000415875.2"
#> 
#> [[189]]$gene_id
#> [1] "ENSG00000223837"
#> 
#> [[189]]$version
#> [1] 2
#> 
#> [[189]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[189]]$strand
#> [1] 1
#> 
#> [[189]]$seq_region_name
#> [1] "6"
#> 
#> [[189]]$source
#> [1] "havana"
#> 
#> [[189]]$description
#> [1] "BRD2 intronic transcript 1 (non-protein coding) [Source:HGNC Symbol;Acc:41311]"
#> 
#> [[189]]$feature_type
#> [1] "gene"
#> 
#> [[189]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[189]]$biotype
#> [1] "sense_intronic"
#> 
#> 
#> [[190]]
#> [[190]]$source
#> [1] "havana"
#> 
#> [[190]]$seq_region_name
#> [1] "6"
#> 
#> [[190]]$description
#> NULL
#> 
#> [[190]]$biotype
#> [1] "antisense"
#> 
#> [[190]]$feature_type
#> [1] "gene"
#> 
#> [[190]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[190]]$start
#> [1] 32939842
#> 
#> [[190]]$canonical_transcript
#> [1] "ENST00000580587.1"
#> 
#> [[190]]$external_name
#> [1] "XXbac-BPG181M17.6"
#> 
#> [[190]]$id
#> [1] "ENSG00000263756"
#> 
#> [[190]]$end
#> [1] 32940630
#> 
#> [[190]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[190]]$version
#> [1] 1
#> 
#> [[190]]$strand
#> [1] -1
#> 
#> [[190]]$gene_id
#> [1] "ENSG00000263756"
#> 
#> 
#> [[191]]
#> [[191]]$start
#> [1] 33047228
#> 
#> [[191]]$id
#> [1] "ENSG00000224796"
#> 
#> [[191]]$end
#> [1] 33047637
#> 
#> [[191]]$canonical_transcript
#> [1] "ENST00000439737.1"
#> 
#> [[191]]$external_name
#> [1] "RPL32P1"
#> 
#> [[191]]$gene_id
#> [1] "ENSG00000224796"
#> 
#> [[191]]$strand
#> [1] 1
#> 
#> [[191]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[191]]$version
#> [1] 1
#> 
#> [[191]]$seq_region_name
#> [1] "6"
#> 
#> [[191]]$source
#> [1] "havana"
#> 
#> [[191]]$description
#> [1] "ribosomal protein L32 pseudogene 1 [Source:HGNC Symbol;Acc:10339]"
#> 
#> [[191]]$feature_type
#> [1] "gene"
#> 
#> [[191]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[191]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[192]]
#> [[192]]$gene_id
#> [1] "ENSG00000231461"
#> 
#> [[192]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[192]]$strand
#> [1] -1
#> 
#> [[192]]$version
#> [1] 1
#> 
#> [[192]]$end
#> [1] 33065072
#> 
#> [[192]]$id
#> [1] "ENSG00000231461"
#> 
#> [[192]]$canonical_transcript
#> [1] "ENST00000433582.1"
#> 
#> [[192]]$external_name
#> [1] "HLA-DPA2"
#> 
#> [[192]]$start
#> [1] 33059530
#> 
#> [[192]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[192]]$feature_type
#> [1] "gene"
#> 
#> [[192]]$biotype
#> [1] "pseudogene"
#> 
#> [[192]]$description
#> [1] "major histocompatibility complex, class II, DP alpha 2 (pseudogene) [Source:HGNC Symbol;Acc:4939]"
#> 
#> [[192]]$seq_region_name
#> [1] "6"
#> 
#> [[192]]$source
#> [1] "havana"
#> 
#> 
#> [[193]]
#> [[193]]$seq_region_name
#> [1] "6"
#> 
#> [[193]]$source
#> [1] "havana"
#> 
#> [[193]]$feature_type
#> [1] "gene"
#> 
#> [[193]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[193]]$biotype
#> [1] "pseudogene"
#> 
#> [[193]]$description
#> [1] "collagen, type XI, alpha 2 pseudogene 1 [Source:HGNC Symbol;Acc:13947]"
#> 
#> [[193]]$end
#> [1] 33075107
#> 
#> [[193]]$id
#> [1] "ENSG00000228688"
#> 
#> [[193]]$canonical_transcript
#> [1] "ENST00000441798.1"
#> 
#> [[193]]$external_name
#> [1] "COL11A2P1"
#> 
#> [[193]]$start
#> [1] 33071571
#> 
#> [[193]]$gene_id
#> [1] "ENSG00000228688"
#> 
#> [[193]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[193]]$strand
#> [1] -1
#> 
#> [[193]]$version
#> [1] 1
#> 
#> 
#> [[194]]
#> [[194]]$gene_id
#> [1] "ENSG00000224557"
#> 
#> [[194]]$strand
#> [1] 1
#> 
#> [[194]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[194]]$version
#> [1] 3
#> 
#> [[194]]$end
#> [1] 33102442
#> 
#> [[194]]$id
#> [1] "ENSG00000224557"
#> 
#> [[194]]$external_name
#> [1] "HLA-DPB2"
#> 
#> [[194]]$canonical_transcript
#> [1] "ENST00000435499.2"
#> 
#> [[194]]$start
#> [1] 33080228
#> 
#> [[194]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[194]]$feature_type
#> [1] "gene"
#> 
#> [[194]]$biotype
#> [1] "pseudogene"
#> 
#> [[194]]$description
#> [1] "major histocompatibility complex, class II, DP beta 2 (pseudogene) [Source:HGNC Symbol;Acc:4941]"
#> 
#> [[194]]$seq_region_name
#> [1] "6"
#> 
#> [[194]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[195]]
#> [[195]]$start
#> [1] 33098993
#> 
#> [[195]]$id
#> [1] "ENSG00000237398"
#> 
#> [[195]]$end
#> [1] 33111102
#> 
#> [[195]]$canonical_transcript
#> [1] "ENST00000454398.1"
#> 
#> [[195]]$external_name
#> [1] "HLA-DPA3"
#> 
#> [[195]]$gene_id
#> [1] "ENSG00000237398"
#> 
#> [[195]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[195]]$version
#> [1] 1
#> 
#> [[195]]$strand
#> [1] -1
#> 
#> [[195]]$seq_region_name
#> [1] "6"
#> 
#> [[195]]$source
#> [1] "havana"
#> 
#> [[195]]$description
#> [1] "major histocompatibility complex, class II, DP alpha 3 (pseudogene) [Source:HGNC Symbol;Acc:19393]"
#> 
#> [[195]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[195]]$feature_type
#> [1] "gene"
#> 
#> [[195]]$biotype
#> [1] "pseudogene"
#> 
#> 
#> [[196]]
#> [[196]]$start
#> [1] 33112560
#> 
#> [[196]]$canonical_transcript
#> [1] "ENST00000414398.1"
#> 
#> [[196]]$external_name
#> [1] "HCG24"
#> 
#> [[196]]$end
#> [1] 33115544
#> 
#> [[196]]$id
#> [1] "ENSG00000230313"
#> 
#> [[196]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[196]]$strand
#> [1] 1
#> 
#> [[196]]$version
#> [1] 1
#> 
#> [[196]]$gene_id
#> [1] "ENSG00000230313"
#> 
#> [[196]]$source
#> [1] "havana"
#> 
#> [[196]]$seq_region_name
#> [1] "6"
#> 
#> [[196]]$description
#> [1] "HLA complex group 24 (non-protein coding) [Source:HGNC Symbol;Acc:23500]"
#> 
#> [[196]]$biotype
#> [1] "antisense"
#> 
#> [[196]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[196]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[197]]
#> [[197]]$seq_region_name
#> [1] "6"
#> 
#> [[197]]$source
#> [1] "ensembl"
#> 
#> [[197]]$description
#> [1] "RNA, Ro-associated Y4 pseudogene 10 [Source:HGNC Symbol;Acc:34060]"
#> 
#> [[197]]$feature_type
#> [1] "gene"
#> 
#> [[197]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[197]]$biotype
#> [1] "misc_RNA"
#> 
#> [[197]]$start
#> [1] 33167378
#> 
#> [[197]]$id
#> [1] "ENSG00000202441"
#> 
#> [[197]]$end
#> [1] 33167472
#> 
#> [[197]]$external_name
#> [1] "RNY4P10"
#> 
#> [[197]]$canonical_transcript
#> [1] "ENST00000365571.1"
#> 
#> [[197]]$gene_id
#> [1] "ENSG00000202441"
#> 
#> [[197]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[197]]$version
#> [1] 1
#> 
#> [[197]]$strand
#> [1] 1
#> 
#> 
#> [[198]]
#> [[198]]$description
#> [1] "microRNA 219-1 [Source:HGNC Symbol;Acc:31597]"
#> 
#> [[198]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[198]]$feature_type
#> [1] "gene"
#> 
#> [[198]]$biotype
#> [1] "miRNA"
#> 
#> [[198]]$seq_region_name
#> [1] "6"
#> 
#> [[198]]$source
#> [1] "ensembl"
#> 
#> [[198]]$gene_id
#> [1] "ENSG00000199036"
#> 
#> [[198]]$version
#> [1] 1
#> 
#> [[198]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[198]]$strand
#> [1] 1
#> 
#> [[198]]$start
#> [1] 33175612
#> 
#> [[198]]$id
#> [1] "ENSG00000199036"
#> 
#> [[198]]$end
#> [1] 33175721
#> 
#> [[198]]$canonical_transcript
#> [1] "ENST00000362166.1"
#> 
#> [[198]]$external_name
#> [1] "MIR219-1"
#> 
#> 
#> [[199]]
#> [[199]]$biotype
#> [1] "pseudogene"
#> 
#> [[199]]$feature_type
#> [1] "gene"
#> 
#> [[199]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[199]]$description
#> [1] "zinc finger protein 70 pseudogene 1 [Source:HGNC Symbol;Acc:13846]"
#> 
#> [[199]]$source
#> [1] "havana"
#> 
#> [[199]]$seq_region_name
#> [1] "6"
#> 
#> [[199]]$strand
#> [1] 1
#> 
#> [[199]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[199]]$version
#> [1] 1
#> 
#> [[199]]$gene_id
#> [1] "ENSG00000225463"
#> 
#> [[199]]$external_name
#> [1] "ZNF70P1"
#> 
#> [[199]]$canonical_transcript
#> [1] "ENST00000417480.1"
#> 
#> [[199]]$end
#> [1] 33184105
#> 
#> [[199]]$id
#> [1] "ENSG00000225463"
#> 
#> [[199]]$start
#> [1] 33183482
#> 
#> 
#> [[200]]
#> [[200]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[200]]$version
#> [1] 1
#> 
#> [[200]]$strand
#> [1] 1
#> 
#> [[200]]$gene_id
#> [1] "ENSG00000223457"
#> 
#> [[200]]$start
#> [1] 33205576
#> 
#> [[200]]$canonical_transcript
#> [1] "ENST00000432008.1"
#> 
#> [[200]]$external_name
#> [1] "HTATSF1P"
#> 
#> [[200]]$end
#> [1] 33207467
#> 
#> [[200]]$id
#> [1] "ENSG00000223457"
#> 
#> [[200]]$description
#> NULL
#> 
#> [[200]]$biotype
#> [1] "pseudogene"
#> 
#> [[200]]$feature_type
#> [1] "gene"
#> 
#> [[200]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[200]]$source
#> [1] "havana"
#> 
#> [[200]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[201]]
#> [[201]]$seq_region_name
#> [1] "6"
#> 
#> [[201]]$source
#> [1] "havana"
#> 
#> [[201]]$description
#> NULL
#> 
#> [[201]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[201]]$feature_type
#> [1] "gene"
#> 
#> [[201]]$biotype
#> [1] "processed_transcript"
#> 
#> [[201]]$start
#> [1] 33213852
#> 
#> [[201]]$end
#> [1] 33214633
#> 
#> [[201]]$id
#> [1] "ENSG00000272217"
#> 
#> [[201]]$canonical_transcript
#> [1] "ENST00000606432.1"
#> 
#> [[201]]$external_name
#> [1] "XXbac-BPG157A10.21"
#> 
#> [[201]]$gene_id
#> [1] "ENSG00000272217"
#> 
#> [[201]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[201]]$version
#> [1] 1
#> 
#> [[201]]$strand
#> [1] -1
#> 
#> 
#> [[202]]
#> [[202]]$source
#> [1] "havana"
#> 
#> [[202]]$seq_region_name
#> [1] "6"
#> 
#> [[202]]$biotype
#> [1] "antisense"
#> 
#> [[202]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[202]]$feature_type
#> [1] "gene"
#> 
#> [[202]]$description
#> [1] "HLA complex group 25 (non-protein coding) [Source:HGNC Symbol;Acc:20196]"
#> 
#> [[202]]$external_name
#> [1] "HCG25"
#> 
#> [[202]]$canonical_transcript
#> [1] "ENST00000450514.1"
#> 
#> [[202]]$end
#> [1] 33222766
#> 
#> [[202]]$id
#> [1] "ENSG00000232940"
#> 
#> [[202]]$start
#> [1] 33217311
#> 
#> [[202]]$strand
#> [1] 1
#> 
#> [[202]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[202]]$version
#> [1] 1
#> 
#> [[202]]$gene_id
#> [1] "ENSG00000232940"
#> 
#> 
#> [[203]]
#> [[203]]$feature_type
#> [1] "gene"
#> 
#> [[203]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[203]]$biotype
#> [1] "pseudogene"
#> 
#> [[203]]$description
#> [1] "myosin, light chain 8, pseudogene [Source:HGNC Symbol;Acc:7589]"
#> 
#> [[203]]$seq_region_name
#> [1] "6"
#> 
#> [[203]]$source
#> [1] "ensembl_havana"
#> 
#> [[203]]$gene_id
#> [1] "ENSG00000229596"
#> 
#> [[203]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[203]]$strand
#> [1] 1
#> 
#> [[203]]$version
#> [1] 2
#> 
#> [[203]]$end
#> [1] 33307272
#> 
#> [[203]]$id
#> [1] "ENSG00000229596"
#> 
#> [[203]]$canonical_transcript
#> [1] "ENST00000431574.1"
#> 
#> [[203]]$external_name
#> [1] "MYL8P"
#> 
#> [[203]]$start
#> [1] 33306755
#> 
#> 
#> [[204]]
#> [[204]]$seq_region_name
#> [1] "6"
#> 
#> [[204]]$source
#> [1] "ensembl_havana"
#> 
#> [[204]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[204]]$feature_type
#> [1] "gene"
#> 
#> [[204]]$biotype
#> [1] "pseudogene"
#> 
#> [[204]]$description
#> [1] "lysophospholipase II pseudogene 1 [Source:HGNC Symbol;Acc:21069]"
#> 
#> [[204]]$end
#> [1] 33334020
#> 
#> [[204]]$id
#> [1] "ENSG00000228285"
#> 
#> [[204]]$canonical_transcript
#> [1] "ENST00000447689.1"
#> 
#> [[204]]$external_name
#> [1] "LYPLA2P1"
#> 
#> [[204]]$start
#> [1] 33333325
#> 
#> [[204]]$gene_id
#> [1] "ENSG00000228285"
#> 
#> [[204]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[204]]$version
#> [1] 2
#> 
#> [[204]]$strand
#> [1] -1
#> 
#> 
#> [[205]]
#> [[205]]$description
#> [1] "ribosomal protein L35a pseudogene 4 [Source:HGNC Symbol;Acc:33460]"
#> 
#> [[205]]$biotype
#> [1] "pseudogene"
#> 
#> [[205]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[205]]$feature_type
#> [1] "gene"
#> 
#> [[205]]$source
#> [1] "havana"
#> 
#> [[205]]$seq_region_name
#> [1] "6"
#> 
#> [[205]]$version
#> [1] 1
#> 
#> [[205]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[205]]$strand
#> [1] -1
#> 
#> [[205]]$gene_id
#> [1] "ENSG00000225644"
#> 
#> [[205]]$start
#> [1] 33357151
#> 
#> [[205]]$canonical_transcript
#> [1] "ENST00000412007.1"
#> 
#> [[205]]$external_name
#> [1] "RPL35AP4"
#> 
#> [[205]]$id
#> [1] "ENSG00000225644"
#> 
#> [[205]]$end
#> [1] 33357252
#> 
#> 
#> [[206]]
#> [[206]]$seq_region_name
#> [1] "6"
#> 
#> [[206]]$source
#> [1] "ensembl_havana"
#> 
#> [[206]]$description
#> [1] "ribosomal protein L12 pseudogene 1 [Source:HGNC Symbol;Acc:13976]"
#> 
#> [[206]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[206]]$feature_type
#> [1] "gene"
#> 
#> [[206]]$biotype
#> [1] "pseudogene"
#> 
#> [[206]]$start
#> [1] 33367836
#> 
#> [[206]]$id
#> [1] "ENSG00000204194"
#> 
#> [[206]]$end
#> [1] 33368333
#> 
#> [[206]]$canonical_transcript
#> [1] "ENST00000374520.1"
#> 
#> [[206]]$external_name
#> [1] "RPL12P1"
#> 
#> [[206]]$gene_id
#> [1] "ENSG00000204194"
#> 
#> [[206]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[206]]$strand
#> [1] -1
#> 
#> [[206]]$version
#> [1] 1
#> 
#> 
#> [[207]]
#> [[207]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[207]]$strand
#> [1] 1
#> 
#> [[207]]$version
#> [1] 1
#> 
#> [[207]]$gene_id
#> [1] "ENSG00000264085"
#> 
#> [[207]]$canonical_transcript
#> [1] "ENST00000579078.1"
#> 
#> [[207]]$external_name
#> [1] "MIR5004"
#> 
#> [[207]]$end
#> [1] 33406214
#> 
#> [[207]]$id
#> [1] "ENSG00000264085"
#> 
#> [[207]]$start
#> [1] 33406108
#> 
#> [[207]]$biotype
#> [1] "miRNA"
#> 
#> [[207]]$feature_type
#> [1] "gene"
#> 
#> [[207]]$logic_name
#> [1] "ncrna_homo_sapiens_37"
#> 
#> [[207]]$description
#> [1] "microRNA 5004 [Source:HGNC Symbol;Acc:43532]"
#> 
#> [[207]]$source
#> [1] "ensembl"
#> 
#> [[207]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[208]]
#> [[208]]$start
#> [1] 31654739
#> 
#> [[208]]$id
#> [1] "ENSG00000204422"
#> 
#> [[208]]$end
#> [1] 31681849
#> 
#> [[208]]$external_name
#> [1] "XXbac-BPG32J3.20"
#> 
#> [[208]]$canonical_transcript
#> [1] "ENST00000461287.1"
#> 
#> [[208]]$gene_id
#> [1] "ENSG00000204422"
#> 
#> [[208]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[208]]$strand
#> [1] -1
#> 
#> [[208]]$version
#> [1] 7
#> 
#> [[208]]$seq_region_name
#> [1] "6"
#> 
#> [[208]]$source
#> [1] "havana"
#> 
#> [[208]]$description
#> [1] "Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:U3KQA8]"
#> 
#> [[208]]$feature_type
#> [1] "gene"
#> 
#> [[208]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[208]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[209]]
#> [[209]]$source
#> [1] "ensembl_havana"
#> 
#> [[209]]$seq_region_name
#> [1] "6"
#> 
#> [[209]]$description
#> [1] "lymphocyte antigen 6 complex, locus G6F [Source:HGNC Symbol;Acc:13933]"
#> 
#> [[209]]$biotype
#> [1] "protein_coding"
#> 
#> [[209]]$feature_type
#> [1] "gene"
#> 
#> [[209]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[209]]$start
#> [1] 31674640
#> 
#> [[209]]$canonical_transcript
#> [1] "ENST00000375832.4"
#> 
#> [[209]]$external_name
#> [1] "LY6G6F"
#> 
#> [[209]]$end
#> [1] 31685581
#> 
#> [[209]]$id
#> [1] "ENSG00000204424"
#> 
#> [[209]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[209]]$version
#> [1] 8
#> 
#> [[209]]$strand
#> [1] 1
#> 
#> [[209]]$gene_id
#> [1] "ENSG00000204424"
#> 
#> 
#> [[210]]
#> [[210]]$start
#> [1] 33239787
#> 
#> [[210]]$id
#> [1] "ENSG00000231500"
#> 
#> [[210]]$end
#> [1] 33244287
#> 
#> [[210]]$external_name
#> [1] "RPS18"
#> 
#> [[210]]$canonical_transcript
#> [1] "ENST00000439602.2"
#> 
#> [[210]]$gene_id
#> [1] "ENSG00000231500"
#> 
#> [[210]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[210]]$strand
#> [1] 1
#> 
#> [[210]]$version
#> [1] 2
#> 
#> [[210]]$seq_region_name
#> [1] "6"
#> 
#> [[210]]$source
#> [1] "ensembl_havana"
#> 
#> [[210]]$description
#> [1] "ribosomal protein S18 [Source:HGNC Symbol;Acc:10401]"
#> 
#> [[210]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[210]]$feature_type
#> [1] "gene"
#> 
#> [[210]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[211]]
#> [[211]]$description
#> [1] "UDP-Gal:betaGlcNAc beta 1,3-galactosyltransferase, polypeptide 4 [Source:HGNC Symbol;Acc:919]"
#> 
#> [[211]]$biotype
#> [1] "protein_coding"
#> 
#> [[211]]$feature_type
#> [1] "gene"
#> 
#> [[211]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[211]]$source
#> [1] "ensembl_havana"
#> 
#> [[211]]$seq_region_name
#> [1] "6"
#> 
#> [[211]]$strand
#> [1] 1
#> 
#> [[211]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[211]]$version
#> [1] 2
#> 
#> [[211]]$gene_id
#> [1] "ENSG00000235863"
#> 
#> [[211]]$start
#> [1] 33244917
#> 
#> [[211]]$canonical_transcript
#> [1] "ENST00000451237.1"
#> 
#> [[211]]$external_name
#> [1] "B3GALT4"
#> 
#> [[211]]$id
#> [1] "ENSG00000235863"
#> 
#> [[211]]$end
#> [1] 33252609
#> 
#> 
#> [[212]]
#> [[212]]$seq_region_name
#> [1] "6"
#> 
#> [[212]]$source
#> [1] "ensembl_havana"
#> 
#> [[212]]$description
#> [1] "WD repeat domain 46 [Source:HGNC Symbol;Acc:13923]"
#> 
#> [[212]]$feature_type
#> [1] "gene"
#> 
#> [[212]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[212]]$biotype
#> [1] "protein_coding"
#> 
#> [[212]]$start
#> [1] 33246885
#> 
#> [[212]]$id
#> [1] "ENSG00000227057"
#> 
#> [[212]]$end
#> [1] 33257304
#> 
#> [[212]]$canonical_transcript
#> [1] "ENST00000374617.4"
#> 
#> [[212]]$external_name
#> [1] "WDR46"
#> 
#> [[212]]$gene_id
#> [1] "ENSG00000227057"
#> 
#> [[212]]$strand
#> [1] -1
#> 
#> [[212]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[212]]$version
#> [1] 3
#> 
#> 
#> [[213]]
#> [[213]]$description
#> [1] "HCG43720, isoform CRA_a; Lymphocyte antigen 6 complex locus protein G6f; Megakaryocyte-enhanced gene transcript 1 protein; Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:Q9NZJ1]"
#> 
#> [[213]]$feature_type
#> [1] "gene"
#> 
#> [[213]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[213]]$biotype
#> [1] "protein_coding"
#> 
#> [[213]]$seq_region_name
#> [1] "6"
#> 
#> [[213]]$source
#> [1] "havana"
#> 
#> [[213]]$gene_id
#> [1] "ENSG00000250641"
#> 
#> [[213]]$strand
#> [1] 1
#> 
#> [[213]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[213]]$version
#> [1] 1
#> 
#> [[213]]$start
#> [1] 31674681
#> 
#> [[213]]$id
#> [1] "ENSG00000250641"
#> 
#> [[213]]$end
#> [1] 31685695
#> 
#> [[213]]$external_name
#> [1] "MEGT1"
#> 
#> [[213]]$canonical_transcript
#> [1] "ENST00000503322.1"
#> 
#> 
#> [[214]]
#> [[214]]$start
#> [1] 31679548
#> 
#> [[214]]$id
#> [1] "ENSG00000255552"
#> 
#> [[214]]$end
#> [1] 31681842
#> 
#> [[214]]$external_name
#> [1] "LY6G6E"
#> 
#> [[214]]$canonical_transcript
#> [1] "ENST00000383418.4"
#> 
#> [[214]]$gene_id
#> [1] "ENSG00000255552"
#> 
#> [[214]]$version
#> [1] 3
#> 
#> [[214]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[214]]$strand
#> [1] -1
#> 
#> [[214]]$seq_region_name
#> [1] "6"
#> 
#> [[214]]$source
#> [1] "ensembl_havana"
#> 
#> [[214]]$description
#> [1] "lymphocyte antigen 6 complex, locus G6E (pseudogene) [Source:HGNC Symbol;Acc:13934]"
#> 
#> [[214]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[214]]$feature_type
#> [1] "gene"
#> 
#> [[214]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[215]]
#> [[215]]$seq_region_name
#> [1] "6"
#> 
#> [[215]]$source
#> [1] "ensembl_havana"
#> 
#> [[215]]$feature_type
#> [1] "gene"
#> 
#> [[215]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[215]]$biotype
#> [1] "protein_coding"
#> 
#> [[215]]$description
#> [1] "lymphocyte antigen 6 complex, locus G6D [Source:HGNC Symbol;Acc:13935]"
#> 
#> [[215]]$id
#> [1] "ENSG00000244355"
#> 
#> [[215]]$end
#> [1] 31685581
#> 
#> [[215]]$external_name
#> [1] "LY6G6D"
#> 
#> [[215]]$canonical_transcript
#> [1] "ENST00000375825.3"
#> 
#> [[215]]$start
#> [1] 31683133
#> 
#> [[215]]$gene_id
#> [1] "ENSG00000244355"
#> 
#> [[215]]$version
#> [1] 3
#> 
#> [[215]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[215]]$strand
#> [1] 1
#> 
#> 
#> [[216]]
#> [[216]]$canonical_transcript
#> [1] "ENST00000375806.2"
#> 
#> [[216]]$external_name
#> [1] "C6orf25"
#> 
#> [[216]]$id
#> [1] "ENSG00000204420"
#> 
#> [[216]]$end
#> [1] 31694491
#> 
#> [[216]]$start
#> [1] 31686371
#> 
#> [[216]]$version
#> [1] 4
#> 
#> [[216]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[216]]$strand
#> [1] 1
#> 
#> [[216]]$gene_id
#> [1] "ENSG00000204420"
#> 
#> [[216]]$source
#> [1] "ensembl_havana"
#> 
#> [[216]]$seq_region_name
#> [1] "6"
#> 
#> [[216]]$biotype
#> [1] "protein_coding"
#> 
#> [[216]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[216]]$feature_type
#> [1] "gene"
#> 
#> [[216]]$description
#> [1] "chromosome 6 open reading frame 25 [Source:HGNC Symbol;Acc:13937]"
#> 
#> 
#> [[217]]
#> [[217]]$seq_region_name
#> [1] "6"
#> 
#> [[217]]$source
#> [1] "ensembl_havana"
#> 
#> [[217]]$description
#> [1] "prefoldin subunit 6 [Source:HGNC Symbol;Acc:4926]"
#> 
#> [[217]]$feature_type
#> [1] "gene"
#> 
#> [[217]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[217]]$biotype
#> [1] "protein_coding"
#> 
#> [[217]]$start
#> [1] 33257079
#> 
#> [[217]]$end
#> [1] 33266178
#> 
#> [[217]]$id
#> [1] "ENSG00000204220"
#> 
#> [[217]]$canonical_transcript
#> [1] "ENST00000395131.1"
#> 
#> [[217]]$external_name
#> [1] "PFDN6"
#> 
#> [[217]]$gene_id
#> [1] "ENSG00000204220"
#> 
#> [[217]]$strand
#> [1] 1
#> 
#> [[217]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[217]]$version
#> [1] 5
#> 
#> 
#> [[218]]
#> [[218]]$source
#> [1] "ensembl_havana"
#> 
#> [[218]]$seq_region_name
#> [1] "6"
#> 
#> [[218]]$biotype
#> [1] "protein_coding"
#> 
#> [[218]]$feature_type
#> [1] "gene"
#> 
#> [[218]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[218]]$description
#> [1] "ral guanine nucleotide dissociation stimulator-like 2 [Source:HGNC Symbol;Acc:9769]"
#> 
#> [[218]]$canonical_transcript
#> [1] "ENST00000497454.1"
#> 
#> [[218]]$external_name
#> [1] "RGL2"
#> 
#> [[218]]$id
#> [1] "ENSG00000237441"
#> 
#> [[218]]$end
#> [1] 33267101
#> 
#> [[218]]$start
#> [1] 33259431
#> 
#> [[218]]$strand
#> [1] -1
#> 
#> [[218]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[218]]$version
#> [1] 5
#> 
#> [[218]]$gene_id
#> [1] "ENSG00000237441"
#> 
#> 
#> [[219]]
#> [[219]]$seq_region_name
#> [1] "6"
#> 
#> [[219]]$source
#> [1] "ensembl_havana"
#> 
#> [[219]]$description
#> [1] "lymphocyte antigen 6 complex, locus G6C [Source:HGNC Symbol;Acc:13936]"
#> 
#> [[219]]$feature_type
#> [1] "gene"
#> 
#> [[219]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[219]]$biotype
#> [1] "protein_coding"
#> 
#> [[219]]$start
#> [1] 31686425
#> 
#> [[219]]$end
#> [1] 31689622
#> 
#> [[219]]$id
#> [1] "ENSG00000204421"
#> 
#> [[219]]$canonical_transcript
#> [1] "ENST00000375819.2"
#> 
#> [[219]]$external_name
#> [1] "LY6G6C"
#> 
#> [[219]]$gene_id
#> [1] "ENSG00000204421"
#> 
#> [[219]]$strand
#> [1] -1
#> 
#> [[219]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[219]]$version
#> [1] 2
#> 
#> 
#> [[220]]
#> [[220]]$canonical_transcript
#> [1] "ENST00000375789.2"
#> 
#> [[220]]$external_name
#> [1] "DDAH2"
#> 
#> [[220]]$end
#> [1] 31698394
#> 
#> [[220]]$id
#> [1] "ENSG00000213722"
#> 
#> [[220]]$start
#> [1] 31694815
#> 
#> [[220]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[220]]$strand
#> [1] -1
#> 
#> [[220]]$version
#> [1] 4
#> 
#> [[220]]$gene_id
#> [1] "ENSG00000213722"
#> 
#> [[220]]$source
#> [1] "ensembl_havana"
#> 
#> [[220]]$seq_region_name
#> [1] "6"
#> 
#> [[220]]$biotype
#> [1] "protein_coding"
#> 
#> [[220]]$feature_type
#> [1] "gene"
#> 
#> [[220]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[220]]$description
#> [1] "dimethylarginine dimethylaminohydrolase 2 [Source:HGNC Symbol;Acc:2716]"
#> 
#> 
#> [[221]]
#> [[221]]$start
#> [1] 31698358
#> 
#> [[221]]$canonical_transcript
#> [1] "ENST00000375780.2"
#> 
#> [[221]]$external_name
#> [1] "CLIC1"
#> 
#> [[221]]$id
#> [1] "ENSG00000213719"
#> 
#> [[221]]$end
#> [1] 31707540
#> 
#> [[221]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[221]]$strand
#> [1] -1
#> 
#> [[221]]$version
#> [1] 4
#> 
#> [[221]]$gene_id
#> [1] "ENSG00000213719"
#> 
#> [[221]]$source
#> [1] "ensembl_havana"
#> 
#> [[221]]$seq_region_name
#> [1] "6"
#> 
#> [[221]]$description
#> [1] "chloride intracellular channel 1 [Source:HGNC Symbol;Acc:2062]"
#> 
#> [[221]]$biotype
#> [1] "protein_coding"
#> 
#> [[221]]$feature_type
#> [1] "gene"
#> 
#> [[221]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> 
#> [[222]]
#> [[222]]$gene_id
#> [1] "ENSG00000231925"
#> 
#> [[222]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[222]]$strand
#> [1] -1
#> 
#> [[222]]$version
#> [1] 7
#> 
#> [[222]]$start
#> [1] 33267471
#> 
#> [[222]]$end
#> [1] 33282164
#> 
#> [[222]]$id
#> [1] "ENSG00000231925"
#> 
#> [[222]]$canonical_transcript
#> [1] "ENST00000426633.2"
#> 
#> [[222]]$external_name
#> [1] "TAPBP"
#> 
#> [[222]]$description
#> [1] "TAP binding protein (tapasin) [Source:HGNC Symbol;Acc:11566]"
#> 
#> [[222]]$feature_type
#> [1] "gene"
#> 
#> [[222]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[222]]$biotype
#> [1] "protein_coding"
#> 
#> [[222]]$seq_region_name
#> [1] "6"
#> 
#> [[222]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[223]]
#> [[223]]$seq_region_name
#> [1] "6"
#> 
#> [[223]]$source
#> [1] "ensembl_havana"
#> 
#> [[223]]$description
#> [1] "mutS homolog 5 [Source:HGNC Symbol;Acc:7328]"
#> 
#> [[223]]$feature_type
#> [1] "gene"
#> 
#> [[223]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[223]]$biotype
#> [1] "protein_coding"
#> 
#> [[223]]$start
#> [1] 31707725
#> 
#> [[223]]$end
#> [1] 31732622
#> 
#> [[223]]$id
#> [1] "ENSG00000204410"
#> 
#> [[223]]$canonical_transcript
#> [1] "ENST00000375703.3"
#> 
#> [[223]]$external_name
#> [1] "MSH5"
#> 
#> [[223]]$gene_id
#> [1] "ENSG00000204410"
#> 
#> [[223]]$strand
#> [1] 1
#> 
#> [[223]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[223]]$version
#> [1] 10
#> 
#> 
#> [[224]]
#> [[224]]$description
#> [1] "zinc finger and BTB domain containing 22 [Source:HGNC Symbol;Acc:13085]"
#> 
#> [[224]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[224]]$feature_type
#> [1] "gene"
#> 
#> [[224]]$biotype
#> [1] "protein_coding"
#> 
#> [[224]]$seq_region_name
#> [1] "6"
#> 
#> [[224]]$source
#> [1] "ensembl_havana"
#> 
#> [[224]]$gene_id
#> [1] "ENSG00000236104"
#> 
#> [[224]]$version
#> [1] 2
#> 
#> [[224]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[224]]$strand
#> [1] -1
#> 
#> [[224]]$start
#> [1] 33282183
#> 
#> [[224]]$id
#> [1] "ENSG00000236104"
#> 
#> [[224]]$end
#> [1] 33285719
#> 
#> [[224]]$external_name
#> [1] "ZBTB22"
#> 
#> [[224]]$canonical_transcript
#> [1] "ENST00000431845.2"
#> 
#> 
#> [[225]]
#> [[225]]$start
#> [1] 33286335
#> 
#> [[225]]$end
#> [1] 33297046
#> 
#> [[225]]$id
#> [1] "ENSG00000204209"
#> 
#> [[225]]$canonical_transcript
#> [1] "ENST00000374542.5"
#> 
#> [[225]]$external_name
#> [1] "DAXX"
#> 
#> [[225]]$gene_id
#> [1] "ENSG00000204209"
#> 
#> [[225]]$version
#> [1] 6
#> 
#> [[225]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[225]]$strand
#> [1] -1
#> 
#> [[225]]$seq_region_name
#> [1] "6"
#> 
#> [[225]]$source
#> [1] "ensembl_havana"
#> 
#> [[225]]$description
#> [1] "death-domain associated protein [Source:HGNC Symbol;Acc:2681]"
#> 
#> [[225]]$feature_type
#> [1] "gene"
#> 
#> [[225]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[225]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[226]]
#> [[226]]$seq_region_name
#> [1] "6"
#> 
#> [[226]]$source
#> [1] "ensembl_havana"
#> 
#> [[226]]$description
#> [1] "kinesin family member C1 [Source:HGNC Symbol;Acc:6389]"
#> 
#> [[226]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[226]]$feature_type
#> [1] "gene"
#> 
#> [[226]]$biotype
#> [1] "protein_coding"
#> 
#> [[226]]$start
#> [1] 33359313
#> 
#> [[226]]$end
#> [1] 33377701
#> 
#> [[226]]$id
#> [1] "ENSG00000237649"
#> 
#> [[226]]$external_name
#> [1] "KIFC1"
#> 
#> [[226]]$canonical_transcript
#> [1] "ENST00000428849.2"
#> 
#> [[226]]$gene_id
#> [1] "ENSG00000237649"
#> 
#> [[226]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[226]]$version
#> [1] 3
#> 
#> [[226]]$strand
#> [1] 1
#> 
#> 
#> [[227]]
#> [[227]]$seq_region_name
#> [1] "6"
#> 
#> [[227]]$source
#> [1] "ensembl_havana"
#> 
#> [[227]]$feature_type
#> [1] "gene"
#> 
#> [[227]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[227]]$biotype
#> [1] "protein_coding"
#> 
#> [[227]]$description
#> [1] "PHD finger protein 1 [Source:HGNC Symbol;Acc:8919]"
#> 
#> [[227]]$end
#> [1] 33384230
#> 
#> [[227]]$id
#> [1] "ENSG00000112511"
#> 
#> [[227]]$canonical_transcript
#> [1] "ENST00000374516.3"
#> 
#> [[227]]$external_name
#> [1] "PHF1"
#> 
#> [[227]]$start
#> [1] 33378176
#> 
#> [[227]]$gene_id
#> [1] "ENSG00000112511"
#> 
#> [[227]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[227]]$version
#> [1] 13
#> 
#> [[227]]$strand
#> [1] 1
#> 
#> 
#> [[228]]
#> [[228]]$gene_id
#> [1] "ENSG00000255152"
#> 
#> [[228]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[228]]$strand
#> [1] 1
#> 
#> [[228]]$version
#> [1] 4
#> 
#> [[228]]$end
#> [1] 31732628
#> 
#> [[228]]$id
#> [1] "ENSG00000255152"
#> 
#> [[228]]$canonical_transcript
#> [1] "ENST00000493662.2"
#> 
#> [[228]]$external_name
#> [1] "MSH5-SAPCD1"
#> 
#> [[228]]$start
#> [1] 31707797
#> 
#> [[228]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[228]]$feature_type
#> [1] "gene"
#> 
#> [[228]]$biotype
#> [1] "protein_coding"
#> 
#> [[228]]$description
#> [1] "MSH5-SAPCD1 readthrough (NMD candidate) [Source:HGNC Symbol;Acc:41994]"
#> 
#> [[228]]$seq_region_name
#> [1] "6"
#> 
#> [[228]]$source
#> [1] "havana"
#> 
#> 
#> [[229]]
#> [[229]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[229]]$version
#> [1] 11
#> 
#> [[229]]$strand
#> [1] -1
#> 
#> [[229]]$gene_id
#> [1] "ENSG00000112514"
#> 
#> [[229]]$external_name
#> [1] "CUTA"
#> 
#> [[229]]$canonical_transcript
#> [1] "ENST00000374500.5"
#> 
#> [[229]]$end
#> [1] 33386094
#> 
#> [[229]]$id
#> [1] "ENSG00000112514"
#> 
#> [[229]]$start
#> [1] 33384219
#> 
#> [[229]]$biotype
#> [1] "protein_coding"
#> 
#> [[229]]$feature_type
#> [1] "gene"
#> 
#> [[229]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[229]]$description
#> [1] "cutA divalent cation tolerance homolog (E. coli) [Source:HGNC Symbol;Acc:21101]"
#> 
#> [[229]]$source
#> [1] "ensembl_havana"
#> 
#> [[229]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[230]]
#> [[230]]$gene_id
#> [1] "ENSG00000228727"
#> 
#> [[230]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[230]]$version
#> [1] 4
#> 
#> [[230]]$strand
#> [1] 1
#> 
#> [[230]]$id
#> [1] "ENSG00000228727"
#> 
#> [[230]]$end
#> [1] 31732628
#> 
#> [[230]]$external_name
#> [1] "SAPCD1"
#> 
#> [[230]]$canonical_transcript
#> [1] "ENST00000415669.2"
#> 
#> [[230]]$start
#> [1] 31730576
#> 
#> [[230]]$feature_type
#> [1] "gene"
#> 
#> [[230]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[230]]$biotype
#> [1] "protein_coding"
#> 
#> [[230]]$description
#> [1] "suppressor APC domain containing 1 [Source:HGNC Symbol;Acc:13938]"
#> 
#> [[230]]$seq_region_name
#> [1] "6"
#> 
#> [[230]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[231]]
#> [[231]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[231]]$feature_type
#> [1] "gene"
#> 
#> [[231]]$biotype
#> [1] "protein_coding"
#> 
#> [[231]]$description
#> [1] "von Willebrand factor A domain containing 7 [Source:HGNC Symbol;Acc:13939]"
#> 
#> [[231]]$seq_region_name
#> [1] "6"
#> 
#> [[231]]$source
#> [1] "ensembl_havana"
#> 
#> [[231]]$gene_id
#> [1] "ENSG00000204396"
#> 
#> [[231]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[231]]$version
#> [1] 6
#> 
#> [[231]]$strand
#> [1] -1
#> 
#> [[231]]$end
#> [1] 31745108
#> 
#> [[231]]$id
#> [1] "ENSG00000204396"
#> 
#> [[231]]$external_name
#> [1] "VWA7"
#> 
#> [[231]]$canonical_transcript
#> [1] "ENST00000375688.4"
#> 
#> [[231]]$start
#> [1] 31733367
#> 
#> 
#> [[232]]
#> [[232]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[232]]$feature_type
#> [1] "gene"
#> 
#> [[232]]$biotype
#> [1] "protein_coding"
#> 
#> [[232]]$description
#> [1] "valyl-tRNA synthetase [Source:HGNC Symbol;Acc:12651]"
#> 
#> [[232]]$seq_region_name
#> [1] "6"
#> 
#> [[232]]$source
#> [1] "ensembl_havana"
#> 
#> [[232]]$gene_id
#> [1] "ENSG00000204394"
#> 
#> [[232]]$version
#> [1] 8
#> 
#> [[232]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[232]]$strand
#> [1] -1
#> 
#> [[232]]$id
#> [1] "ENSG00000204394"
#> 
#> [[232]]$end
#> [1] 31763730
#> 
#> [[232]]$canonical_transcript
#> [1] "ENST00000375663.3"
#> 
#> [[232]]$external_name
#> [1] "VARS"
#> 
#> [[232]]$start
#> [1] 31745295
#> 
#> 
#> [[233]]
#> [[233]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[233]]$strand
#> [1] 1
#> 
#> [[233]]$version
#> [1] 8
#> 
#> [[233]]$gene_id
#> [1] "ENSG00000197283"
#> 
#> [[233]]$canonical_transcript
#> [1] "ENST00000418600.2"
#> 
#> [[233]]$external_name
#> [1] "SYNGAP1"
#> 
#> [[233]]$id
#> [1] "ENSG00000197283"
#> 
#> [[233]]$end
#> [1] 33421466
#> 
#> [[233]]$start
#> [1] 33387847
#> 
#> [[233]]$biotype
#> [1] "protein_coding"
#> 
#> [[233]]$feature_type
#> [1] "gene"
#> 
#> [[233]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[233]]$description
#> [1] "synaptic Ras GTPase activating protein 1 [Source:HGNC Symbol;Acc:11497]"
#> 
#> [[233]]$source
#> [1] "ensembl_havana"
#> 
#> [[233]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[234]]
#> [[234]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[234]]$strand
#> [1] -1
#> 
#> [[234]]$version
#> [1] 6
#> 
#> [[234]]$gene_id
#> [1] "ENSG00000204392"
#> 
#> [[234]]$external_name
#> [1] "LSM2"
#> 
#> [[234]]$canonical_transcript
#> [1] "ENST00000375661.5"
#> 
#> [[234]]$id
#> [1] "ENSG00000204392"
#> 
#> [[234]]$end
#> [1] 31774761
#> 
#> [[234]]$start
#> [1] 31765173
#> 
#> [[234]]$biotype
#> [1] "protein_coding"
#> 
#> [[234]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[234]]$feature_type
#> [1] "gene"
#> 
#> [[234]]$description
#> [1] "LSM2 homolog, U6 small nuclear RNA associated (S. cerevisiae) [Source:HGNC Symbol;Acc:13940]"
#> 
#> [[234]]$source
#> [1] "ensembl_havana"
#> 
#> [[234]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[235]]
#> [[235]]$seq_region_name
#> [1] "6"
#> 
#> [[235]]$source
#> [1] "ensembl_havana"
#> 
#> [[235]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[235]]$feature_type
#> [1] "gene"
#> 
#> [[235]]$biotype
#> [1] "protein_coding"
#> 
#> [[235]]$description
#> [1] "zinc finger and BTB domain containing 9 [Source:HGNC Symbol;Acc:28323]"
#> 
#> [[235]]$id
#> [1] "ENSG00000213588"
#> 
#> [[235]]$end
#> [1] 33425325
#> 
#> [[235]]$canonical_transcript
#> [1] "ENST00000395064.2"
#> 
#> [[235]]$external_name
#> [1] "ZBTB9"
#> 
#> [[235]]$start
#> [1] 33422356
#> 
#> [[235]]$gene_id
#> [1] "ENSG00000213588"
#> 
#> [[235]]$version
#> [1] 4
#> 
#> [[235]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[235]]$strand
#> [1] 1
#> 
#> 
#> [[236]]
#> [[236]]$strand
#> [1] -1
#> 
#> [[236]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[236]]$version
#> [1] 5
#> 
#> [[236]]$gene_id
#> [1] "ENSG00000198704"
#> 
#> [[236]]$canonical_transcript
#> [1] "ENST00000361902.1"
#> 
#> [[236]]$external_name
#> [1] "GPX6"
#> 
#> [[236]]$id
#> [1] "ENSG00000198704"
#> 
#> [[236]]$end
#> [1] 28495992
#> 
#> [[236]]$start
#> [1] 28471073
#> 
#> [[236]]$biotype
#> [1] "protein_coding"
#> 
#> [[236]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[236]]$feature_type
#> [1] "gene"
#> 
#> [[236]]$description
#> [1] "glutathione peroxidase 6 (olfactory) [Source:HGNC Symbol;Acc:4558]"
#> 
#> [[236]]$source
#> [1] "ensembl_havana"
#> 
#> [[236]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[237]]
#> [[237]]$version
#> [1] 8
#> 
#> [[237]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[237]]$strand
#> [1] -1
#> 
#> [[237]]$gene_id
#> [1] "ENSG00000204390"
#> 
#> [[237]]$start
#> [1] 31777396
#> 
#> [[237]]$canonical_transcript
#> [1] "ENST00000375654.4"
#> 
#> [[237]]$external_name
#> [1] "HSPA1L"
#> 
#> [[237]]$id
#> [1] "ENSG00000204390"
#> 
#> [[237]]$end
#> [1] 31783437
#> 
#> [[237]]$description
#> [1] "heat shock 70kDa protein 1-like [Source:HGNC Symbol;Acc:5234]"
#> 
#> [[237]]$biotype
#> [1] "protein_coding"
#> 
#> [[237]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[237]]$feature_type
#> [1] "gene"
#> 
#> [[237]]$source
#> [1] "ensembl_havana"
#> 
#> [[237]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[238]]
#> [[238]]$strand
#> [1] 1
#> 
#> [[238]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[238]]$version
#> [1] 8
#> 
#> [[238]]$gene_id
#> [1] "ENSG00000204389"
#> 
#> [[238]]$external_name
#> [1] "HSPA1A"
#> 
#> [[238]]$canonical_transcript
#> [1] "ENST00000375651.5"
#> 
#> [[238]]$end
#> [1] 31785723
#> 
#> [[238]]$id
#> [1] "ENSG00000204389"
#> 
#> [[238]]$start
#> [1] 31783291
#> 
#> [[238]]$biotype
#> [1] "protein_coding"
#> 
#> [[238]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[238]]$feature_type
#> [1] "gene"
#> 
#> [[238]]$description
#> [1] "heat shock 70kDa protein 1A [Source:HGNC Symbol;Acc:5232]"
#> 
#> [[238]]$source
#> [1] "ensembl_havana"
#> 
#> [[238]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[239]]
#> [[239]]$description
#> [1] "glutathione peroxidase 5 (epididymal androgen-related protein) [Source:HGNC Symbol;Acc:4557]"
#> 
#> [[239]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[239]]$feature_type
#> [1] "gene"
#> 
#> [[239]]$biotype
#> [1] "protein_coding"
#> 
#> [[239]]$seq_region_name
#> [1] "6"
#> 
#> [[239]]$source
#> [1] "ensembl_havana"
#> 
#> [[239]]$gene_id
#> [1] "ENSG00000224586"
#> 
#> [[239]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[239]]$version
#> [1] 2
#> 
#> [[239]]$strand
#> [1] 1
#> 
#> [[239]]$start
#> [1] 28493702
#> 
#> [[239]]$end
#> [1] 28502729
#> 
#> [[239]]$id
#> [1] "ENSG00000224586"
#> 
#> [[239]]$external_name
#> [1] "GPX5"
#> 
#> [[239]]$canonical_transcript
#> [1] "ENST00000412168.2"
#> 
#> 
#> [[240]]
#> [[240]]$version
#> [1] 2
#> 
#> [[240]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[240]]$strand
#> [1] -1
#> 
#> [[240]]$gene_id
#> [1] "ENSG00000232040"
#> 
#> [[240]]$external_name
#> [1] "SCAND3"
#> 
#> [[240]]$canonical_transcript
#> [1] "ENST00000452236.2"
#> 
#> [[240]]$end
#> [1] 28583989
#> 
#> [[240]]$id
#> [1] "ENSG00000232040"
#> 
#> [[240]]$start
#> [1] 28539407
#> 
#> [[240]]$biotype
#> [1] "protein_coding"
#> 
#> [[240]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[240]]$feature_type
#> [1] "gene"
#> 
#> [[240]]$description
#> [1] "SCAN domain containing 3 [Source:HGNC Symbol;Acc:13851]"
#> 
#> [[240]]$source
#> [1] "ensembl_havana"
#> 
#> [[240]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[241]]
#> [[241]]$description
#> [1] "heat shock 70kDa protein 1B [Source:HGNC Symbol;Acc:5233]"
#> 
#> [[241]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[241]]$feature_type
#> [1] "gene"
#> 
#> [[241]]$biotype
#> [1] "protein_coding"
#> 
#> [[241]]$seq_region_name
#> [1] "6"
#> 
#> [[241]]$source
#> [1] "ensembl_havana"
#> 
#> [[241]]$gene_id
#> [1] "ENSG00000204388"
#> 
#> [[241]]$version
#> [1] 5
#> 
#> [[241]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[241]]$strand
#> [1] 1
#> 
#> [[241]]$start
#> [1] 31795512
#> 
#> [[241]]$end
#> [1] 31798031
#> 
#> [[241]]$id
#> [1] "ENSG00000204388"
#> 
#> [[241]]$external_name
#> [1] "HSPA1B"
#> 
#> [[241]]$canonical_transcript
#> [1] "ENST00000375650.3"
#> 
#> 
#> [[242]]
#> [[242]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[242]]$version
#> [1] 8
#> 
#> [[242]]$strand
#> [1] 1
#> 
#> [[242]]$gene_id
#> [1] "ENSG00000204387"
#> 
#> [[242]]$start
#> [1] 31802385
#> 
#> [[242]]$external_name
#> [1] "C6orf48"
#> 
#> [[242]]$canonical_transcript
#> [1] "ENST00000375640.3"
#> 
#> [[242]]$id
#> [1] "ENSG00000204387"
#> 
#> [[242]]$end
#> [1] 31807541
#> 
#> [[242]]$description
#> [1] "chromosome 6 open reading frame 48 [Source:HGNC Symbol;Acc:19078]"
#> 
#> [[242]]$biotype
#> [1] "protein_coding"
#> 
#> [[242]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[242]]$feature_type
#> [1] "gene"
#> 
#> [[242]]$source
#> [1] "ensembl_havana"
#> 
#> [[242]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[243]]
#> [[243]]$description
#> [1] "sialidase 1 (lysosomal sialidase) [Source:HGNC Symbol;Acc:7758]"
#> 
#> [[243]]$biotype
#> [1] "protein_coding"
#> 
#> [[243]]$feature_type
#> [1] "gene"
#> 
#> [[243]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[243]]$source
#> [1] "ensembl_havana"
#> 
#> [[243]]$seq_region_name
#> [1] "6"
#> 
#> [[243]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[243]]$version
#> [1] 6
#> 
#> [[243]]$strand
#> [1] -1
#> 
#> [[243]]$gene_id
#> [1] "ENSG00000204386"
#> 
#> [[243]]$start
#> [1] 31825436
#> 
#> [[243]]$canonical_transcript
#> [1] "ENST00000375631.4"
#> 
#> [[243]]$external_name
#> [1] "NEU1"
#> 
#> [[243]]$id
#> [1] "ENSG00000204386"
#> 
#> [[243]]$end
#> [1] 31830683
#> 
#> 
#> [[244]]
#> [[244]]$source
#> [1] "ensembl_havana"
#> 
#> [[244]]$seq_region_name
#> [1] "6"
#> 
#> [[244]]$biotype
#> [1] "protein_coding"
#> 
#> [[244]]$feature_type
#> [1] "gene"
#> 
#> [[244]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[244]]$description
#> [1] "tripartite motif containing 27 [Source:HGNC Symbol;Acc:9975]"
#> 
#> [[244]]$canonical_transcript
#> [1] "ENST00000377199.3"
#> 
#> [[244]]$external_name
#> [1] "TRIM27"
#> 
#> [[244]]$end
#> [1] 28891766
#> 
#> [[244]]$id
#> [1] "ENSG00000204713"
#> 
#> [[244]]$start
#> [1] 28870779
#> 
#> [[244]]$version
#> [1] 6
#> 
#> [[244]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[244]]$strand
#> [1] -1
#> 
#> [[244]]$gene_id
#> [1] "ENSG00000204713"
#> 
#> 
#> [[245]]
#> [[245]]$source
#> [1] "havana"
#> 
#> [[245]]$seq_region_name
#> [1] "6"
#> 
#> [[245]]$description
#> [1] "chromosome 6 open reading frame 100 [Source:HGNC Symbol;Acc:21195]"
#> 
#> [[245]]$biotype
#> [1] "protein_coding"
#> 
#> [[245]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[245]]$feature_type
#> [1] "gene"
#> 
#> [[245]]$start
#> [1] 28911654
#> 
#> [[245]]$external_name
#> [1] "C6orf100"
#> 
#> [[245]]$canonical_transcript
#> [1] "ENST00000377186.3"
#> 
#> [[245]]$id
#> [1] "ENSG00000204709"
#> 
#> [[245]]$end
#> [1] 28912314
#> 
#> [[245]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[245]]$version
#> [1] 4
#> 
#> [[245]]$strand
#> [1] 1
#> 
#> [[245]]$gene_id
#> [1] "ENSG00000204709"
#> 
#> 
#> [[246]]
#> [[246]]$start
#> [1] 28962562
#> 
#> [[246]]$external_name
#> [1] "ZNF311"
#> 
#> [[246]]$canonical_transcript
#> [1] "ENST00000377179.3"
#> 
#> [[246]]$end
#> [1] 28973093
#> 
#> [[246]]$id
#> [1] "ENSG00000197935"
#> 
#> [[246]]$version
#> [1] 6
#> 
#> [[246]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[246]]$strand
#> [1] -1
#> 
#> [[246]]$gene_id
#> [1] "ENSG00000197935"
#> 
#> [[246]]$source
#> [1] "ensembl_havana"
#> 
#> [[246]]$seq_region_name
#> [1] "6"
#> 
#> [[246]]$description
#> [1] "zinc finger protein 311 [Source:HGNC Symbol;Acc:13847]"
#> 
#> [[246]]$biotype
#> [1] "protein_coding"
#> 
#> [[246]]$feature_type
#> [1] "gene"
#> 
#> [[246]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> 
#> [[247]]
#> [[247]]$source
#> [1] "ensembl_havana"
#> 
#> [[247]]$seq_region_name
#> [1] "6"
#> 
#> [[247]]$description
#> [1] "olfactory receptor, family 2, subfamily W, member 1 [Source:HGNC Symbol;Acc:8281]"
#> 
#> [[247]]$biotype
#> [1] "protein_coding"
#> 
#> [[247]]$feature_type
#> [1] "gene"
#> 
#> [[247]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[247]]$start
#> [1] 29011990
#> 
#> [[247]]$external_name
#> [1] "OR2W1"
#> 
#> [[247]]$canonical_transcript
#> [1] "ENST00000377175.1"
#> 
#> [[247]]$id
#> [1] "ENSG00000204704"
#> 
#> [[247]]$end
#> [1] 29013017
#> 
#> [[247]]$strand
#> [1] -1
#> 
#> [[247]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[247]]$version
#> [1] 2
#> 
#> [[247]]$gene_id
#> [1] "ENSG00000204704"
#> 
#> 
#> [[248]]
#> [[248]]$description
#> [1] "olfactory receptor, family 2, subfamily B, member 3 [Source:HGNC Symbol;Acc:8238]"
#> 
#> [[248]]$biotype
#> [1] "protein_coding"
#> 
#> [[248]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[248]]$feature_type
#> [1] "gene"
#> 
#> [[248]]$source
#> [1] "ensembl_havana"
#> 
#> [[248]]$seq_region_name
#> [1] "6"
#> 
#> [[248]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[248]]$strand
#> [1] -1
#> 
#> [[248]]$version
#> [1] 3
#> 
#> [[248]]$gene_id
#> [1] "ENSG00000204703"
#> 
#> [[248]]$start
#> [1] 29053985
#> 
#> [[248]]$canonical_transcript
#> [1] "ENST00000377173.2"
#> 
#> [[248]]$external_name
#> [1] "OR2B3"
#> 
#> [[248]]$end
#> [1] 29055090
#> 
#> [[248]]$id
#> [1] "ENSG00000204703"
#> 
#> 
#> [[249]]
#> [[249]]$description
#> [1] "olfactory receptor, family 2, subfamily J, member 1 (gene/pseudogene) [Source:HGNC Symbol;Acc:8259]"
#> 
#> [[249]]$biotype
#> [1] "protein_coding"
#> 
#> [[249]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[249]]$feature_type
#> [1] "gene"
#> 
#> [[249]]$source
#> [1] "havana"
#> 
#> [[249]]$seq_region_name
#> [1] "6"
#> 
#> [[249]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[249]]$version
#> [1] 4
#> 
#> [[249]]$strand
#> [1] 1
#> 
#> [[249]]$gene_id
#> [1] "ENSG00000204702"
#> 
#> [[249]]$start
#> [1] 29068386
#> 
#> [[249]]$canonical_transcript
#> [1] "ENST00000377171.3"
#> 
#> [[249]]$external_name
#> [1] "OR2J1"
#> 
#> [[249]]$id
#> [1] "ENSG00000204702"
#> 
#> [[249]]$end
#> [1] 29069658
#> 
#> 
#> [[250]]
#> [[250]]$seq_region_name
#> [1] "6"
#> 
#> [[250]]$source
#> [1] "ensembl_havana"
#> 
#> [[250]]$description
#> [1] "olfactory receptor, family 2, subfamily J, member 3 [Source:HGNC Symbol;Acc:8261]"
#> 
#> [[250]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[250]]$feature_type
#> [1] "gene"
#> 
#> [[250]]$biotype
#> [1] "protein_coding"
#> 
#> [[250]]$start
#> [1] 29079668
#> 
#> [[250]]$end
#> [1] 29080603
#> 
#> [[250]]$id
#> [1] "ENSG00000204701"
#> 
#> [[250]]$canonical_transcript
#> [1] "ENST00000377169.1"
#> 
#> [[250]]$external_name
#> [1] "OR2J3"
#> 
#> [[250]]$gene_id
#> [1] "ENSG00000204701"
#> 
#> [[250]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[250]]$version
#> [1] 1
#> 
#> [[250]]$strand
#> [1] 1
#> 
#> 
#> [[251]]
#> [[251]]$external_name
#> [1] "SLC44A4"
#> 
#> [[251]]$canonical_transcript
#> [1] "ENST00000229729.6"
#> 
#> [[251]]$end
#> [1] 31846823
#> 
#> [[251]]$id
#> [1] "ENSG00000204385"
#> 
#> [[251]]$start
#> [1] 31830969
#> 
#> [[251]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[251]]$strand
#> [1] -1
#> 
#> [[251]]$version
#> [1] 6
#> 
#> [[251]]$gene_id
#> [1] "ENSG00000204385"
#> 
#> [[251]]$source
#> [1] "ensembl_havana"
#> 
#> [[251]]$seq_region_name
#> [1] "6"
#> 
#> [[251]]$biotype
#> [1] "protein_coding"
#> 
#> [[251]]$feature_type
#> [1] "gene"
#> 
#> [[251]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[251]]$description
#> [1] "solute carrier family 44, member 4 [Source:HGNC Symbol;Acc:13941]"
#> 
#> 
#> [[252]]
#> [[252]]$end
#> [1] 29142351
#> 
#> [[252]]$id
#> [1] "ENSG00000204700"
#> 
#> [[252]]$canonical_transcript
#> [1] "ENST00000377167.2"
#> 
#> [[252]]$external_name
#> [1] "OR2J2"
#> 
#> [[252]]$start
#> [1] 29141311
#> 
#> [[252]]$gene_id
#> [1] "ENSG00000204700"
#> 
#> [[252]]$strand
#> [1] 1
#> 
#> [[252]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[252]]$version
#> [1] 3
#> 
#> [[252]]$seq_region_name
#> [1] "6"
#> 
#> [[252]]$source
#> [1] "ensembl_havana"
#> 
#> [[252]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[252]]$feature_type
#> [1] "gene"
#> 
#> [[252]]$biotype
#> [1] "protein_coding"
#> 
#> [[252]]$description
#> [1] "olfactory receptor, family 2, subfamily J, member 2 [Source:HGNC Symbol;Acc:8260]"
#> 
#> 
#> [[253]]
#> [[253]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[253]]$feature_type
#> [1] "gene"
#> 
#> [[253]]$biotype
#> [1] "protein_coding"
#> 
#> [[253]]$description
#> [1] "olfactory receptor, family 14, subfamily J, member 1 [Source:HGNC Symbol;Acc:13971]"
#> 
#> [[253]]$seq_region_name
#> [1] "6"
#> 
#> [[253]]$source
#> [1] "ensembl_havana"
#> 
#> [[253]]$gene_id
#> [1] "ENSG00000204695"
#> 
#> [[253]]$version
#> [1] 2
#> 
#> [[253]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[253]]$strand
#> [1] 1
#> 
#> [[253]]$end
#> [1] 29275519
#> 
#> [[253]]$id
#> [1] "ENSG00000204695"
#> 
#> [[253]]$canonical_transcript
#> [1] "ENST00000377160.2"
#> 
#> [[253]]$external_name
#> [1] "OR14J1"
#> 
#> [[253]]$start
#> [1] 29274403
#> 
#> 
#> [[254]]
#> [[254]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[254]]$strand
#> [1] -1
#> 
#> [[254]]$version
#> [1] 3
#> 
#> [[254]]$gene_id
#> [1] "ENSG00000243729"
#> 
#> [[254]]$start
#> [1] 29323007
#> 
#> [[254]]$external_name
#> [1] "OR5V1"
#> 
#> [[254]]$canonical_transcript
#> [1] "ENST00000377154.1"
#> 
#> [[254]]$end
#> [1] 29399744
#> 
#> [[254]]$id
#> [1] "ENSG00000243729"
#> 
#> [[254]]$description
#> [1] "olfactory receptor, family 5, subfamily V, member 1 [Source:HGNC Symbol;Acc:13972]"
#> 
#> [[254]]$biotype
#> [1] "protein_coding"
#> 
#> [[254]]$feature_type
#> [1] "gene"
#> 
#> [[254]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[254]]$source
#> [1] "ensembl_havana"
#> 
#> [[254]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[255]]
#> [[255]]$start
#> [1] 31847536
#> 
#> [[255]]$id
#> [1] "ENSG00000204371"
#> 
#> [[255]]$end
#> [1] 31865464
#> 
#> [[255]]$external_name
#> [1] "EHMT2"
#> 
#> [[255]]$canonical_transcript
#> [1] "ENST00000375537.4"
#> 
#> [[255]]$gene_id
#> [1] "ENSG00000204371"
#> 
#> [[255]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[255]]$strand
#> [1] -1
#> 
#> [[255]]$version
#> [1] 7
#> 
#> [[255]]$seq_region_name
#> [1] "6"
#> 
#> [[255]]$source
#> [1] "ensembl_havana"
#> 
#> [[255]]$description
#> [1] "euchromatic histone-lysine N-methyltransferase 2 [Source:HGNC Symbol;Acc:14129]"
#> 
#> [[255]]$feature_type
#> [1] "gene"
#> 
#> [[255]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[255]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[256]]
#> [[256]]$gene_id
#> [1] "ENSG00000112462"
#> 
#> [[256]]$strand
#> [1] -1
#> 
#> [[256]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[256]]$version
#> [1] 8
#> 
#> [[256]]$id
#> [1] "ENSG00000112462"
#> 
#> [[256]]$end
#> [1] 29343068
#> 
#> [[256]]$external_name
#> [1] "OR12D3"
#> 
#> [[256]]$canonical_transcript
#> [1] "ENST00000396806.3"
#> 
#> [[256]]$start
#> [1] 29341200
#> 
#> [[256]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[256]]$feature_type
#> [1] "gene"
#> 
#> [[256]]$biotype
#> [1] "protein_coding"
#> 
#> [[256]]$description
#> [1] "olfactory receptor, family 12, subfamily D, member 3 [Source:HGNC Symbol;Acc:13963]"
#> 
#> [[256]]$seq_region_name
#> [1] "6"
#> 
#> [[256]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[257]]
#> [[257]]$start
#> [1] 29364416
#> 
#> [[257]]$external_name
#> [1] "OR12D2"
#> 
#> [[257]]$canonical_transcript
#> [1] "ENST00000383555.2"
#> 
#> [[257]]$end
#> [1] 29365448
#> 
#> [[257]]$id
#> [1] "ENSG00000168787"
#> 
#> [[257]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[257]]$version
#> [1] 4
#> 
#> [[257]]$strand
#> [1] 1
#> 
#> [[257]]$gene_id
#> [1] "ENSG00000168787"
#> 
#> [[257]]$source
#> [1] "ensembl_havana"
#> 
#> [[257]]$seq_region_name
#> [1] "6"
#> 
#> [[257]]$description
#> [1] "olfactory receptor, family 12, subfamily D, member 2 (gene/pseudogene) [Source:HGNC Symbol;Acc:8178]"
#> 
#> [[257]]$biotype
#> [1] "protein_coding"
#> 
#> [[257]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[257]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[258]]
#> [[258]]$start
#> [1] 29385057
#> 
#> [[258]]$external_name
#> [1] "OR12D1"
#> 
#> [[258]]$canonical_transcript
#> [1] "ENST00000514827.1"
#> 
#> [[258]]$id
#> [1] "ENSG00000251608"
#> 
#> [[258]]$end
#> [1] 29386003
#> 
#> [[258]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[258]]$strand
#> [1] 1
#> 
#> [[258]]$version
#> [1] 1
#> 
#> [[258]]$gene_id
#> [1] "ENSG00000251608"
#> 
#> [[258]]$source
#> [1] "havana"
#> 
#> [[258]]$seq_region_name
#> [1] "6"
#> 
#> [[258]]$description
#> [1] "olfactory receptor, family 12, subfamily D, member 1 (gene/pseudogene) [Source:HGNC Symbol;Acc:8177]"
#> 
#> [[258]]$biotype
#> [1] "polymorphic_pseudogene"
#> 
#> [[258]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[258]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[259]]
#> [[259]]$gene_id
#> [1] "ENSG00000204694"
#> 
#> [[259]]$strand
#> [1] -1
#> 
#> [[259]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[259]]$version
#> [1] 6
#> 
#> [[259]]$id
#> [1] "ENSG00000204694"
#> 
#> [[259]]$end
#> [1] 29424848
#> 
#> [[259]]$canonical_transcript
#> [1] "ENST00000377149.1"
#> 
#> [[259]]$external_name
#> [1] "OR11A1"
#> 
#> [[259]]$start
#> [1] 29393281
#> 
#> [[259]]$feature_type
#> [1] "gene"
#> 
#> [[259]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[259]]$biotype
#> [1] "protein_coding"
#> 
#> [[259]]$description
#> [1] "olfactory receptor, family 11, subfamily A, member 1 [Source:HGNC Symbol;Acc:8176]"
#> 
#> [[259]]$seq_region_name
#> [1] "6"
#> 
#> [[259]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[260]]
#> [[260]]$external_name
#> [1] "OR10C1"
#> 
#> [[260]]$canonical_transcript
#> [1] "ENST00000444197.2"
#> 
#> [[260]]$end
#> [1] 29408731
#> 
#> [[260]]$id
#> [1] "ENSG00000206474"
#> 
#> [[260]]$start
#> [1] 29407083
#> 
#> [[260]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[260]]$version
#> [1] 6
#> 
#> [[260]]$strand
#> [1] 1
#> 
#> [[260]]$gene_id
#> [1] "ENSG00000206474"
#> 
#> [[260]]$source
#> [1] "ensembl_havana"
#> 
#> [[260]]$seq_region_name
#> [1] "6"
#> 
#> [[260]]$biotype
#> [1] "protein_coding"
#> 
#> [[260]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[260]]$feature_type
#> [1] "gene"
#> 
#> [[260]]$description
#> [1] "olfactory receptor, family 10, subfamily C, member 1 (gene/pseudogene) [Source:HGNC Symbol;Acc:8165]"
#> 
#> 
#> [[261]]
#> [[261]]$start
#> [1] 29424958
#> 
#> [[261]]$canonical_transcript
#> [1] "ENST00000377136.1"
#> 
#> [[261]]$external_name
#> [1] "OR2H1"
#> 
#> [[261]]$end
#> [1] 29432105
#> 
#> [[261]]$id
#> [1] "ENSG00000204688"
#> 
#> [[261]]$strand
#> [1] 1
#> 
#> [[261]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[261]]$version
#> [1] 5
#> 
#> [[261]]$gene_id
#> [1] "ENSG00000204688"
#> 
#> [[261]]$source
#> [1] "ensembl_havana"
#> 
#> [[261]]$seq_region_name
#> [1] "6"
#> 
#> [[261]]$description
#> [1] "olfactory receptor, family 2, subfamily H, member 1 [Source:HGNC Symbol;Acc:8252]"
#> 
#> [[261]]$biotype
#> [1] "protein_coding"
#> 
#> [[261]]$feature_type
#> [1] "gene"
#> 
#> [[261]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> 
#> [[262]]
#> [[262]]$gene_id
#> [1] "ENSG00000204687"
#> 
#> [[262]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[262]]$strand
#> [1] -1
#> 
#> [[262]]$version
#> [1] 3
#> 
#> [[262]]$start
#> [1] 29454474
#> 
#> [[262]]$end
#> [1] 29455738
#> 
#> [[262]]$id
#> [1] "ENSG00000204687"
#> 
#> [[262]]$canonical_transcript
#> [1] "ENST00000377127.3"
#> 
#> [[262]]$external_name
#> [1] "MAS1L"
#> 
#> [[262]]$description
#> [1] "MAS1 oncogene-like [Source:HGNC Symbol;Acc:13961]"
#> 
#> [[262]]$feature_type
#> [1] "gene"
#> 
#> [[262]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[262]]$biotype
#> [1] "protein_coding"
#> 
#> [[262]]$seq_region_name
#> [1] "6"
#> 
#> [[262]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[263]]
#> [[263]]$gene_id
#> [1] "ENSG00000213886"
#> 
#> [[263]]$version
#> [1] 3
#> 
#> [[263]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[263]]$strand
#> [1] -1
#> 
#> [[263]]$start
#> [1] 29523292
#> 
#> [[263]]$end
#> [1] 29527702
#> 
#> [[263]]$id
#> [1] "ENSG00000213886"
#> 
#> [[263]]$canonical_transcript
#> [1] "ENST00000377050.4"
#> 
#> [[263]]$external_name
#> [1] "UBD"
#> 
#> [[263]]$description
#> [1] "ubiquitin D [Source:HGNC Symbol;Acc:18795]"
#> 
#> [[263]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[263]]$feature_type
#> [1] "gene"
#> 
#> [[263]]$biotype
#> [1] "protein_coding"
#> 
#> [[263]]$seq_region_name
#> [1] "6"
#> 
#> [[263]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[264]]
#> [[264]]$version
#> [1] 6
#> 
#> [[264]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[264]]$strand
#> [1] -1
#> 
#> [[264]]$gene_id
#> [1] "ENSG00000204681"
#> 
#> [[264]]$start
#> [1] 29523406
#> 
#> [[264]]$canonical_transcript
#> [1] "ENST00000377034.4"
#> 
#> [[264]]$external_name
#> [1] "GABBR1"
#> 
#> [[264]]$end
#> [1] 29601753
#> 
#> [[264]]$id
#> [1] "ENSG00000204681"
#> 
#> [[264]]$description
#> [1] "gamma-aminobutyric acid (GABA) B receptor, 1 [Source:HGNC Symbol;Acc:4070]"
#> 
#> [[264]]$biotype
#> [1] "protein_coding"
#> 
#> [[264]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[264]]$feature_type
#> [1] "gene"
#> 
#> [[264]]$source
#> [1] "ensembl_havana"
#> 
#> [[264]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[265]]
#> [[265]]$feature_type
#> [1] "gene"
#> 
#> [[265]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[265]]$biotype
#> [1] "protein_coding"
#> 
#> [[265]]$description
#> [1] "complement component 2 [Source:HGNC Symbol;Acc:1248]"
#> 
#> [[265]]$seq_region_name
#> [1] "6"
#> 
#> [[265]]$source
#> [1] "ensembl_havana"
#> 
#> [[265]]$gene_id
#> [1] "ENSG00000166278"
#> 
#> [[265]]$version
#> [1] 10
#> 
#> [[265]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[265]]$strand
#> [1] 1
#> 
#> [[265]]$end
#> [1] 31913449
#> 
#> [[265]]$id
#> [1] "ENSG00000166278"
#> 
#> [[265]]$canonical_transcript
#> [1] "ENST00000299367.5"
#> 
#> [[265]]$external_name
#> [1] "C2"
#> 
#> [[265]]$start
#> [1] 31865562
#> 
#> 
#> [[266]]
#> [[266]]$description
#> [1] "olfactory receptor, family 2, subfamily H, member 2 [Source:HGNC Symbol;Acc:8253]"
#> 
#> [[266]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[266]]$feature_type
#> [1] "gene"
#> 
#> [[266]]$biotype
#> [1] "protein_coding"
#> 
#> [[266]]$seq_region_name
#> [1] "6"
#> 
#> [[266]]$source
#> [1] "ensembl_havana"
#> 
#> [[266]]$gene_id
#> [1] "ENSG00000204657"
#> 
#> [[266]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[266]]$version
#> [1] 2
#> 
#> [[266]]$strand
#> [1] 1
#> 
#> [[266]]$start
#> [1] 29555683
#> 
#> [[266]]$id
#> [1] "ENSG00000204657"
#> 
#> [[266]]$end
#> [1] 29556745
#> 
#> [[266]]$canonical_transcript
#> [1] "ENST00000383640.2"
#> 
#> [[266]]$external_name
#> [1] "OR2H2"
#> 
#> 
#> [[267]]
#> [[267]]$seq_region_name
#> [1] "6"
#> 
#> [[267]]$source
#> [1] "ensembl_havana"
#> 
#> [[267]]$description
#> [1] "myelin oligodendrocyte glycoprotein [Source:HGNC Symbol;Acc:7197]"
#> 
#> [[267]]$feature_type
#> [1] "gene"
#> 
#> [[267]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[267]]$biotype
#> [1] "protein_coding"
#> 
#> [[267]]$start
#> [1] 29624758
#> 
#> [[267]]$id
#> [1] "ENSG00000204655"
#> 
#> [[267]]$end
#> [1] 29640149
#> 
#> [[267]]$external_name
#> [1] "MOG"
#> 
#> [[267]]$canonical_transcript
#> [1] "ENST00000376898.3"
#> 
#> [[267]]$gene_id
#> [1] "ENSG00000204655"
#> 
#> [[267]]$version
#> [1] 7
#> 
#> [[267]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[267]]$strand
#> [1] 1
#> 
#> 
#> [[268]]
#> [[268]]$biotype
#> [1] "protein_coding"
#> 
#> [[268]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[268]]$feature_type
#> [1] "gene"
#> 
#> [[268]]$description
#> [1] "zinc finger and BTB domain containing 12 [Source:HGNC Symbol;Acc:19066]"
#> 
#> [[268]]$source
#> [1] "ensembl_havana"
#> 
#> [[268]]$seq_region_name
#> [1] "6"
#> 
#> [[268]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[268]]$version
#> [1] 3
#> 
#> [[268]]$strand
#> [1] -1
#> 
#> [[268]]$gene_id
#> [1] "ENSG00000204366"
#> 
#> [[268]]$canonical_transcript
#> [1] "ENST00000375527.2"
#> 
#> [[268]]$external_name
#> [1] "ZBTB12"
#> 
#> [[268]]$id
#> [1] "ENSG00000204366"
#> 
#> [[268]]$end
#> [1] 31869769
#> 
#> [[268]]$start
#> [1] 31867384
#> 
#> 
#> [[269]]
#> [[269]]$version
#> [1] 4
#> 
#> [[269]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[269]]$strand
#> [1] 1
#> 
#> [[269]]$gene_id
#> [1] "ENSG00000243649"
#> 
#> [[269]]$start
#> [1] 31895475
#> 
#> [[269]]$external_name
#> [1] "CFB"
#> 
#> [[269]]$canonical_transcript
#> [1] "ENST00000425368.2"
#> 
#> [[269]]$id
#> [1] "ENSG00000243649"
#> 
#> [[269]]$end
#> [1] 31919861
#> 
#> [[269]]$description
#> [1] "complement factor B [Source:HGNC Symbol;Acc:1037]"
#> 
#> [[269]]$biotype
#> [1] "protein_coding"
#> 
#> [[269]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[269]]$feature_type
#> [1] "gene"
#> 
#> [[269]]$source
#> [1] "ensembl_havana"
#> 
#> [[269]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[270]]
#> [[270]]$seq_region_name
#> [1] "6"
#> 
#> [[270]]$source
#> [1] "havana"
#> 
#> [[270]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[270]]$feature_type
#> [1] "gene"
#> 
#> [[270]]$biotype
#> [1] "protein_coding"
#> 
#> [[270]]$description
#> [1] "Complement factor B; Uncharacterized protein; cDNA FLJ55673, highly similar to Complement factor B   [Source:UniProtKB/TrEMBL;Acc:B4E1Z4]"
#> 
#> [[270]]$id
#> [1] "ENSG00000244255"
#> 
#> [[270]]$end
#> [1] 31919825
#> 
#> [[270]]$canonical_transcript
#> [1] "ENST00000456570.1"
#> 
#> [[270]]$external_name
#> [1] "CFB"
#> 
#> [[270]]$start
#> [1] 31895475
#> 
#> [[270]]$gene_id
#> [1] "ENSG00000244255"
#> 
#> [[270]]$strand
#> [1] 1
#> 
#> [[270]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[270]]$version
#> [1] 1
#> 
#> 
#> [[271]]
#> [[271]]$external_name
#> [1] "NELFE"
#> 
#> [[271]]$canonical_transcript
#> [1] "ENST00000375429.3"
#> 
#> [[271]]$end
#> [1] 31926887
#> 
#> [[271]]$id
#> [1] "ENSG00000204356"
#> 
#> [[271]]$start
#> [1] 31919864
#> 
#> [[271]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[271]]$version
#> [1] 7
#> 
#> [[271]]$strand
#> [1] -1
#> 
#> [[271]]$gene_id
#> [1] "ENSG00000204356"
#> 
#> [[271]]$source
#> [1] "ensembl_havana"
#> 
#> [[271]]$seq_region_name
#> [1] "6"
#> 
#> [[271]]$biotype
#> [1] "protein_coding"
#> 
#> [[271]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[271]]$feature_type
#> [1] "gene"
#> 
#> [[271]]$description
#> [1] "negative elongation factor complex member E [Source:HGNC Symbol;Acc:13974]"
#> 
#> 
#> [[272]]
#> [[272]]$strand
#> [1] -1
#> 
#> [[272]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[272]]$version
#> [1] 5
#> 
#> [[272]]$gene_id
#> [1] "ENSG00000204644"
#> 
#> [[272]]$canonical_transcript
#> [1] "ENST00000488757.1"
#> 
#> [[272]]$external_name
#> [1] "ZFP57"
#> 
#> [[272]]$id
#> [1] "ENSG00000204644"
#> 
#> [[272]]$end
#> [1] 29648887
#> 
#> [[272]]$start
#> [1] 29640169
#> 
#> [[272]]$biotype
#> [1] "protein_coding"
#> 
#> [[272]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[272]]$feature_type
#> [1] "gene"
#> 
#> [[272]]$description
#> [1] "ZFP57 zinc finger protein [Source:HGNC Symbol;Acc:18791]"
#> 
#> [[272]]$source
#> [1] "ensembl_havana"
#> 
#> [[272]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[273]]
#> [[273]]$source
#> [1] "ensembl_havana"
#> 
#> [[273]]$seq_region_name
#> [1] "6"
#> 
#> [[273]]$description
#> [1] "major histocompatibility complex, class I, F [Source:HGNC Symbol;Acc:4963]"
#> 
#> [[273]]$biotype
#> [1] "protein_coding"
#> 
#> [[273]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[273]]$feature_type
#> [1] "gene"
#> 
#> [[273]]$start
#> [1] 29690552
#> 
#> [[273]]$canonical_transcript
#> [1] "ENST00000259951.7"
#> 
#> [[273]]$external_name
#> [1] "HLA-F"
#> 
#> [[273]]$id
#> [1] "ENSG00000204642"
#> 
#> [[273]]$end
#> [1] 29706305
#> 
#> [[273]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[273]]$version
#> [1] 9
#> 
#> [[273]]$strand
#> [1] 1
#> 
#> [[273]]$gene_id
#> [1] "ENSG00000204642"
#> 
#> 
#> [[274]]
#> [[274]]$biotype
#> [1] "protein_coding"
#> 
#> [[274]]$feature_type
#> [1] "gene"
#> 
#> [[274]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[274]]$description
#> [1] "superkiller viralicidic activity 2-like (S. cerevisiae) [Source:HGNC Symbol;Acc:10898]"
#> 
#> [[274]]$source
#> [1] "ensembl_havana"
#> 
#> [[274]]$seq_region_name
#> [1] "6"
#> 
#> [[274]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[274]]$version
#> [1] 7
#> 
#> [[274]]$strand
#> [1] 1
#> 
#> [[274]]$gene_id
#> [1] "ENSG00000204351"
#> 
#> [[274]]$canonical_transcript
#> [1] "ENST00000375394.2"
#> 
#> [[274]]$external_name
#> [1] "SKIV2L"
#> 
#> [[274]]$id
#> [1] "ENSG00000204351"
#> 
#> [[274]]$end
#> [1] 31937532
#> 
#> [[274]]$start
#> [1] 31926857
#> 
#> 
#> [[275]]
#> [[275]]$start
#> [1] 29794744
#> 
#> [[275]]$id
#> [1] "ENSG00000204632"
#> 
#> [[275]]$end
#> [1] 29798902
#> 
#> [[275]]$external_name
#> [1] "HLA-G"
#> 
#> [[275]]$canonical_transcript
#> [1] "ENST00000428701.1"
#> 
#> [[275]]$gene_id
#> [1] "ENSG00000204632"
#> 
#> [[275]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[275]]$version
#> [1] 7
#> 
#> [[275]]$strand
#> [1] 1
#> 
#> [[275]]$seq_region_name
#> [1] "6"
#> 
#> [[275]]$source
#> [1] "ensembl_havana"
#> 
#> [[275]]$description
#> [1] "major histocompatibility complex, class I, G [Source:HGNC Symbol;Acc:4964]"
#> 
#> [[275]]$feature_type
#> [1] "gene"
#> 
#> [[275]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[275]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[276]]
#> [[276]]$gene_id
#> [1] "ENSG00000206503"
#> 
#> [[276]]$strand
#> [1] 1
#> 
#> [[276]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[276]]$version
#> [1] 7
#> 
#> [[276]]$id
#> [1] "ENSG00000206503"
#> 
#> [[276]]$end
#> [1] 29913661
#> 
#> [[276]]$canonical_transcript
#> [1] "ENST00000396634.1"
#> 
#> [[276]]$external_name
#> [1] "HLA-A"
#> 
#> [[276]]$start
#> [1] 29909037
#> 
#> [[276]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[276]]$feature_type
#> [1] "gene"
#> 
#> [[276]]$biotype
#> [1] "protein_coding"
#> 
#> [[276]]$description
#> [1] "major histocompatibility complex, class I, A [Source:HGNC Symbol;Acc:4931]"
#> 
#> [[276]]$seq_region_name
#> [1] "6"
#> 
#> [[276]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[277]]
#> [[277]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[277]]$version
#> [1] 5
#> 
#> [[277]]$strand
#> [1] -1
#> 
#> [[277]]$gene_id
#> [1] "ENSG00000204348"
#> 
#> [[277]]$start
#> [1] 31937587
#> 
#> [[277]]$external_name
#> [1] "DXO"
#> 
#> [[277]]$canonical_transcript
#> [1] "ENST00000375349.3"
#> 
#> [[277]]$id
#> [1] "ENSG00000204348"
#> 
#> [[277]]$end
#> [1] 31940069
#> 
#> [[277]]$description
#> [1] "decapping exoribonuclease [Source:HGNC Symbol;Acc:2992]"
#> 
#> [[277]]$biotype
#> [1] "protein_coding"
#> 
#> [[277]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[277]]$feature_type
#> [1] "gene"
#> 
#> [[277]]$source
#> [1] "ensembl_havana"
#> 
#> [[277]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[278]]
#> [[278]]$id
#> [1] "ENSG00000204344"
#> 
#> [[278]]$end
#> [1] 31950598
#> 
#> [[278]]$canonical_transcript
#> [1] "ENST00000375333.2"
#> 
#> [[278]]$external_name
#> [1] "STK19"
#> 
#> [[278]]$start
#> [1] 31938868
#> 
#> [[278]]$gene_id
#> [1] "ENSG00000204344"
#> 
#> [[278]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[278]]$strand
#> [1] 1
#> 
#> [[278]]$version
#> [1] 10
#> 
#> [[278]]$seq_region_name
#> [1] "6"
#> 
#> [[278]]$source
#> [1] "ensembl_havana"
#> 
#> [[278]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[278]]$feature_type
#> [1] "gene"
#> 
#> [[278]]$biotype
#> [1] "protein_coding"
#> 
#> [[278]]$description
#> [1] "serine/threonine kinase 19 [Source:HGNC Symbol;Acc:11398]"
#> 
#> 
#> [[279]]
#> [[279]]$external_name
#> [1] "ZNRD1"
#> 
#> [[279]]$canonical_transcript
#> [1] "ENST00000332435.5"
#> 
#> [[279]]$id
#> [1] "ENSG00000066379"
#> 
#> [[279]]$end
#> [1] 30032686
#> 
#> [[279]]$start
#> [1] 30026676
#> 
#> [[279]]$strand
#> [1] 1
#> 
#> [[279]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[279]]$version
#> [1] 10
#> 
#> [[279]]$gene_id
#> [1] "ENSG00000066379"
#> 
#> [[279]]$source
#> [1] "ensembl_havana"
#> 
#> [[279]]$seq_region_name
#> [1] "6"
#> 
#> [[279]]$biotype
#> [1] "protein_coding"
#> 
#> [[279]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[279]]$feature_type
#> [1] "gene"
#> 
#> [[279]]$description
#> [1] "zinc ribbon domain containing 1 [Source:HGNC Symbol;Acc:13182]"
#> 
#> 
#> [[280]]
#> [[280]]$start
#> [1] 30034486
#> 
#> [[280]]$canonical_transcript
#> [1] "ENST00000376772.3"
#> 
#> [[280]]$external_name
#> [1] "PPP1R11"
#> 
#> [[280]]$end
#> [1] 30038110
#> 
#> [[280]]$id
#> [1] "ENSG00000204619"
#> 
#> [[280]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[280]]$strand
#> [1] 1
#> 
#> [[280]]$version
#> [1] 3
#> 
#> [[280]]$gene_id
#> [1] "ENSG00000204619"
#> 
#> [[280]]$source
#> [1] "ensembl_havana"
#> 
#> [[280]]$seq_region_name
#> [1] "6"
#> 
#> [[280]]$description
#> [1] "protein phosphatase 1, regulatory (inhibitor) subunit 11 [Source:HGNC Symbol;Acc:9285]"
#> 
#> [[280]]$biotype
#> [1] "protein_coding"
#> 
#> [[280]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[280]]$feature_type
#> [1] "gene"
#> 
#> 
#> [[281]]
#> [[281]]$start
#> [1] 31949801
#> 
#> [[281]]$id
#> [1] "ENSG00000244731"
#> 
#> [[281]]$end
#> [1] 31970458
#> 
#> [[281]]$external_name
#> [1] "C4A"
#> 
#> [[281]]$canonical_transcript
#> [1] "ENST00000428956.2"
#> 
#> [[281]]$gene_id
#> [1] "ENSG00000244731"
#> 
#> [[281]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[281]]$strand
#> [1] 1
#> 
#> [[281]]$version
#> [1] 3
#> 
#> [[281]]$seq_region_name
#> [1] "6"
#> 
#> [[281]]$source
#> [1] "ensembl_havana"
#> 
#> [[281]]$description
#> [1] "complement component 4A (Rodgers blood group) [Source:HGNC Symbol;Acc:1323]"
#> 
#> [[281]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[281]]$feature_type
#> [1] "gene"
#> 
#> [[281]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[282]]
#> [[282]]$id
#> [1] "ENSG00000204618"
#> 
#> [[282]]$end
#> [1] 30043664
#> 
#> [[282]]$canonical_transcript
#> [1] "ENST00000244360.6"
#> 
#> [[282]]$external_name
#> [1] "RNF39"
#> 
#> [[282]]$start
#> [1] 30038043
#> 
#> [[282]]$gene_id
#> [1] "ENSG00000204618"
#> 
#> [[282]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[282]]$version
#> [1] 4
#> 
#> [[282]]$strand
#> [1] -1
#> 
#> [[282]]$seq_region_name
#> [1] "6"
#> 
#> [[282]]$source
#> [1] "ensembl_havana"
#> 
#> [[282]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[282]]$feature_type
#> [1] "gene"
#> 
#> [[282]]$biotype
#> [1] "protein_coding"
#> 
#> [[282]]$description
#> [1] "ring finger protein 39 [Source:HGNC Symbol;Acc:18064]"
#> 
#> 
#> [[283]]
#> [[283]]$biotype
#> [1] "protein_coding"
#> 
#> [[283]]$feature_type
#> [1] "gene"
#> 
#> [[283]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[283]]$description
#> [1] "tripartite motif containing 31 [Source:HGNC Symbol;Acc:16289]"
#> 
#> [[283]]$source
#> [1] "ensembl_havana"
#> 
#> [[283]]$seq_region_name
#> [1] "6"
#> 
#> [[283]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[283]]$strand
#> [1] -1
#> 
#> [[283]]$version
#> [1] 6
#> 
#> [[283]]$gene_id
#> [1] "ENSG00000204616"
#> 
#> [[283]]$canonical_transcript
#> [1] "ENST00000376734.3"
#> 
#> [[283]]$external_name
#> [1] "TRIM31"
#> 
#> [[283]]$id
#> [1] "ENSG00000204616"
#> 
#> [[283]]$end
#> [1] 30080883
#> 
#> [[283]]$start
#> [1] 30070674
#> 
#> 
#> [[284]]
#> [[284]]$source
#> [1] "ensembl_havana"
#> 
#> [[284]]$seq_region_name
#> [1] "6"
#> 
#> [[284]]$description
#> [1] "tripartite motif containing 40 [Source:HGNC Symbol;Acc:18736]"
#> 
#> [[284]]$biotype
#> [1] "protein_coding"
#> 
#> [[284]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[284]]$feature_type
#> [1] "gene"
#> 
#> [[284]]$start
#> [1] 30103885
#> 
#> [[284]]$canonical_transcript
#> [1] "ENST00000307859.4"
#> 
#> [[284]]$external_name
#> [1] "TRIM40"
#> 
#> [[284]]$id
#> [1] "ENSG00000204614"
#> 
#> [[284]]$end
#> [1] 30116512
#> 
#> [[284]]$strand
#> [1] 1
#> 
#> [[284]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[284]]$version
#> [1] 4
#> 
#> [[284]]$gene_id
#> [1] "ENSG00000204614"
#> 
#> 
#> [[285]]
#> [[285]]$gene_id
#> [1] "ENSG00000204613"
#> 
#> [[285]]$strand
#> [1] -1
#> 
#> [[285]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[285]]$version
#> [1] 6
#> 
#> [[285]]$start
#> [1] 30119722
#> 
#> [[285]]$end
#> [1] 30128711
#> 
#> [[285]]$id
#> [1] "ENSG00000204613"
#> 
#> [[285]]$canonical_transcript
#> [1] "ENST00000449742.2"
#> 
#> [[285]]$external_name
#> [1] "TRIM10"
#> 
#> [[285]]$description
#> [1] "tripartite motif containing 10 [Source:HGNC Symbol;Acc:10072]"
#> 
#> [[285]]$feature_type
#> [1] "gene"
#> 
#> [[285]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[285]]$biotype
#> [1] "protein_coding"
#> 
#> [[285]]$seq_region_name
#> [1] "6"
#> 
#> [[285]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[286]]
#> [[286]]$start
#> [1] 30130993
#> 
#> [[286]]$id
#> [1] "ENSG00000204610"
#> 
#> [[286]]$end
#> [1] 30140473
#> 
#> [[286]]$external_name
#> [1] "TRIM15"
#> 
#> [[286]]$canonical_transcript
#> [1] "ENST00000376694.4"
#> 
#> [[286]]$gene_id
#> [1] "ENSG00000204610"
#> 
#> [[286]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[286]]$version
#> [1] 8
#> 
#> [[286]]$strand
#> [1] 1
#> 
#> [[286]]$seq_region_name
#> [1] "6"
#> 
#> [[286]]$source
#> [1] "ensembl_havana"
#> 
#> [[286]]$description
#> [1] "tripartite motif containing 15 [Source:HGNC Symbol;Acc:16284]"
#> 
#> [[286]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[286]]$feature_type
#> [1] "gene"
#> 
#> [[286]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[287]]
#> [[287]]$description
#> [1] "tripartite motif containing 26 [Source:HGNC Symbol;Acc:12962]"
#> 
#> [[287]]$biotype
#> [1] "protein_coding"
#> 
#> [[287]]$feature_type
#> [1] "gene"
#> 
#> [[287]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[287]]$source
#> [1] "ensembl_havana"
#> 
#> [[287]]$seq_region_name
#> [1] "6"
#> 
#> [[287]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[287]]$strand
#> [1] -1
#> 
#> [[287]]$version
#> [1] 4
#> 
#> [[287]]$gene_id
#> [1] "ENSG00000234127"
#> 
#> [[287]]$start
#> [1] 30152232
#> 
#> [[287]]$external_name
#> [1] "TRIM26"
#> 
#> [[287]]$canonical_transcript
#> [1] "ENST00000454678.2"
#> 
#> [[287]]$id
#> [1] "ENSG00000234127"
#> 
#> [[287]]$end
#> [1] 30181204
#> 
#> 
#> [[288]]
#> [[288]]$end
#> [1] 30311506
#> 
#> [[288]]$id
#> [1] "ENSG00000204599"
#> 
#> [[288]]$external_name
#> [1] "TRIM39"
#> 
#> [[288]]$canonical_transcript
#> [1] "ENST00000376656.4"
#> 
#> [[288]]$start
#> [1] 30294256
#> 
#> [[288]]$gene_id
#> [1] "ENSG00000204599"
#> 
#> [[288]]$version
#> [1] 10
#> 
#> [[288]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[288]]$strand
#> [1] 1
#> 
#> [[288]]$seq_region_name
#> [1] "6"
#> 
#> [[288]]$source
#> [1] "ensembl_havana"
#> 
#> [[288]]$feature_type
#> [1] "gene"
#> 
#> [[288]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[288]]$biotype
#> [1] "protein_coding"
#> 
#> [[288]]$description
#> [1] "tripartite motif containing 39 [Source:HGNC Symbol;Acc:10065]"
#> 
#> 
#> [[289]]
#> [[289]]$feature_type
#> [1] "gene"
#> 
#> [[289]]$logic_name
#> [1] "ensembl_homo_sapiens_37"
#> 
#> [[289]]$biotype
#> [1] "protein_coding"
#> 
#> [[289]]$description
#> [1] "Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:M0R2J5]"
#> 
#> [[289]]$seq_region_name
#> [1] "6"
#> 
#> [[289]]$source
#> [1] "ensembl"
#> 
#> [[289]]$gene_id
#> [1] "ENSG00000268923"
#> 
#> [[289]]$version
#> [1] 1
#> 
#> [[289]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[289]]$strand
#> [1] -1
#> 
#> [[289]]$id
#> [1] "ENSG00000268923"
#> 
#> [[289]]$end
#> [1] 31974881
#> 
#> [[289]]$external_name
#> [1] "AL645922.1"
#> 
#> [[289]]$canonical_transcript
#> [1] "ENST00000594256.1"
#> 
#> [[289]]$start
#> [1] 31973945
#> 
#> 
#> [[290]]
#> [[290]]$description
#> [1] "complement component 4B (Chido blood group) [Source:HGNC Symbol;Acc:1324]"
#> 
#> [[290]]$biotype
#> [1] "protein_coding"
#> 
#> [[290]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[290]]$feature_type
#> [1] "gene"
#> 
#> [[290]]$source
#> [1] "ensembl_havana"
#> 
#> [[290]]$seq_region_name
#> [1] "6"
#> 
#> [[290]]$strand
#> [1] 1
#> 
#> [[290]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[290]]$version
#> [1] 4
#> 
#> [[290]]$gene_id
#> [1] "ENSG00000224389"
#> 
#> [[290]]$start
#> [1] 31982539
#> 
#> [[290]]$external_name
#> [1] "C4B"
#> 
#> [[290]]$canonical_transcript
#> [1] "ENST00000435363.2"
#> 
#> [[290]]$id
#> [1] "ENSG00000224389"
#> 
#> [[290]]$end
#> [1] 32003195
#> 
#> 
#> [[291]]
#> [[291]]$gene_id
#> [1] "ENSG00000248167"
#> 
#> [[291]]$version
#> [1] 3
#> 
#> [[291]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[291]]$strand
#> [1] 1
#> 
#> [[291]]$id
#> [1] "ENSG00000248167"
#> 
#> [[291]]$end
#> [1] 30314631
#> 
#> [[291]]$canonical_transcript
#> [1] "ENST00000513556.1"
#> 
#> [[291]]$external_name
#> [1] "TRIM39-RPP21"
#> 
#> [[291]]$start
#> [1] 30297359
#> 
#> [[291]]$feature_type
#> [1] "gene"
#> 
#> [[291]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[291]]$biotype
#> [1] "protein_coding"
#> 
#> [[291]]$description
#> [1] "TRIM39-RPP21 readthrough [Source:HGNC Symbol;Acc:38845]"
#> 
#> [[291]]$seq_region_name
#> [1] "6"
#> 
#> [[291]]$source
#> [1] "havana"
#> 
#> 
#> [[292]]
#> [[292]]$version
#> [1] 1
#> 
#> [[292]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[292]]$strand
#> [1] 1
#> 
#> [[292]]$gene_id
#> [1] "ENSG00000241370"
#> 
#> [[292]]$canonical_transcript
#> [1] "ENST00000433076.2"
#> 
#> [[292]]$external_name
#> [1] "RPP21"
#> 
#> [[292]]$end
#> [1] 30314661
#> 
#> [[292]]$id
#> [1] "ENSG00000241370"
#> 
#> [[292]]$start
#> [1] 30312908
#> 
#> [[292]]$biotype
#> [1] "protein_coding"
#> 
#> [[292]]$feature_type
#> [1] "gene"
#> 
#> [[292]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[292]]$description
#> [1] "ribonuclease P/MRP 21kDa subunit [Source:HGNC Symbol;Acc:21300]"
#> 
#> [[292]]$source
#> [1] "ensembl_havana"
#> 
#> [[292]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[293]]
#> [[293]]$seq_region_name
#> [1] "6"
#> 
#> [[293]]$source
#> [1] "ensembl_havana"
#> 
#> [[293]]$description
#> [1] "major histocompatibility complex, class I, E [Source:HGNC Symbol;Acc:4962]"
#> 
#> [[293]]$feature_type
#> [1] "gene"
#> 
#> [[293]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[293]]$biotype
#> [1] "protein_coding"
#> 
#> [[293]]$start
#> [1] 30457244
#> 
#> [[293]]$id
#> [1] "ENSG00000204592"
#> 
#> [[293]]$end
#> [1] 30461982
#> 
#> [[293]]$external_name
#> [1] "HLA-E"
#> 
#> [[293]]$canonical_transcript
#> [1] "ENST00000376630.4"
#> 
#> [[293]]$gene_id
#> [1] "ENSG00000204592"
#> 
#> [[293]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[293]]$version
#> [1] 5
#> 
#> [[293]]$strand
#> [1] 1
#> 
#> 
#> [[294]]
#> [[294]]$gene_id
#> [1] "ENSG00000204590"
#> 
#> [[294]]$version
#> [1] 8
#> 
#> [[294]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[294]]$strand
#> [1] -1
#> 
#> [[294]]$id
#> [1] "ENSG00000204590"
#> 
#> [[294]]$end
#> [1] 30524951
#> 
#> [[294]]$external_name
#> [1] "GNL1"
#> 
#> [[294]]$canonical_transcript
#> [1] "ENST00000376621.3"
#> 
#> [[294]]$start
#> [1] 30509154
#> 
#> [[294]]$feature_type
#> [1] "gene"
#> 
#> [[294]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[294]]$biotype
#> [1] "protein_coding"
#> 
#> [[294]]$description
#> [1] "guanine nucleotide binding protein-like 1 [Source:HGNC Symbol;Acc:4413]"
#> 
#> [[294]]$seq_region_name
#> [1] "6"
#> 
#> [[294]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[295]]
#> [[295]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[295]]$version
#> [1] 2
#> 
#> [[295]]$strand
#> [1] 1
#> 
#> [[295]]$gene_id
#> [1] "ENSG00000231852"
#> 
#> [[295]]$canonical_transcript
#> [1] "ENST00000418967.2"
#> 
#> [[295]]$external_name
#> [1] "CYP21A2"
#> 
#> [[295]]$id
#> [1] "ENSG00000231852"
#> 
#> [[295]]$end
#> [1] 32009447
#> 
#> [[295]]$start
#> [1] 32006042
#> 
#> [[295]]$biotype
#> [1] "protein_coding"
#> 
#> [[295]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[295]]$feature_type
#> [1] "gene"
#> 
#> [[295]]$description
#> [1] "cytochrome P450, family 21, subfamily A, polypeptide 2 [Source:HGNC Symbol;Acc:2600]"
#> 
#> [[295]]$source
#> [1] "ensembl_havana"
#> 
#> [[295]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[296]]
#> [[296]]$seq_region_name
#> [1] "6"
#> 
#> [[296]]$source
#> [1] "ensembl_havana"
#> 
#> [[296]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[296]]$feature_type
#> [1] "gene"
#> 
#> [[296]]$biotype
#> [1] "protein_coding"
#> 
#> [[296]]$description
#> [1] "proline rich 3 [Source:HGNC Symbol;Acc:21149]"
#> 
#> [[296]]$end
#> [1] 30531500
#> 
#> [[296]]$id
#> [1] "ENSG00000204576"
#> 
#> [[296]]$external_name
#> [1] "PRR3"
#> 
#> [[296]]$canonical_transcript
#> [1] "ENST00000376560.3"
#> 
#> [[296]]$start
#> [1] 30524663
#> 
#> [[296]]$gene_id
#> [1] "ENSG00000204576"
#> 
#> [[296]]$strand
#> [1] 1
#> 
#> [[296]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[296]]$version
#> [1] 7
#> 
#> 
#> [[297]]
#> [[297]]$start
#> [1] 32008931
#> 
#> [[297]]$end
#> [1] 32083111
#> 
#> [[297]]$id
#> [1] "ENSG00000168477"
#> 
#> [[297]]$external_name
#> [1] "TNXB"
#> 
#> [[297]]$canonical_transcript
#> [1] "ENST00000451343.1"
#> 
#> [[297]]$gene_id
#> [1] "ENSG00000168477"
#> 
#> [[297]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[297]]$strand
#> [1] -1
#> 
#> [[297]]$version
#> [1] 13
#> 
#> [[297]]$seq_region_name
#> [1] "6"
#> 
#> [[297]]$source
#> [1] "ensembl_havana"
#> 
#> [[297]]$description
#> [1] "tenascin XB [Source:HGNC Symbol;Acc:11976]"
#> 
#> [[297]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[297]]$feature_type
#> [1] "gene"
#> 
#> [[297]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[298]]
#> [[298]]$biotype
#> [1] "protein_coding"
#> 
#> [[298]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[298]]$feature_type
#> [1] "gene"
#> 
#> [[298]]$description
#> [1] "ATP-binding cassette, sub-family F (GCN20), member 1 [Source:HGNC Symbol;Acc:70]"
#> 
#> [[298]]$source
#> [1] "ensembl_havana"
#> 
#> [[298]]$seq_region_name
#> [1] "6"
#> 
#> [[298]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[298]]$strand
#> [1] 1
#> 
#> [[298]]$version
#> [1] 8
#> 
#> [[298]]$gene_id
#> [1] "ENSG00000204574"
#> 
#> [[298]]$canonical_transcript
#> [1] "ENST00000326195.8"
#> 
#> [[298]]$external_name
#> [1] "ABCF1"
#> 
#> [[298]]$end
#> [1] 30564956
#> 
#> [[298]]$id
#> [1] "ENSG00000204574"
#> 
#> [[298]]$start
#> [1] 30539153
#> 
#> 
#> [[299]]
#> [[299]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[299]]$version
#> [1] 5
#> 
#> [[299]]$strand
#> [1] -1
#> 
#> [[299]]$gene_id
#> [1] "ENSG00000204569"
#> 
#> [[299]]$start
#> [1] 30568177
#> 
#> [[299]]$external_name
#> [1] "PPP1R10"
#> 
#> [[299]]$canonical_transcript
#> [1] "ENST00000376511.2"
#> 
#> [[299]]$id
#> [1] "ENSG00000204569"
#> 
#> [[299]]$end
#> [1] 30586389
#> 
#> [[299]]$description
#> [1] "protein phosphatase 1, regulatory subunit 10 [Source:HGNC Symbol;Acc:9284]"
#> 
#> [[299]]$biotype
#> [1] "protein_coding"
#> 
#> [[299]]$feature_type
#> [1] "gene"
#> 
#> [[299]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[299]]$source
#> [1] "ensembl_havana"
#> 
#> [[299]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[300]]
#> [[300]]$start
#> [1] 32065953
#> 
#> [[300]]$end
#> [1] 32096030
#> 
#> [[300]]$id
#> [1] "ENSG00000213676"
#> 
#> [[300]]$external_name
#> [1] "ATF6B"
#> 
#> [[300]]$canonical_transcript
#> [1] "ENST00000375203.3"
#> 
#> [[300]]$gene_id
#> [1] "ENSG00000213676"
#> 
#> [[300]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[300]]$strand
#> [1] -1
#> 
#> [[300]]$version
#> [1] 6
#> 
#> [[300]]$seq_region_name
#> [1] "6"
#> 
#> [[300]]$source
#> [1] "ensembl_havana"
#> 
#> [[300]]$description
#> [1] "activating transcription factor 6 beta [Source:HGNC Symbol;Acc:2349]"
#> 
#> [[300]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[300]]$feature_type
#> [1] "gene"
#> 
#> [[300]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[301]]
#> [[301]]$id
#> [1] "ENSG00000204568"
#> 
#> [[301]]$end
#> [1] 30594172
#> 
#> [[301]]$canonical_transcript
#> [1] "ENST00000259873.4"
#> 
#> [[301]]$external_name
#> [1] "MRPS18B"
#> 
#> [[301]]$start
#> [1] 30585486
#> 
#> [[301]]$gene_id
#> [1] "ENSG00000204568"
#> 
#> [[301]]$strand
#> [1] 1
#> 
#> [[301]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[301]]$version
#> [1] 7
#> 
#> [[301]]$seq_region_name
#> [1] "6"
#> 
#> [[301]]$source
#> [1] "ensembl_havana"
#> 
#> [[301]]$feature_type
#> [1] "gene"
#> 
#> [[301]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[301]]$biotype
#> [1] "protein_coding"
#> 
#> [[301]]$description
#> [1] "mitochondrial ribosomal protein S18B [Source:HGNC Symbol;Acc:14516]"
#> 
#> 
#> [[302]]
#> [[302]]$start
#> [1] 30594619
#> 
#> [[302]]$canonical_transcript
#> [1] "ENST00000330083.5"
#> 
#> [[302]]$external_name
#> [1] "ATAT1"
#> 
#> [[302]]$id
#> [1] "ENSG00000137343"
#> 
#> [[302]]$end
#> [1] 30614600
#> 
#> [[302]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[302]]$strand
#> [1] 1
#> 
#> [[302]]$version
#> [1] 13
#> 
#> [[302]]$gene_id
#> [1] "ENSG00000137343"
#> 
#> [[302]]$source
#> [1] "ensembl_havana"
#> 
#> [[302]]$seq_region_name
#> [1] "6"
#> 
#> [[302]]$description
#> [1] "alpha tubulin acetyltransferase 1 [Source:HGNC Symbol;Acc:21186]"
#> 
#> [[302]]$biotype
#> [1] "protein_coding"
#> 
#> [[302]]$feature_type
#> [1] "gene"
#> 
#> [[302]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> 
#> [[303]]
#> [[303]]$source
#> [1] "ensembl_havana"
#> 
#> [[303]]$seq_region_name
#> [1] "6"
#> 
#> [[303]]$biotype
#> [1] "protein_coding"
#> 
#> [[303]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[303]]$feature_type
#> [1] "gene"
#> 
#> [[303]]$description
#> [1] "FK506 binding protein like [Source:HGNC Symbol;Acc:13949]"
#> 
#> [[303]]$external_name
#> [1] "FKBPL"
#> 
#> [[303]]$canonical_transcript
#> [1] "ENST00000375156.3"
#> 
#> [[303]]$id
#> [1] "ENSG00000204315"
#> 
#> [[303]]$end
#> [1] 32098068
#> 
#> [[303]]$start
#> [1] 32096484
#> 
#> [[303]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[303]]$strand
#> [1] -1
#> 
#> [[303]]$version
#> [1] 3
#> 
#> [[303]]$gene_id
#> [1] "ENSG00000204315"
#> 
#> 
#> [[304]]
#> [[304]]$external_name
#> [1] "PRRT1"
#> 
#> [[304]]$canonical_transcript
#> [1] "ENST00000211413.5"
#> 
#> [[304]]$id
#> [1] "ENSG00000204314"
#> 
#> [[304]]$end
#> [1] 32122150
#> 
#> [[304]]$start
#> [1] 32116136
#> 
#> [[304]]$strand
#> [1] -1
#> 
#> [[304]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[304]]$version
#> [1] 6
#> 
#> [[304]]$gene_id
#> [1] "ENSG00000204314"
#> 
#> [[304]]$source
#> [1] "ensembl_havana"
#> 
#> [[304]]$seq_region_name
#> [1] "6"
#> 
#> [[304]]$biotype
#> [1] "protein_coding"
#> 
#> [[304]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[304]]$feature_type
#> [1] "gene"
#> 
#> [[304]]$description
#> [1] "proline-rich transmembrane protein 1 [Source:HGNC Symbol;Acc:13943]"
#> 
#> 
#> [[305]]
#> [[305]]$seq_region_name
#> [1] "6"
#> 
#> [[305]]$source
#> [1] "ensembl_havana"
#> 
#> [[305]]$feature_type
#> [1] "gene"
#> 
#> [[305]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[305]]$biotype
#> [1] "protein_coding"
#> 
#> [[305]]$description
#> [1] "palmitoyl-protein thioesterase 2 [Source:HGNC Symbol;Acc:9326]"
#> 
#> [[305]]$id
#> [1] "ENSG00000221988"
#> 
#> [[305]]$end
#> [1] 32134011
#> 
#> [[305]]$canonical_transcript
#> [1] "ENST00000361568.2"
#> 
#> [[305]]$external_name
#> [1] "PPT2"
#> 
#> [[305]]$start
#> [1] 32121218
#> 
#> [[305]]$gene_id
#> [1] "ENSG00000221988"
#> 
#> [[305]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[305]]$strand
#> [1] 1
#> 
#> [[305]]$version
#> [1] 8
#> 
#> 
#> [[306]]
#> [[306]]$seq_region_name
#> [1] "6"
#> 
#> [[306]]$source
#> [1] "havana"
#> 
#> [[306]]$description
#> [1] "PPT2-EGFL8 readthrough (NMD candidate) [Source:HGNC Symbol;Acc:48343]"
#> 
#> [[306]]$feature_type
#> [1] "gene"
#> 
#> [[306]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[306]]$biotype
#> [1] "protein_coding"
#> 
#> [[306]]$start
#> [1] 32121622
#> 
#> [[306]]$end
#> [1] 32139755
#> 
#> [[306]]$id
#> [1] "ENSG00000258388"
#> 
#> [[306]]$canonical_transcript
#> [1] "ENST00000422437.1"
#> 
#> [[306]]$external_name
#> [1] "PPT2-EGFL8"
#> 
#> [[306]]$gene_id
#> [1] "ENSG00000258388"
#> 
#> [[306]]$strand
#> [1] 1
#> 
#> [[306]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[306]]$version
#> [1] 3
#> 
#> 
#> [[307]]
#> [[307]]$gene_id
#> [1] "ENSG00000241404"
#> 
#> [[307]]$strand
#> [1] 1
#> 
#> [[307]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[307]]$version
#> [1] 2
#> 
#> [[307]]$start
#> [1] 32132360
#> 
#> [[307]]$id
#> [1] "ENSG00000241404"
#> 
#> [[307]]$end
#> [1] 32136058
#> 
#> [[307]]$canonical_transcript
#> [1] "ENST00000395512.1"
#> 
#> [[307]]$external_name
#> [1] "EGFL8"
#> 
#> [[307]]$description
#> [1] "EGF-like-domain, multiple 8 [Source:HGNC Symbol;Acc:13944]"
#> 
#> [[307]]$feature_type
#> [1] "gene"
#> 
#> [[307]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[307]]$biotype
#> [1] "protein_coding"
#> 
#> [[307]]$seq_region_name
#> [1] "6"
#> 
#> [[307]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[308]]
#> [[308]]$source
#> [1] "ensembl_havana"
#> 
#> [[308]]$seq_region_name
#> [1] "6"
#> 
#> [[308]]$description
#> [1] "chromosome 6 open reading frame 136 [Source:HGNC Symbol;Acc:21301]"
#> 
#> [[308]]$biotype
#> [1] "protein_coding"
#> 
#> [[308]]$feature_type
#> [1] "gene"
#> 
#> [[308]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[308]]$start
#> [1] 30614816
#> 
#> [[308]]$external_name
#> [1] "C6orf136"
#> 
#> [[308]]$canonical_transcript
#> [1] "ENST00000293604.6"
#> 
#> [[308]]$id
#> [1] "ENSG00000204564"
#> 
#> [[308]]$end
#> [1] 30620987
#> 
#> [[308]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[308]]$strand
#> [1] 1
#> 
#> [[308]]$version
#> [1] 7
#> 
#> [[308]]$gene_id
#> [1] "ENSG00000204564"
#> 
#> 
#> [[309]]
#> [[309]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[309]]$strand
#> [1] -1
#> 
#> [[309]]$version
#> [1] 6
#> 
#> [[309]]$gene_id
#> [1] "ENSG00000204310"
#> 
#> [[309]]$start
#> [1] 32135989
#> 
#> [[309]]$canonical_transcript
#> [1] "ENST00000395499.1"
#> 
#> [[309]]$external_name
#> [1] "AGPAT1"
#> 
#> [[309]]$id
#> [1] "ENSG00000204310"
#> 
#> [[309]]$end
#> [1] 32145873
#> 
#> [[309]]$description
#> [1] "1-acylglycerol-3-phosphate O-acyltransferase 1 [Source:HGNC Symbol;Acc:324]"
#> 
#> [[309]]$biotype
#> [1] "protein_coding"
#> 
#> [[309]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[309]]$feature_type
#> [1] "gene"
#> 
#> [[309]]$source
#> [1] "ensembl_havana"
#> 
#> [[309]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[310]]
#> [[310]]$end
#> [1] 32151930
#> 
#> [[310]]$id
#> [1] "ENSG00000204308"
#> 
#> [[310]]$external_name
#> [1] "RNF5"
#> 
#> [[310]]$canonical_transcript
#> [1] "ENST00000375094.3"
#> 
#> [[310]]$start
#> [1] 32146131
#> 
#> [[310]]$gene_id
#> [1] "ENSG00000204308"
#> 
#> [[310]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[310]]$strand
#> [1] 1
#> 
#> [[310]]$version
#> [1] 6
#> 
#> [[310]]$seq_region_name
#> [1] "6"
#> 
#> [[310]]$source
#> [1] "ensembl_havana"
#> 
#> [[310]]$feature_type
#> [1] "gene"
#> 
#> [[310]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[310]]$biotype
#> [1] "protein_coding"
#> 
#> [[310]]$description
#> [1] "ring finger protein 5, E3 ubiquitin protein ligase [Source:HGNC Symbol;Acc:10068]"
#> 
#> 
#> [[311]]
#> [[311]]$feature_type
#> [1] "gene"
#> 
#> [[311]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[311]]$biotype
#> [1] "protein_coding"
#> 
#> [[311]]$description
#> [1] "DEAH (Asp-Glu-Ala-His) box polypeptide 16 [Source:HGNC Symbol;Acc:2739]"
#> 
#> [[311]]$seq_region_name
#> [1] "6"
#> 
#> [[311]]$source
#> [1] "ensembl_havana"
#> 
#> [[311]]$gene_id
#> [1] "ENSG00000204560"
#> 
#> [[311]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[311]]$version
#> [1] 5
#> 
#> [[311]]$strand
#> [1] -1
#> 
#> [[311]]$end
#> [1] 30640814
#> 
#> [[311]]$id
#> [1] "ENSG00000204560"
#> 
#> [[311]]$canonical_transcript
#> [1] "ENST00000376442.3"
#> 
#> [[311]]$external_name
#> [1] "DHX16"
#> 
#> [[311]]$start
#> [1] 30620896
#> 
#> 
#> [[312]]
#> [[312]]$source
#> [1] "ensembl_havana"
#> 
#> [[312]]$seq_region_name
#> [1] "6"
#> 
#> [[312]]$description
#> [1] "advanced glycosylation end product-specific receptor [Source:HGNC Symbol;Acc:320]"
#> 
#> [[312]]$biotype
#> [1] "protein_coding"
#> 
#> [[312]]$feature_type
#> [1] "gene"
#> 
#> [[312]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[312]]$start
#> [1] 32148745
#> 
#> [[312]]$external_name
#> [1] "AGER"
#> 
#> [[312]]$canonical_transcript
#> [1] "ENST00000375076.4"
#> 
#> [[312]]$id
#> [1] "ENSG00000204305"
#> 
#> [[312]]$end
#> [1] 32152101
#> 
#> [[312]]$version
#> [1] 9
#> 
#> [[312]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[312]]$strand
#> [1] -1
#> 
#> [[312]]$gene_id
#> [1] "ENSG00000204305"
#> 
#> 
#> [[313]]
#> [[313]]$source
#> [1] "ensembl_havana"
#> 
#> [[313]]$seq_region_name
#> [1] "6"
#> 
#> [[313]]$description
#> [1] "protein phosphatase 1, regulatory subunit 18 [Source:HGNC Symbol;Acc:29413]"
#> 
#> [[313]]$biotype
#> [1] "protein_coding"
#> 
#> [[313]]$feature_type
#> [1] "gene"
#> 
#> [[313]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[313]]$start
#> [1] 30644166
#> 
#> [[313]]$canonical_transcript
#> [1] "ENST00000274853.3"
#> 
#> [[313]]$external_name
#> [1] "PPP1R18"
#> 
#> [[313]]$end
#> [1] 30655672
#> 
#> [[313]]$id
#> [1] "ENSG00000146112"
#> 
#> [[313]]$strand
#> [1] -1
#> 
#> [[313]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[313]]$version
#> [1] 7
#> 
#> [[313]]$gene_id
#> [1] "ENSG00000146112"
#> 
#> 
#> [[314]]
#> [[314]]$external_name
#> [1] "NRM"
#> 
#> [[314]]$canonical_transcript
#> [1] "ENST00000259953.4"
#> 
#> [[314]]$end
#> [1] 30659197
#> 
#> [[314]]$id
#> [1] "ENSG00000137404"
#> 
#> [[314]]$start
#> [1] 30655824
#> 
#> [[314]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[314]]$version
#> [1] 10
#> 
#> [[314]]$strand
#> [1] -1
#> 
#> [[314]]$gene_id
#> [1] "ENSG00000137404"
#> 
#> [[314]]$source
#> [1] "ensembl_havana"
#> 
#> [[314]]$seq_region_name
#> [1] "6"
#> 
#> [[314]]$biotype
#> [1] "protein_coding"
#> 
#> [[314]]$feature_type
#> [1] "gene"
#> 
#> [[314]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[314]]$description
#> [1] "nurim (nuclear envelope membrane protein) [Source:HGNC Symbol;Acc:8003]"
#> 
#> 
#> [[315]]
#> [[315]]$start
#> [1] 30667584
#> 
#> [[315]]$id
#> [1] "ENSG00000137337"
#> 
#> [[315]]$end
#> [1] 30685666
#> 
#> [[315]]$canonical_transcript
#> [1] "ENST00000376406.3"
#> 
#> [[315]]$external_name
#> [1] "MDC1"
#> 
#> [[315]]$gene_id
#> [1] "ENSG00000137337"
#> 
#> [[315]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[315]]$version
#> [1] 10
#> 
#> [[315]]$strand
#> [1] -1
#> 
#> [[315]]$seq_region_name
#> [1] "6"
#> 
#> [[315]]$source
#> [1] "ensembl_havana"
#> 
#> [[315]]$description
#> [1] "mediator of DNA-damage checkpoint 1 [Source:HGNC Symbol;Acc:21163]"
#> 
#> [[315]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[315]]$feature_type
#> [1] "gene"
#> 
#> [[315]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[316]]
#> [[316]]$source
#> [1] "ensembl_havana"
#> 
#> [[316]]$seq_region_name
#> [1] "6"
#> 
#> [[316]]$biotype
#> [1] "protein_coding"
#> 
#> [[316]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[316]]$feature_type
#> [1] "gene"
#> 
#> [[316]]$description
#> [1] "pre-B-cell leukemia homeobox 2 [Source:HGNC Symbol;Acc:8633]"
#> 
#> [[316]]$external_name
#> [1] "PBX2"
#> 
#> [[316]]$canonical_transcript
#> [1] "ENST00000375050.4"
#> 
#> [[316]]$end
#> [1] 32157963
#> 
#> [[316]]$id
#> [1] "ENSG00000204304"
#> 
#> [[316]]$start
#> [1] 32152512
#> 
#> [[316]]$strand
#> [1] -1
#> 
#> [[316]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[316]]$version
#> [1] 7
#> 
#> [[316]]$gene_id
#> [1] "ENSG00000204304"
#> 
#> 
#> [[317]]
#> [[317]]$gene_id
#> [1] "ENSG00000213654"
#> 
#> [[317]]$strand
#> [1] -1
#> 
#> [[317]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[317]]$version
#> [1] 5
#> 
#> [[317]]$id
#> [1] "ENSG00000213654"
#> 
#> [[317]]$end
#> [1] 32163300
#> 
#> [[317]]$external_name
#> [1] "GPSM3"
#> 
#> [[317]]$canonical_transcript
#> [1] "ENST00000375040.3"
#> 
#> [[317]]$start
#> [1] 32158543
#> 
#> [[317]]$feature_type
#> [1] "gene"
#> 
#> [[317]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[317]]$biotype
#> [1] "protein_coding"
#> 
#> [[317]]$description
#> [1] "G-protein signaling modulator 3 [Source:HGNC Symbol;Acc:13945]"
#> 
#> [[317]]$seq_region_name
#> [1] "6"
#> 
#> [[317]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[318]]
#> [[318]]$start
#> [1] 32162620
#> 
#> [[318]]$end
#> [1] 32191844
#> 
#> [[318]]$id
#> [1] "ENSG00000204301"
#> 
#> [[318]]$canonical_transcript
#> [1] "ENST00000375023.3"
#> 
#> [[318]]$external_name
#> [1] "NOTCH4"
#> 
#> [[318]]$gene_id
#> [1] "ENSG00000204301"
#> 
#> [[318]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[318]]$strand
#> [1] -1
#> 
#> [[318]]$version
#> [1] 5
#> 
#> [[318]]$seq_region_name
#> [1] "6"
#> 
#> [[318]]$source
#> [1] "ensembl_havana"
#> 
#> [[318]]$description
#> [1] "notch 4 [Source:HGNC Symbol;Acc:7884]"
#> 
#> [[318]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[318]]$feature_type
#> [1] "gene"
#> 
#> [[318]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[319]]
#> [[319]]$description
#> [1] "tubulin, beta class I [Source:HGNC Symbol;Acc:20778]"
#> 
#> [[319]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[319]]$feature_type
#> [1] "gene"
#> 
#> [[319]]$biotype
#> [1] "protein_coding"
#> 
#> [[319]]$seq_region_name
#> [1] "6"
#> 
#> [[319]]$source
#> [1] "ensembl_havana"
#> 
#> [[319]]$gene_id
#> [1] "ENSG00000196230"
#> 
#> [[319]]$strand
#> [1] 1
#> 
#> [[319]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[319]]$version
#> [1] 8
#> 
#> [[319]]$start
#> [1] 30687978
#> 
#> [[319]]$end
#> [1] 30693203
#> 
#> [[319]]$id
#> [1] "ENSG00000196230"
#> 
#> [[319]]$canonical_transcript
#> [1] "ENST00000327892.8"
#> 
#> [[319]]$external_name
#> [1] "TUBB"
#> 
#> 
#> [[320]]
#> [[320]]$strand
#> [1] -1
#> 
#> [[320]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[320]]$version
#> [1] 7
#> 
#> [[320]]$gene_id
#> [1] "ENSG00000204296"
#> 
#> [[320]]$external_name
#> [1] "C6orf10"
#> 
#> [[320]]$canonical_transcript
#> [1] "ENST00000447241.2"
#> 
#> [[320]]$end
#> [1] 32339684
#> 
#> [[320]]$id
#> [1] "ENSG00000204296"
#> 
#> [[320]]$start
#> [1] 32256303
#> 
#> [[320]]$biotype
#> [1] "protein_coding"
#> 
#> [[320]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[320]]$feature_type
#> [1] "gene"
#> 
#> [[320]]$description
#> [1] "chromosome 6 open reading frame 10 [Source:HGNC Symbol;Acc:13922]"
#> 
#> [[320]]$source
#> [1] "ensembl_havana"
#> 
#> [[320]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[321]]
#> [[321]]$external_name
#> [1] "BTNL2"
#> 
#> [[321]]$canonical_transcript
#> [1] "ENST00000454136.3"
#> 
#> [[321]]$id
#> [1] "ENSG00000204290"
#> 
#> [[321]]$end
#> [1] 32374905
#> 
#> [[321]]$start
#> [1] 32361740
#> 
#> [[321]]$strand
#> [1] -1
#> 
#> [[321]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[321]]$version
#> [1] 6
#> 
#> [[321]]$gene_id
#> [1] "ENSG00000204290"
#> 
#> [[321]]$source
#> [1] "ensembl_havana"
#> 
#> [[321]]$seq_region_name
#> [1] "6"
#> 
#> [[321]]$biotype
#> [1] "protein_coding"
#> 
#> [[321]]$feature_type
#> [1] "gene"
#> 
#> [[321]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[321]]$description
#> [1] "butyrophilin-like 2 (MHC class II associated) [Source:HGNC Symbol;Acc:1142]"
#> 
#> 
#> [[322]]
#> [[322]]$description
#> [1] "major histocompatibility complex, class II, DR alpha [Source:HGNC Symbol;Acc:4947]"
#> 
#> [[322]]$biotype
#> [1] "protein_coding"
#> 
#> [[322]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[322]]$feature_type
#> [1] "gene"
#> 
#> [[322]]$source
#> [1] "ensembl_havana"
#> 
#> [[322]]$seq_region_name
#> [1] "6"
#> 
#> [[322]]$strand
#> [1] 1
#> 
#> [[322]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[322]]$version
#> [1] 9
#> 
#> [[322]]$gene_id
#> [1] "ENSG00000204287"
#> 
#> [[322]]$start
#> [1] 32407619
#> 
#> [[322]]$canonical_transcript
#> [1] "ENST00000395388.2"
#> 
#> [[322]]$external_name
#> [1] "HLA-DRA"
#> 
#> [[322]]$end
#> [1] 32412823
#> 
#> [[322]]$id
#> [1] "ENSG00000204287"
#> 
#> 
#> [[323]]
#> [[323]]$source
#> [1] "ensembl_havana"
#> 
#> [[323]]$seq_region_name
#> [1] "6"
#> 
#> [[323]]$biotype
#> [1] "protein_coding"
#> 
#> [[323]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[323]]$feature_type
#> [1] "gene"
#> 
#> [[323]]$description
#> [1] "flotillin 1 [Source:HGNC Symbol;Acc:3757]"
#> 
#> [[323]]$canonical_transcript
#> [1] "ENST00000376389.3"
#> 
#> [[323]]$external_name
#> [1] "FLOT1"
#> 
#> [[323]]$end
#> [1] 30710510
#> 
#> [[323]]$id
#> [1] "ENSG00000137312"
#> 
#> [[323]]$start
#> [1] 30695486
#> 
#> [[323]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[323]]$strand
#> [1] -1
#> 
#> [[323]]$version
#> [1] 10
#> 
#> [[323]]$gene_id
#> [1] "ENSG00000137312"
#> 
#> 
#> [[324]]
#> [[324]]$feature_type
#> [1] "gene"
#> 
#> [[324]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[324]]$biotype
#> [1] "protein_coding"
#> 
#> [[324]]$description
#> [1] "major histocompatibility complex, class II, DR beta 5 [Source:HGNC Symbol;Acc:4953]"
#> 
#> [[324]]$seq_region_name
#> [1] "6"
#> 
#> [[324]]$source
#> [1] "ensembl_havana"
#> 
#> [[324]]$gene_id
#> [1] "ENSG00000198502"
#> 
#> [[324]]$version
#> [1] 5
#> 
#> [[324]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[324]]$strand
#> [1] -1
#> 
#> [[324]]$end
#> [1] 32498064
#> 
#> [[324]]$id
#> [1] "ENSG00000198502"
#> 
#> [[324]]$external_name
#> [1] "HLA-DRB5"
#> 
#> [[324]]$canonical_transcript
#> [1] "ENST00000374975.3"
#> 
#> [[324]]$start
#> [1] 32485120
#> 
#> 
#> [[325]]
#> [[325]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[325]]$version
#> [1] 6
#> 
#> [[325]]$strand
#> [1] -1
#> 
#> [[325]]$gene_id
#> [1] "ENSG00000196126"
#> 
#> [[325]]$canonical_transcript
#> [1] "ENST00000360004.5"
#> 
#> [[325]]$external_name
#> [1] "HLA-DRB1"
#> 
#> [[325]]$id
#> [1] "ENSG00000196126"
#> 
#> [[325]]$end
#> [1] 32557625
#> 
#> [[325]]$start
#> [1] 32546546
#> 
#> [[325]]$biotype
#> [1] "protein_coding"
#> 
#> [[325]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[325]]$feature_type
#> [1] "gene"
#> 
#> [[325]]$description
#> [1] "major histocompatibility complex, class II, DR beta 1 [Source:HGNC Symbol;Acc:4948]"
#> 
#> [[325]]$source
#> [1] "ensembl_havana"
#> 
#> [[325]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[326]]
#> [[326]]$description
#> [1] "major histocompatibility complex, class II, DQ alpha 1 [Source:HGNC Symbol;Acc:4942]"
#> 
#> [[326]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[326]]$feature_type
#> [1] "gene"
#> 
#> [[326]]$biotype
#> [1] "protein_coding"
#> 
#> [[326]]$seq_region_name
#> [1] "6"
#> 
#> [[326]]$source
#> [1] "ensembl_havana"
#> 
#> [[326]]$gene_id
#> [1] "ENSG00000196735"
#> 
#> [[326]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[326]]$strand
#> [1] 1
#> 
#> [[326]]$version
#> [1] 7
#> 
#> [[326]]$start
#> [1] 32595956
#> 
#> [[326]]$end
#> [1] 32614839
#> 
#> [[326]]$id
#> [1] "ENSG00000196735"
#> 
#> [[326]]$canonical_transcript
#> [1] "ENST00000343139.5"
#> 
#> [[326]]$external_name
#> [1] "HLA-DQA1"
#> 
#> 
#> [[327]]
#> [[327]]$gene_id
#> [1] "ENSG00000137331"
#> 
#> [[327]]$strand
#> [1] -1
#> 
#> [[327]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[327]]$version
#> [1] 11
#> 
#> [[327]]$end
#> [1] 30712331
#> 
#> [[327]]$id
#> [1] "ENSG00000137331"
#> 
#> [[327]]$external_name
#> [1] "IER3"
#> 
#> [[327]]$canonical_transcript
#> [1] "ENST00000259874.5"
#> 
#> [[327]]$start
#> [1] 30710976
#> 
#> [[327]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[327]]$feature_type
#> [1] "gene"
#> 
#> [[327]]$biotype
#> [1] "protein_coding"
#> 
#> [[327]]$description
#> [1] "immediate early response 3 [Source:HGNC Symbol;Acc:5392]"
#> 
#> [[327]]$seq_region_name
#> [1] "6"
#> 
#> [[327]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[328]]
#> [[328]]$start
#> [1] 32627244
#> 
#> [[328]]$end
#> [1] 32636160
#> 
#> [[328]]$id
#> [1] "ENSG00000179344"
#> 
#> [[328]]$canonical_transcript
#> [1] "ENST00000374943.4"
#> 
#> [[328]]$external_name
#> [1] "HLA-DQB1"
#> 
#> [[328]]$gene_id
#> [1] "ENSG00000179344"
#> 
#> [[328]]$version
#> [1] 12
#> 
#> [[328]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[328]]$strand
#> [1] -1
#> 
#> [[328]]$seq_region_name
#> [1] "6"
#> 
#> [[328]]$source
#> [1] "ensembl_havana"
#> 
#> [[328]]$description
#> [1] "major histocompatibility complex, class II, DQ beta 1 [Source:HGNC Symbol;Acc:4944]"
#> 
#> [[328]]$feature_type
#> [1] "gene"
#> 
#> [[328]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[328]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[329]]
#> [[329]]$source
#> [1] "ensembl_havana"
#> 
#> [[329]]$seq_region_name
#> [1] "6"
#> 
#> [[329]]$biotype
#> [1] "protein_coding"
#> 
#> [[329]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[329]]$feature_type
#> [1] "gene"
#> 
#> [[329]]$description
#> [1] "major histocompatibility complex, class II, DQ alpha 2 [Source:HGNC Symbol;Acc:4943]"
#> 
#> [[329]]$canonical_transcript
#> [1] "ENST00000374940.3"
#> 
#> [[329]]$external_name
#> [1] "HLA-DQA2"
#> 
#> [[329]]$end
#> [1] 32714992
#> 
#> [[329]]$id
#> [1] "ENSG00000237541"
#> 
#> [[329]]$start
#> [1] 32709119
#> 
#> [[329]]$version
#> [1] 3
#> 
#> [[329]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[329]]$strand
#> [1] 1
#> 
#> [[329]]$gene_id
#> [1] "ENSG00000237541"
#> 
#> 
#> [[330]]
#> [[330]]$source
#> [1] "ensembl_havana"
#> 
#> [[330]]$seq_region_name
#> [1] "6"
#> 
#> [[330]]$biotype
#> [1] "protein_coding"
#> 
#> [[330]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[330]]$feature_type
#> [1] "gene"
#> 
#> [[330]]$description
#> [1] "major histocompatibility complex, class II, DQ beta 2 [Source:HGNC Symbol;Acc:4945]"
#> 
#> [[330]]$canonical_transcript
#> [1] "ENST00000411527.1"
#> 
#> [[330]]$external_name
#> [1] "HLA-DQB2"
#> 
#> [[330]]$end
#> [1] 32731311
#> 
#> [[330]]$id
#> [1] "ENSG00000232629"
#> 
#> [[330]]$start
#> [1] 32723875
#> 
#> [[330]]$strand
#> [1] -1
#> 
#> [[330]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[330]]$version
#> [1] 4
#> 
#> [[330]]$gene_id
#> [1] "ENSG00000232629"
#> 
#> 
#> [[331]]
#> [[331]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[331]]$version
#> [1] 2
#> 
#> [[331]]$strand
#> [1] -1
#> 
#> [[331]]$gene_id
#> [1] "ENSG00000241106"
#> 
#> [[331]]$external_name
#> [1] "HLA-DOB"
#> 
#> [[331]]$canonical_transcript
#> [1] "ENST00000438763.2"
#> 
#> [[331]]$end
#> [1] 32784825
#> 
#> [[331]]$id
#> [1] "ENSG00000241106"
#> 
#> [[331]]$start
#> [1] 32780540
#> 
#> [[331]]$biotype
#> [1] "protein_coding"
#> 
#> [[331]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[331]]$feature_type
#> [1] "gene"
#> 
#> [[331]]$description
#> [1] "major histocompatibility complex, class II, DO beta [Source:HGNC Symbol;Acc:4937]"
#> 
#> [[331]]$source
#> [1] "ensembl_havana"
#> 
#> [[331]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[332]]
#> [[332]]$gene_id
#> [1] "ENSG00000204580"
#> 
#> [[332]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[332]]$version
#> [1] 7
#> 
#> [[332]]$strand
#> [1] 1
#> 
#> [[332]]$end
#> [1] 30867933
#> 
#> [[332]]$id
#> [1] "ENSG00000204580"
#> 
#> [[332]]$canonical_transcript
#> [1] "ENST00000376575.3"
#> 
#> [[332]]$external_name
#> [1] "DDR1"
#> 
#> [[332]]$start
#> [1] 30844198
#> 
#> [[332]]$feature_type
#> [1] "gene"
#> 
#> [[332]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[332]]$biotype
#> [1] "protein_coding"
#> 
#> [[332]]$description
#> [1] "discoidin domain receptor tyrosine kinase 1 [Source:HGNC Symbol;Acc:2730]"
#> 
#> [[332]]$seq_region_name
#> [1] "6"
#> 
#> [[332]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[333]]
#> [[333]]$biotype
#> [1] "protein_coding"
#> 
#> [[333]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[333]]$feature_type
#> [1] "gene"
#> 
#> [[333]]$description
#> [1] "Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:E7ENX8]"
#> 
#> [[333]]$source
#> [1] "havana"
#> 
#> [[333]]$seq_region_name
#> [1] "6"
#> 
#> [[333]]$version
#> [1] 1
#> 
#> [[333]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[333]]$strand
#> [1] -1
#> 
#> [[333]]$gene_id
#> [1] "ENSG00000250264"
#> 
#> [[333]]$canonical_transcript
#> [1] "ENST00000452392.2"
#> 
#> [[333]]$external_name
#> [1] "TAP2"
#> 
#> [[333]]$id
#> [1] "ENSG00000250264"
#> 
#> [[333]]$end
#> [1] 32806599
#> 
#> [[333]]$start
#> [1] 32781544
#> 
#> 
#> [[334]]
#> [[334]]$start
#> [1] 32789610
#> 
#> [[334]]$id
#> [1] "ENSG00000204267"
#> 
#> [[334]]$end
#> [1] 32806557
#> 
#> [[334]]$canonical_transcript
#> [1] "ENST00000374899.4"
#> 
#> [[334]]$external_name
#> [1] "TAP2"
#> 
#> [[334]]$gene_id
#> [1] "ENSG00000204267"
#> 
#> [[334]]$strand
#> [1] -1
#> 
#> [[334]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[334]]$version
#> [1] 9
#> 
#> [[334]]$seq_region_name
#> [1] "6"
#> 
#> [[334]]$source
#> [1] "ensembl_havana"
#> 
#> [[334]]$description
#> [1] "transporter 2, ATP-binding cassette, sub-family B (MDR/TAP) [Source:HGNC Symbol;Acc:44]"
#> 
#> [[334]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[334]]$feature_type
#> [1] "gene"
#> 
#> [[334]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[335]]
#> [[335]]$source
#> [1] "ensembl_havana"
#> 
#> [[335]]$seq_region_name
#> [1] "6"
#> 
#> [[335]]$description
#> [1] "proteasome (prosome, macropain) subunit, beta type, 8 [Source:HGNC Symbol;Acc:9545]"
#> 
#> [[335]]$biotype
#> [1] "protein_coding"
#> 
#> [[335]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[335]]$feature_type
#> [1] "gene"
#> 
#> [[335]]$start
#> [1] 32808494
#> 
#> [[335]]$external_name
#> [1] "PSMB8"
#> 
#> [[335]]$canonical_transcript
#> [1] "ENST00000374882.3"
#> 
#> [[335]]$end
#> [1] 32812480
#> 
#> [[335]]$id
#> [1] "ENSG00000204264"
#> 
#> [[335]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[335]]$version
#> [1] 4
#> 
#> [[335]]$strand
#> [1] -1
#> 
#> [[335]]$gene_id
#> [1] "ENSG00000204264"
#> 
#> 
#> [[336]]
#> [[336]]$end
#> [1] 32827362
#> 
#> [[336]]$id
#> [1] "ENSG00000240065"
#> 
#> [[336]]$canonical_transcript
#> [1] "ENST00000374859.2"
#> 
#> [[336]]$external_name
#> [1] "PSMB9"
#> 
#> [[336]]$start
#> [1] 32811913
#> 
#> [[336]]$gene_id
#> [1] "ENSG00000240065"
#> 
#> [[336]]$version
#> [1] 3
#> 
#> [[336]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[336]]$strand
#> [1] 1
#> 
#> [[336]]$seq_region_name
#> [1] "6"
#> 
#> [[336]]$source
#> [1] "ensembl_havana"
#> 
#> [[336]]$feature_type
#> [1] "gene"
#> 
#> [[336]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[336]]$biotype
#> [1] "protein_coding"
#> 
#> [[336]]$description
#> [1] "proteasome (prosome, macropain) subunit, beta type, 9 [Source:HGNC Symbol;Acc:9546]"
#> 
#> 
#> [[337]]
#> [[337]]$seq_region_name
#> [1] "6"
#> 
#> [[337]]$source
#> [1] "ensembl_havana"
#> 
#> [[337]]$feature_type
#> [1] "gene"
#> 
#> [[337]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[337]]$biotype
#> [1] "protein_coding"
#> 
#> [[337]]$description
#> [1] "transporter 1, ATP-binding cassette, sub-family B (MDR/TAP) [Source:HGNC Symbol;Acc:43]"
#> 
#> [[337]]$end
#> [1] 32821755
#> 
#> [[337]]$id
#> [1] "ENSG00000168394"
#> 
#> [[337]]$external_name
#> [1] "TAP1"
#> 
#> [[337]]$canonical_transcript
#> [1] "ENST00000354258.4"
#> 
#> [[337]]$start
#> [1] 32812986
#> 
#> [[337]]$gene_id
#> [1] "ENSG00000168394"
#> 
#> [[337]]$version
#> [1] 9
#> 
#> [[337]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[337]]$strand
#> [1] -1
#> 
#> 
#> [[338]]
#> [[338]]$gene_id
#> [1] "ENSG00000242574"
#> 
#> [[338]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[338]]$version
#> [1] 4
#> 
#> [[338]]$strand
#> [1] -1
#> 
#> [[338]]$end
#> [1] 32908847
#> 
#> [[338]]$id
#> [1] "ENSG00000242574"
#> 
#> [[338]]$external_name
#> [1] "HLA-DMB"
#> 
#> [[338]]$canonical_transcript
#> [1] "ENST00000418107.2"
#> 
#> [[338]]$start
#> [1] 32902406
#> 
#> [[338]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[338]]$feature_type
#> [1] "gene"
#> 
#> [[338]]$biotype
#> [1] "protein_coding"
#> 
#> [[338]]$description
#> [1] "major histocompatibility complex, class II, DM beta [Source:HGNC Symbol;Acc:4935]"
#> 
#> [[338]]$seq_region_name
#> [1] "6"
#> 
#> [[338]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[339]]
#> [[339]]$description
#> [1] "Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:F6UB75]"
#> 
#> [[339]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[339]]$feature_type
#> [1] "gene"
#> 
#> [[339]]$biotype
#> [1] "protein_coding"
#> 
#> [[339]]$seq_region_name
#> [1] "6"
#> 
#> [[339]]$source
#> [1] "havana"
#> 
#> [[339]]$gene_id
#> [1] "ENSG00000248993"
#> 
#> [[339]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[339]]$version
#> [1] 1
#> 
#> [[339]]$strand
#> [1] -1
#> 
#> [[339]]$start
#> [1] 32905141
#> 
#> [[339]]$id
#> [1] "ENSG00000248993"
#> 
#> [[339]]$end
#> [1] 32920899
#> 
#> [[339]]$canonical_transcript
#> [1] "ENST00000429234.1"
#> 
#> [[339]]$external_name
#> [1] "XXbac-BPG181M17.5"
#> 
#> 
#> [[340]]
#> [[340]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[340]]$strand
#> [1] -1
#> 
#> [[340]]$version
#> [1] 10
#> 
#> [[340]]$gene_id
#> [1] "ENSG00000204257"
#> 
#> [[340]]$start
#> [1] 32916390
#> 
#> [[340]]$canonical_transcript
#> [1] "ENST00000374843.4"
#> 
#> [[340]]$external_name
#> [1] "HLA-DMA"
#> 
#> [[340]]$end
#> [1] 32936871
#> 
#> [[340]]$id
#> [1] "ENSG00000204257"
#> 
#> [[340]]$description
#> [1] "major histocompatibility complex, class II, DM alpha [Source:HGNC Symbol;Acc:4934]"
#> 
#> [[340]]$biotype
#> [1] "protein_coding"
#> 
#> [[340]]$feature_type
#> [1] "gene"
#> 
#> [[340]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[340]]$source
#> [1] "ensembl_havana"
#> 
#> [[340]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[341]]
#> [[341]]$biotype
#> [1] "protein_coding"
#> 
#> [[341]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[341]]$feature_type
#> [1] "gene"
#> 
#> [[341]]$description
#> [1] "bromodomain containing 2 [Source:HGNC Symbol;Acc:1103]"
#> 
#> [[341]]$source
#> [1] "ensembl_havana"
#> 
#> [[341]]$seq_region_name
#> [1] "6"
#> 
#> [[341]]$strand
#> [1] 1
#> 
#> [[341]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[341]]$version
#> [1] 8
#> 
#> [[341]]$gene_id
#> [1] "ENSG00000204256"
#> 
#> [[341]]$canonical_transcript
#> [1] "ENST00000395289.2"
#> 
#> [[341]]$external_name
#> [1] "BRD2"
#> 
#> [[341]]$end
#> [1] 32949282
#> 
#> [[341]]$id
#> [1] "ENSG00000204256"
#> 
#> [[341]]$start
#> [1] 32936437
#> 
#> 
#> [[342]]
#> [[342]]$gene_id
#> [1] "ENSG00000204252"
#> 
#> [[342]]$strand
#> [1] -1
#> 
#> [[342]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[342]]$version
#> [1] 8
#> 
#> [[342]]$id
#> [1] "ENSG00000204252"
#> 
#> [[342]]$end
#> [1] 32977389
#> 
#> [[342]]$external_name
#> [1] "HLA-DOA"
#> 
#> [[342]]$canonical_transcript
#> [1] "ENST00000229829.5"
#> 
#> [[342]]$start
#> [1] 32971955
#> 
#> [[342]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[342]]$feature_type
#> [1] "gene"
#> 
#> [[342]]$biotype
#> [1] "protein_coding"
#> 
#> [[342]]$description
#> [1] "major histocompatibility complex, class II, DO alpha [Source:HGNC Symbol;Acc:4936]"
#> 
#> [[342]]$seq_region_name
#> [1] "6"
#> 
#> [[342]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[343]]
#> [[343]]$strand
#> [1] -1
#> 
#> [[343]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[343]]$version
#> [1] 3
#> 
#> [[343]]$gene_id
#> [1] "ENSG00000231389"
#> 
#> [[343]]$canonical_transcript
#> [1] "ENST00000419277.1"
#> 
#> [[343]]$external_name
#> [1] "HLA-DPA1"
#> 
#> [[343]]$id
#> [1] "ENSG00000231389"
#> 
#> [[343]]$end
#> [1] 33048552
#> 
#> [[343]]$start
#> [1] 33032346
#> 
#> [[343]]$biotype
#> [1] "protein_coding"
#> 
#> [[343]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[343]]$feature_type
#> [1] "gene"
#> 
#> [[343]]$description
#> [1] "major histocompatibility complex, class II, DP alpha 1 [Source:HGNC Symbol;Acc:4938]"
#> 
#> [[343]]$source
#> [1] "ensembl_havana"
#> 
#> [[343]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[344]]
#> [[344]]$end
#> [1] 33054978
#> 
#> [[344]]$id
#> [1] "ENSG00000223865"
#> 
#> [[344]]$external_name
#> [1] "HLA-DPB1"
#> 
#> [[344]]$canonical_transcript
#> [1] "ENST00000418931.2"
#> 
#> [[344]]$start
#> [1] 33043703
#> 
#> [[344]]$gene_id
#> [1] "ENSG00000223865"
#> 
#> [[344]]$strand
#> [1] 1
#> 
#> [[344]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[344]]$version
#> [1] 6
#> 
#> [[344]]$seq_region_name
#> [1] "6"
#> 
#> [[344]]$source
#> [1] "ensembl_havana"
#> 
#> [[344]]$feature_type
#> [1] "gene"
#> 
#> [[344]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[344]]$biotype
#> [1] "protein_coding"
#> 
#> [[344]]$description
#> [1] "major histocompatibility complex, class II, DP beta 1 [Source:HGNC Symbol;Acc:4940]"
#> 
#> 
#> [[345]]
#> [[345]]$end
#> [1] 30881883
#> 
#> [[345]]$id
#> [1] "ENSG00000213780"
#> 
#> [[345]]$canonical_transcript
#> [1] "ENST00000259895.4"
#> 
#> [[345]]$external_name
#> [1] "GTF2H4"
#> 
#> [[345]]$start
#> [1] 30875961
#> 
#> [[345]]$gene_id
#> [1] "ENSG00000213780"
#> 
#> [[345]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[345]]$version
#> [1] 6
#> 
#> [[345]]$strand
#> [1] 1
#> 
#> [[345]]$seq_region_name
#> [1] "6"
#> 
#> [[345]]$source
#> [1] "ensembl_havana"
#> 
#> [[345]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[345]]$feature_type
#> [1] "gene"
#> 
#> [[345]]$biotype
#> [1] "protein_coding"
#> 
#> [[345]]$description
#> [1] "general transcription factor IIH, polypeptide 4, 52kDa [Source:HGNC Symbol;Acc:4658]"
#> 
#> 
#> [[346]]
#> [[346]]$start
#> [1] 30876019
#> 
#> [[346]]$end
#> [1] 30894236
#> 
#> [[346]]$id
#> [1] "ENSG00000137411"
#> 
#> [[346]]$canonical_transcript
#> [1] "ENST00000541562.1"
#> 
#> [[346]]$external_name
#> [1] "VARS2"
#> 
#> [[346]]$gene_id
#> [1] "ENSG00000137411"
#> 
#> [[346]]$version
#> [1] 12
#> 
#> [[346]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[346]]$strand
#> [1] 1
#> 
#> [[346]]$seq_region_name
#> [1] "6"
#> 
#> [[346]]$source
#> [1] "ensembl_havana"
#> 
#> [[346]]$description
#> [1] "valyl-tRNA synthetase 2, mitochondrial [Source:HGNC Symbol;Acc:21642]"
#> 
#> [[346]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[346]]$feature_type
#> [1] "gene"
#> 
#> [[346]]$biotype
#> [1] "protein_coding"
#> 
#> 
#> [[347]]
#> [[347]]$id
#> [1] "ENSG00000204248"
#> 
#> [[347]]$end
#> [1] 33160276
#> 
#> [[347]]$canonical_transcript
#> [1] "ENST00000374708.4"
#> 
#> [[347]]$external_name
#> [1] "COL11A2"
#> 
#> [[347]]$start
#> [1] 33130458
#> 
#> [[347]]$gene_id
#> [1] "ENSG00000204248"
#> 
#> [[347]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[347]]$version
#> [1] 6
#> 
#> [[347]]$strand
#> [1] -1
#> 
#> [[347]]$seq_region_name
#> [1] "6"
#> 
#> [[347]]$source
#> [1] "ensembl_havana"
#> 
#> [[347]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[347]]$feature_type
#> [1] "gene"
#> 
#> [[347]]$biotype
#> [1] "protein_coding"
#> 
#> [[347]]$description
#> [1] "collagen, type XI, alpha 2 [Source:HGNC Symbol;Acc:2187]"
#> 
#> 
#> [[348]]
#> [[348]]$id
#> [1] "ENSG00000196260"
#> 
#> [[348]]$end
#> [1] 30899952
#> 
#> [[348]]$external_name
#> [1] "SFTA2"
#> 
#> [[348]]$canonical_transcript
#> [1] "ENST00000359086.3"
#> 
#> [[348]]$start
#> [1] 30899130
#> 
#> [[348]]$gene_id
#> [1] "ENSG00000196260"
#> 
#> [[348]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[348]]$strand
#> [1] -1
#> 
#> [[348]]$version
#> [1] 3
#> 
#> [[348]]$seq_region_name
#> [1] "6"
#> 
#> [[348]]$source
#> [1] "ensembl_havana"
#> 
#> [[348]]$feature_type
#> [1] "gene"
#> 
#> [[348]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[348]]$biotype
#> [1] "protein_coding"
#> 
#> [[348]]$description
#> [1] "surfactant associated 2 [Source:HGNC Symbol;Acc:18386]"
#> 
#> 
#> [[349]]
#> [[349]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[349]]$strand
#> [1] 1
#> 
#> [[349]]$version
#> [1] 7
#> 
#> [[349]]$gene_id
#> [1] "ENSG00000168631"
#> 
#> [[349]]$start
#> [1] 30908749
#> 
#> [[349]]$external_name
#> [1] "DPCR1"
#> 
#> [[349]]$canonical_transcript
#> [1] "ENST00000462446.1"
#> 
#> [[349]]$end
#> [1] 30921998
#> 
#> [[349]]$id
#> [1] "ENSG00000168631"
#> 
#> [[349]]$description
#> [1] "diffuse panbronchiolitis critical region 1 [Source:HGNC Symbol;Acc:21666]"
#> 
#> [[349]]$biotype
#> [1] "protein_coding"
#> 
#> [[349]]$feature_type
#> [1] "gene"
#> 
#> [[349]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[349]]$source
#> [1] "ensembl_havana"
#> 
#> [[349]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[350]]
#> [[350]]$seq_region_name
#> [1] "6"
#> 
#> [[350]]$source
#> [1] "ensembl_havana"
#> 
#> [[350]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[350]]$feature_type
#> [1] "gene"
#> 
#> [[350]]$biotype
#> [1] "protein_coding"
#> 
#> [[350]]$description
#> [1] "mucin 21, cell surface associated [Source:HGNC Symbol;Acc:21661]"
#> 
#> [[350]]$end
#> [1] 30957680
#> 
#> [[350]]$id
#> [1] "ENSG00000204544"
#> 
#> [[350]]$external_name
#> [1] "MUC21"
#> 
#> [[350]]$canonical_transcript
#> [1] "ENST00000376296.3"
#> 
#> [[350]]$start
#> [1] 30951495
#> 
#> [[350]]$gene_id
#> [1] "ENSG00000204544"
#> 
#> [[350]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[350]]$strand
#> [1] 1
#> 
#> [[350]]$version
#> [1] 5
#> 
#> 
#> [[351]]
#> [[351]]$seq_region_name
#> [1] "6"
#> 
#> [[351]]$source
#> [1] "havana"
#> 
#> [[351]]$description
#> [1] "mucin 22 [Source:HGNC Symbol;Acc:39755]"
#> 
#> [[351]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[351]]$feature_type
#> [1] "gene"
#> 
#> [[351]]$biotype
#> [1] "protein_coding"
#> 
#> [[351]]$start
#> [1] 30978251
#> 
#> [[351]]$id
#> [1] "ENSG00000261272"
#> 
#> [[351]]$end
#> [1] 31003179
#> 
#> [[351]]$canonical_transcript
#> [1] "ENST00000561890.1"
#> 
#> [[351]]$external_name
#> [1] "MUC22"
#> 
#> [[351]]$gene_id
#> [1] "ENSG00000261272"
#> 
#> [[351]]$strand
#> [1] 1
#> 
#> [[351]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[351]]$version
#> [1] 1
#> 
#> 
#> [[352]]
#> [[352]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[352]]$strand
#> [1] -1
#> 
#> [[352]]$version
#> [1] 2
#> 
#> [[352]]$gene_id
#> [1] "ENSG00000204542"
#> 
#> [[352]]$start
#> [1] 31079000
#> 
#> [[352]]$canonical_transcript
#> [1] "ENST00000259870.3"
#> 
#> [[352]]$external_name
#> [1] "C6orf15"
#> 
#> [[352]]$end
#> [1] 31080336
#> 
#> [[352]]$id
#> [1] "ENSG00000204542"
#> 
#> [[352]]$description
#> [1] "chromosome 6 open reading frame 15 [Source:HGNC Symbol;Acc:13927]"
#> 
#> [[352]]$biotype
#> [1] "protein_coding"
#> 
#> [[352]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[352]]$feature_type
#> [1] "gene"
#> 
#> [[352]]$source
#> [1] "ensembl_havana"
#> 
#> [[352]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[353]]
#> [[353]]$start
#> [1] 31082527
#> 
#> [[353]]$external_name
#> [1] "PSORS1C1"
#> 
#> [[353]]$canonical_transcript
#> [1] "ENST00000259881.9"
#> 
#> [[353]]$end
#> [1] 31107869
#> 
#> [[353]]$id
#> [1] "ENSG00000204540"
#> 
#> [[353]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[353]]$version
#> [1] 6
#> 
#> [[353]]$strand
#> [1] 1
#> 
#> [[353]]$gene_id
#> [1] "ENSG00000204540"
#> 
#> [[353]]$source
#> [1] "ensembl_havana"
#> 
#> [[353]]$seq_region_name
#> [1] "6"
#> 
#> [[353]]$description
#> [1] "psoriasis susceptibility 1 candidate 1 [Source:HGNC Symbol;Acc:17202]"
#> 
#> [[353]]$biotype
#> [1] "protein_coding"
#> 
#> [[353]]$feature_type
#> [1] "gene"
#> 
#> [[353]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> 
#> [[354]]
#> [[354]]$description
#> [1] "corneodesmosin [Source:HGNC Symbol;Acc:1802]"
#> 
#> [[354]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[354]]$feature_type
#> [1] "gene"
#> 
#> [[354]]$biotype
#> [1] "protein_coding"
#> 
#> [[354]]$seq_region_name
#> [1] "6"
#> 
#> [[354]]$source
#> [1] "ensembl_havana"
#> 
#> [[354]]$gene_id
#> [1] "ENSG00000204539"
#> 
#> [[354]]$version
#> [1] 3
#> 
#> [[354]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[354]]$strand
#> [1] -1
#> 
#> [[354]]$start
#> [1] 31082867
#> 
#> [[354]]$id
#> [1] "ENSG00000204539"
#> 
#> [[354]]$end
#> [1] 31088223
#> 
#> [[354]]$canonical_transcript
#> [1] "ENST00000376288.2"
#> 
#> [[354]]$external_name
#> [1] "CDSN"
#> 
#> 
#> [[355]]
#> [[355]]$biotype
#> [1] "protein_coding"
#> 
#> [[355]]$feature_type
#> [1] "gene"
#> 
#> [[355]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[355]]$description
#> [1] "psoriasis susceptibility 1 candidate 2 [Source:HGNC Symbol;Acc:17199]"
#> 
#> [[355]]$source
#> [1] "ensembl_havana"
#> 
#> [[355]]$seq_region_name
#> [1] "6"
#> 
#> [[355]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[355]]$version
#> [1] 3
#> 
#> [[355]]$strand
#> [1] -1
#> 
#> [[355]]$gene_id
#> [1] "ENSG00000204538"
#> 
#> [[355]]$canonical_transcript
#> [1] "ENST00000259845.4"
#> 
#> [[355]]$external_name
#> [1] "PSORS1C2"
#> 
#> [[355]]$end
#> [1] 31107127
#> 
#> [[355]]$id
#> [1] "ENSG00000204538"
#> 
#> [[355]]$start
#> [1] 31105313
#> 
#> 
#> [[356]]
#> [[356]]$source
#> [1] "ensembl_havana"
#> 
#> [[356]]$seq_region_name
#> [1] "6"
#> 
#> [[356]]$biotype
#> [1] "protein_coding"
#> 
#> [[356]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[356]]$feature_type
#> [1] "gene"
#> 
#> [[356]]$description
#> [1] "retinoid X receptor, beta [Source:HGNC Symbol;Acc:10478]"
#> 
#> [[356]]$canonical_transcript
#> [1] "ENST00000374685.4"
#> 
#> [[356]]$external_name
#> [1] "RXRB"
#> 
#> [[356]]$id
#> [1] "ENSG00000204231"
#> 
#> [[356]]$end
#> [1] 33168630
#> 
#> [[356]]$start
#> [1] 33161365
#> 
#> [[356]]$strand
#> [1] -1
#> 
#> [[356]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[356]]$version
#> [1] 6
#> 
#> [[356]]$gene_id
#> [1] "ENSG00000204231"
#> 
#> 
#> [[357]]
#> [[357]]$strand
#> [1] -1
#> 
#> [[357]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[357]]$version
#> [1] 9
#> 
#> [[357]]$gene_id
#> [1] "ENSG00000204536"
#> 
#> [[357]]$start
#> [1] 31110216
#> 
#> [[357]]$external_name
#> [1] "CCHCR1"
#> 
#> [[357]]$canonical_transcript
#> [1] "ENST00000396268.3"
#> 
#> [[357]]$id
#> [1] "ENSG00000204536"
#> 
#> [[357]]$end
#> [1] 31126015
#> 
#> [[357]]$description
#> [1] "coiled-coil alpha-helical rod protein 1 [Source:HGNC Symbol;Acc:13930]"
#> 
#> [[357]]$biotype
#> [1] "protein_coding"
#> 
#> [[357]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[357]]$feature_type
#> [1] "gene"
#> 
#> [[357]]$source
#> [1] "ensembl_havana"
#> 
#> [[357]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[358]]
#> [[358]]$gene_id
#> [1] "ENSG00000112473"
#> 
#> [[358]]$strand
#> [1] 1
#> 
#> [[358]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[358]]$version
#> [1] 12
#> 
#> [[358]]$end
#> [1] 33172216
#> 
#> [[358]]$id
#> [1] "ENSG00000112473"
#> 
#> [[358]]$external_name
#> [1] "SLC39A7"
#> 
#> [[358]]$canonical_transcript
#> [1] "ENST00000374677.3"
#> 
#> [[358]]$start
#> [1] 33168222
#> 
#> [[358]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[358]]$feature_type
#> [1] "gene"
#> 
#> [[358]]$biotype
#> [1] "protein_coding"
#> 
#> [[358]]$description
#> [1] "solute carrier family 39 (zinc transporter), member 7 [Source:HGNC Symbol;Acc:4927]"
#> 
#> [[358]]$seq_region_name
#> [1] "6"
#> 
#> [[358]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[359]]
#> [[359]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[359]]$version
#> [1] 7
#> 
#> [[359]]$strand
#> [1] 1
#> 
#> [[359]]$gene_id
#> [1] "ENSG00000137310"
#> 
#> [[359]]$external_name
#> [1] "TCF19"
#> 
#> [[359]]$canonical_transcript
#> [1] "ENST00000376257.3"
#> 
#> [[359]]$id
#> [1] "ENSG00000137310"
#> 
#> [[359]]$end
#> [1] 31134936
#> 
#> [[359]]$start
#> [1] 31126319
#> 
#> [[359]]$biotype
#> [1] "protein_coding"
#> 
#> [[359]]$feature_type
#> [1] "gene"
#> 
#> [[359]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[359]]$description
#> [1] "transcription factor 19 [Source:HGNC Symbol;Acc:11629]"
#> 
#> [[359]]$source
#> [1] "ensembl_havana"
#> 
#> [[359]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[360]]
#> [[360]]$end
#> [1] 31148508
#> 
#> [[360]]$id
#> [1] "ENSG00000204531"
#> 
#> [[360]]$external_name
#> [1] "POU5F1"
#> 
#> [[360]]$canonical_transcript
#> [1] "ENST00000259915.8"
#> 
#> [[360]]$start
#> [1] 31132119
#> 
#> [[360]]$gene_id
#> [1] "ENSG00000204531"
#> 
#> [[360]]$version
#> [1] 11
#> 
#> [[360]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[360]]$strand
#> [1] -1
#> 
#> [[360]]$seq_region_name
#> [1] "6"
#> 
#> [[360]]$source
#> [1] "ensembl_havana"
#> 
#> [[360]]$feature_type
#> [1] "gene"
#> 
#> [[360]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[360]]$biotype
#> [1] "protein_coding"
#> 
#> [[360]]$description
#> [1] "POU class 5 homeobox 1 [Source:HGNC Symbol;Acc:9221]"
#> 
#> 
#> [[361]]
#> [[361]]$source
#> [1] "havana"
#> 
#> [[361]]$seq_region_name
#> [1] "6"
#> 
#> [[361]]$biotype
#> [1] "protein_coding"
#> 
#> [[361]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[361]]$feature_type
#> [1] "gene"
#> 
#> [[361]]$description
#> [1] "HLA complex group 27 (non-protein coding) [Source:HGNC Symbol;Acc:27366]"
#> 
#> [[361]]$canonical_transcript
#> [1] "ENST00000383331.4"
#> 
#> [[361]]$external_name
#> [1] "HCG27"
#> 
#> [[361]]$id
#> [1] "ENSG00000206344"
#> 
#> [[361]]$end
#> [1] 31171745
#> 
#> [[361]]$start
#> [1] 31165537
#> 
#> [[361]]$version
#> [1] 6
#> 
#> [[361]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[361]]$strand
#> [1] 1
#> 
#> [[361]]$gene_id
#> [1] "ENSG00000206344"
#> 
#> 
#> [[362]]
#> [[362]]$seq_region_name
#> [1] "6"
#> 
#> [[362]]$source
#> [1] "ensembl_havana"
#> 
#> [[362]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[362]]$feature_type
#> [1] "gene"
#> 
#> [[362]]$biotype
#> [1] "protein_coding"
#> 
#> [[362]]$description
#> [1] "major histocompatibility complex, class I, C [Source:HGNC Symbol;Acc:4933]"
#> 
#> [[362]]$id
#> [1] "ENSG00000204525"
#> 
#> [[362]]$end
#> [1] 31239907
#> 
#> [[362]]$canonical_transcript
#> [1] "ENST00000376228.5"
#> 
#> [[362]]$external_name
#> [1] "HLA-C"
#> 
#> [[362]]$start
#> [1] 31236526
#> 
#> [[362]]$gene_id
#> [1] "ENSG00000204525"
#> 
#> [[362]]$strand
#> [1] -1
#> 
#> [[362]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[362]]$version
#> [1] 10
#> 
#> 
#> [[363]]
#> [[363]]$end
#> [1] 31324965
#> 
#> [[363]]$id
#> [1] "ENSG00000234745"
#> 
#> [[363]]$canonical_transcript
#> [1] "ENST00000412585.2"
#> 
#> [[363]]$external_name
#> [1] "HLA-B"
#> 
#> [[363]]$start
#> [1] 31321649
#> 
#> [[363]]$gene_id
#> [1] "ENSG00000234745"
#> 
#> [[363]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[363]]$strand
#> [1] -1
#> 
#> [[363]]$version
#> [1] 5
#> 
#> [[363]]$seq_region_name
#> [1] "6"
#> 
#> [[363]]$source
#> [1] "ensembl_havana"
#> 
#> [[363]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[363]]$feature_type
#> [1] "gene"
#> 
#> [[363]]$biotype
#> [1] "protein_coding"
#> 
#> [[363]]$description
#> [1] "major histocompatibility complex, class I, B [Source:HGNC Symbol;Acc:4932]"
#> 
#> 
#> [[364]]
#> [[364]]$description
#> [1] "MHC class I polypeptide-related sequence A [Source:HGNC Symbol;Acc:7090]"
#> 
#> [[364]]$biotype
#> [1] "protein_coding"
#> 
#> [[364]]$feature_type
#> [1] "gene"
#> 
#> [[364]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[364]]$source
#> [1] "ensembl_havana"
#> 
#> [[364]]$seq_region_name
#> [1] "6"
#> 
#> [[364]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[364]]$version
#> [1] 8
#> 
#> [[364]]$strand
#> [1] 1
#> 
#> [[364]]$gene_id
#> [1] "ENSG00000204520"
#> 
#> [[364]]$start
#> [1] 31371356
#> 
#> [[364]]$canonical_transcript
#> [1] "ENST00000449934.2"
#> 
#> [[364]]$external_name
#> [1] "MICA"
#> 
#> [[364]]$end
#> [1] 31383092
#> 
#> [[364]]$id
#> [1] "ENSG00000204520"
#> 
#> 
#> [[365]]
#> [[365]]$strand
#> [1] 1
#> 
#> [[365]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[365]]$version
#> [1] 5
#> 
#> [[365]]$gene_id
#> [1] "ENSG00000204516"
#> 
#> [[365]]$canonical_transcript
#> [1] "ENST00000252229.6"
#> 
#> [[365]]$external_name
#> [1] "MICB"
#> 
#> [[365]]$end
#> [1] 31478901
#> 
#> [[365]]$id
#> [1] "ENSG00000204516"
#> 
#> [[365]]$start
#> [1] 31462658
#> 
#> [[365]]$biotype
#> [1] "protein_coding"
#> 
#> [[365]]$feature_type
#> [1] "gene"
#> 
#> [[365]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[365]]$description
#> [1] "MHC class I polypeptide-related sequence B [Source:HGNC Symbol;Acc:7091]"
#> 
#> [[365]]$source
#> [1] "ensembl_havana"
#> 
#> [[365]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[366]]
#> [[366]]$seq_region_name
#> [1] "6"
#> 
#> [[366]]$source
#> [1] "ensembl_havana"
#> 
#> [[366]]$description
#> [1] "mitochondrial coiled-coil domain 1 [Source:HGNC Symbol;Acc:20668]"
#> 
#> [[366]]$feature_type
#> [1] "gene"
#> 
#> [[366]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[366]]$biotype
#> [1] "protein_coding"
#> 
#> [[366]]$start
#> [1] 31496494
#> 
#> [[366]]$id
#> [1] "ENSG00000204511"
#> 
#> [[366]]$end
#> [1] 31498009
#> 
#> [[366]]$external_name
#> [1] "MCCD1"
#> 
#> [[366]]$canonical_transcript
#> [1] "ENST00000376191.2"
#> 
#> [[366]]$gene_id
#> [1] "ENSG00000204511"
#> 
#> [[366]]$version
#> [1] 2
#> 
#> [[366]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[366]]$strand
#> [1] 1
#> 
#> 
#> [[367]]
#> [[367]]$source
#> [1] "havana"
#> 
#> [[367]]$seq_region_name
#> [1] "6"
#> 
#> [[367]]$description
#> [1] "ATP6V1G2-DDX39B readthrough (NMD candidate) [Source:HGNC Symbol;Acc:41999]"
#> 
#> [[367]]$biotype
#> [1] "protein_coding"
#> 
#> [[367]]$logic_name
#> [1] "havana_homo_sapiens_37"
#> 
#> [[367]]$feature_type
#> [1] "gene"
#> 
#> [[367]]$start
#> [1] 31497996
#> 
#> [[367]]$external_name
#> [1] "ATP6V1G2-DDX39B"
#> 
#> [[367]]$canonical_transcript
#> [1] "ENST00000376185.1"
#> 
#> [[367]]$end
#> [1] 31514385
#> 
#> [[367]]$id
#> [1] "ENSG00000254870"
#> 
#> [[367]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[367]]$version
#> [1] 1
#> 
#> [[367]]$strand
#> [1] -1
#> 
#> [[367]]$gene_id
#> [1] "ENSG00000254870"
#> 
#> 
#> [[368]]
#> [[368]]$canonical_transcript
#> [1] "ENST00000396172.1"
#> 
#> [[368]]$external_name
#> [1] "DDX39B"
#> 
#> [[368]]$end
#> [1] 31510225
#> 
#> [[368]]$id
#> [1] "ENSG00000198563"
#> 
#> [[368]]$start
#> [1] 31497996
#> 
#> [[368]]$strand
#> [1] -1
#> 
#> [[368]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[368]]$version
#> [1] 9
#> 
#> [[368]]$gene_id
#> [1] "ENSG00000198563"
#> 
#> [[368]]$source
#> [1] "ensembl_havana"
#> 
#> [[368]]$seq_region_name
#> [1] "6"
#> 
#> [[368]]$biotype
#> [1] "protein_coding"
#> 
#> [[368]]$feature_type
#> [1] "gene"
#> 
#> [[368]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[368]]$description
#> [1] "DEAD (Asp-Glu-Ala-Asp) box polypeptide 39B [Source:HGNC Symbol;Acc:13917]"
#> 
#> 
#> [[369]]
#> [[369]]$version
#> [1] 6
#> 
#> [[369]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[369]]$strand
#> [1] -1
#> 
#> [[369]]$gene_id
#> [1] "ENSG00000213760"
#> 
#> [[369]]$canonical_transcript
#> [1] "ENST00000303892.5"
#> 
#> [[369]]$external_name
#> [1] "ATP6V1G2"
#> 
#> [[369]]$id
#> [1] "ENSG00000213760"
#> 
#> [[369]]$end
#> [1] 31516204
#> 
#> [[369]]$start
#> [1] 31512239
#> 
#> [[369]]$biotype
#> [1] "protein_coding"
#> 
#> [[369]]$feature_type
#> [1] "gene"
#> 
#> [[369]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[369]]$description
#> [1] "ATPase, H+ transporting, lysosomal 13kDa, V1 subunit G2 [Source:HGNC Symbol;Acc:862]"
#> 
#> [[369]]$source
#> [1] "ensembl_havana"
#> 
#> [[369]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[370]]
#> [[370]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[370]]$version
#> [1] 6
#> 
#> [[370]]$strand
#> [1] 1
#> 
#> [[370]]$gene_id
#> [1] "ENSG00000204498"
#> 
#> [[370]]$canonical_transcript
#> [1] "ENST00000376148.4"
#> 
#> [[370]]$external_name
#> [1] "NFKBIL1"
#> 
#> [[370]]$end
#> [1] 31526606
#> 
#> [[370]]$id
#> [1] "ENSG00000204498"
#> 
#> [[370]]$start
#> [1] 31514647
#> 
#> [[370]]$biotype
#> [1] "protein_coding"
#> 
#> [[370]]$feature_type
#> [1] "gene"
#> 
#> [[370]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[370]]$description
#> [1] "nuclear factor of kappa light polypeptide gene enhancer in B-cells inhibitor-like 1 [Source:HGNC Symbol;Acc:7800]"
#> 
#> [[370]]$source
#> [1] "ensembl_havana"
#> 
#> [[370]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[371]]
#> [[371]]$seq_region_name
#> [1] "6"
#> 
#> [[371]]$source
#> [1] "ensembl_havana"
#> 
#> [[371]]$description
#> [1] "lymphotoxin alpha [Source:HGNC Symbol;Acc:6709]"
#> 
#> [[371]]$feature_type
#> [1] "gene"
#> 
#> [[371]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[371]]$biotype
#> [1] "protein_coding"
#> 
#> [[371]]$start
#> [1] 31539831
#> 
#> [[371]]$id
#> [1] "ENSG00000226979"
#> 
#> [[371]]$end
#> [1] 31542101
#> 
#> [[371]]$external_name
#> [1] "LTA"
#> 
#> [[371]]$canonical_transcript
#> [1] "ENST00000454783.1"
#> 
#> [[371]]$gene_id
#> [1] "ENSG00000226979"
#> 
#> [[371]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[371]]$version
#> [1] 4
#> 
#> [[371]]$strand
#> [1] 1
#> 
#> 
#> [[372]]
#> [[372]]$source
#> [1] "ensembl_havana"
#> 
#> [[372]]$seq_region_name
#> [1] "6"
#> 
#> [[372]]$description
#> [1] "tumor necrosis factor [Source:HGNC Symbol;Acc:11892]"
#> 
#> [[372]]$biotype
#> [1] "protein_coding"
#> 
#> [[372]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[372]]$feature_type
#> [1] "gene"
#> 
#> [[372]]$start
#> [1] 31543344
#> 
#> [[372]]$external_name
#> [1] "TNF"
#> 
#> [[372]]$canonical_transcript
#> [1] "ENST00000449264.2"
#> 
#> [[372]]$end
#> [1] 31546113
#> 
#> [[372]]$id
#> [1] "ENSG00000232810"
#> 
#> [[372]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[372]]$strand
#> [1] 1
#> 
#> [[372]]$version
#> [1] 3
#> 
#> [[372]]$gene_id
#> [1] "ENSG00000232810"
#> 
#> 
#> [[373]]
#> [[373]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[373]]$version
#> [1] 2
#> 
#> [[373]]$strand
#> [1] -1
#> 
#> [[373]]$gene_id
#> [1] "ENSG00000227507"
#> 
#> [[373]]$start
#> [1] 31548302
#> 
#> [[373]]$canonical_transcript
#> [1] "ENST00000429299.2"
#> 
#> [[373]]$external_name
#> [1] "LTB"
#> 
#> [[373]]$id
#> [1] "ENSG00000227507"
#> 
#> [[373]]$end
#> [1] 31550299
#> 
#> [[373]]$description
#> [1] "lymphotoxin beta (TNF superfamily, member 3) [Source:HGNC Symbol;Acc:6711]"
#> 
#> [[373]]$biotype
#> [1] "protein_coding"
#> 
#> [[373]]$feature_type
#> [1] "gene"
#> 
#> [[373]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[373]]$source
#> [1] "ensembl_havana"
#> 
#> [[373]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[374]]
#> [[374]]$source
#> [1] "ensembl_havana"
#> 
#> [[374]]$seq_region_name
#> [1] "6"
#> 
#> [[374]]$description
#> [1] "leukocyte specific transcript 1 [Source:HGNC Symbol;Acc:14189]"
#> 
#> [[374]]$biotype
#> [1] "protein_coding"
#> 
#> [[374]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[374]]$feature_type
#> [1] "gene"
#> 
#> [[374]]$start
#> [1] 31553901
#> 
#> [[374]]$external_name
#> [1] "LST1"
#> 
#> [[374]]$canonical_transcript
#> [1] "ENST00000376093.2"
#> 
#> [[374]]$id
#> [1] "ENSG00000204482"
#> 
#> [[374]]$end
#> [1] 31556686
#> 
#> [[374]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[374]]$version
#> [1] 6
#> 
#> [[374]]$strand
#> [1] 1
#> 
#> [[374]]$gene_id
#> [1] "ENSG00000204482"
#> 
#> 
#> [[375]]
#> [[375]]$description
#> [1] "natural cytotoxicity triggering receptor 3 [Source:HGNC Symbol;Acc:19077]"
#> 
#> [[375]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[375]]$feature_type
#> [1] "gene"
#> 
#> [[375]]$biotype
#> [1] "protein_coding"
#> 
#> [[375]]$seq_region_name
#> [1] "6"
#> 
#> [[375]]$source
#> [1] "ensembl_havana"
#> 
#> [[375]]$gene_id
#> [1] "ENSG00000204475"
#> 
#> [[375]]$strand
#> [1] -1
#> 
#> [[375]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[375]]$version
#> [1] 5
#> 
#> [[375]]$start
#> [1] 31556672
#> 
#> [[375]]$id
#> [1] "ENSG00000204475"
#> 
#> [[375]]$end
#> [1] 31560762
#> 
#> [[375]]$external_name
#> [1] "NCR3"
#> 
#> [[375]]$canonical_transcript
#> [1] "ENST00000340027.5"
#> 
#> 
#> [[376]]
#> [[376]]$strand
#> [1] 1
#> 
#> [[376]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[376]]$version
#> [1] 8
#> 
#> [[376]]$gene_id
#> [1] "ENSG00000204472"
#> 
#> [[376]]$start
#> [1] 31582961
#> 
#> [[376]]$external_name
#> [1] "AIF1"
#> 
#> [[376]]$canonical_transcript
#> [1] "ENST00000376059.3"
#> 
#> [[376]]$end
#> [1] 31584798
#> 
#> [[376]]$id
#> [1] "ENSG00000204472"
#> 
#> [[376]]$description
#> [1] "allograft inflammatory factor 1 [Source:HGNC Symbol;Acc:352]"
#> 
#> [[376]]$biotype
#> [1] "protein_coding"
#> 
#> [[376]]$feature_type
#> [1] "gene"
#> 
#> [[376]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[376]]$source
#> [1] "ensembl_havana"
#> 
#> [[376]]$seq_region_name
#> [1] "6"
#> 
#> 
#> [[377]]
#> [[377]]$source
#> [1] "ensembl_havana"
#> 
#> [[377]]$seq_region_name
#> [1] "6"
#> 
#> [[377]]$biotype
#> [1] "protein_coding"
#> 
#> [[377]]$feature_type
#> [1] "gene"
#> 
#> [[377]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[377]]$description
#> [1] "proline-rich coiled-coil 2A [Source:HGNC Symbol;Acc:13918]"
#> 
#> [[377]]$canonical_transcript
#> [1] "ENST00000376033.2"
#> 
#> [[377]]$external_name
#> [1] "PRRC2A"
#> 
#> [[377]]$id
#> [1] "ENSG00000204469"
#> 
#> [[377]]$end
#> [1] 31605548
#> 
#> [[377]]$start
#> [1] 31588497
#> 
#> [[377]]$version
#> [1] 8
#> 
#> [[377]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[377]]$strand
#> [1] 1
#> 
#> [[377]]$gene_id
#> [1] "ENSG00000204469"
#> 
#> 
#> [[378]]
#> [[378]]$start
#> [1] 31606805
#> 
#> [[378]]$canonical_transcript
#> [1] "ENST00000375964.6"
#> 
#> [[378]]$external_name
#> [1] "BAG6"
#> 
#> [[378]]$id
#> [1] "ENSG00000204463"
#> 
#> [[378]]$end
#> [1] 31620482
#> 
#> [[378]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[378]]$strand
#> [1] -1
#> 
#> [[378]]$version
#> [1] 8
#> 
#> [[378]]$gene_id
#> [1] "ENSG00000204463"
#> 
#> [[378]]$source
#> [1] "ensembl_havana"
#> 
#> [[378]]$seq_region_name
#> [1] "6"
#> 
#> [[378]]$description
#> [1] "BCL2-associated athanogene 6 [Source:HGNC Symbol;Acc:13919]"
#> 
#> [[378]]$biotype
#> [1] "protein_coding"
#> 
#> [[378]]$feature_type
#> [1] "gene"
#> 
#> [[378]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> 
#> [[379]]
#> [[379]]$description
#> [1] "apolipoprotein M [Source:HGNC Symbol;Acc:13916]"
#> 
#> [[379]]$biotype
#> [1] "protein_coding"
#> 
#> [[379]]$feature_type
#> [1] "gene"
#> 
#> [[379]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[379]]$source
#> [1] "ensembl_havana"
#> 
#> [[379]]$seq_region_name
#> [1] "6"
#> 
#> [[379]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[379]]$version
#> [1] 6
#> 
#> [[379]]$strand
#> [1] 1
#> 
#> [[379]]$gene_id
#> [1] "ENSG00000204444"
#> 
#> [[379]]$start
#> [1] 31620193
#> 
#> [[379]]$external_name
#> [1] "APOM"
#> 
#> [[379]]$canonical_transcript
#> [1] "ENST00000375916.3"
#> 
#> [[379]]$id
#> [1] "ENSG00000204444"
#> 
#> [[379]]$end
#> [1] 31625987
#> 
#> 
#> [[380]]
#> [[380]]$gene_id
#> [1] "ENSG00000204439"
#> 
#> [[380]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[380]]$strand
#> [1] -1
#> 
#> [[380]]$version
#> [1] 3
#> 
#> [[380]]$end
#> [1] 31628549
#> 
#> [[380]]$id
#> [1] "ENSG00000204439"
#> 
#> [[380]]$canonical_transcript
#> [1] "ENST00000375911.1"
#> 
#> [[380]]$external_name
#> [1] "C6orf47"
#> 
#> [[380]]$start
#> [1] 31626075
#> 
#> [[380]]$feature_type
#> [1] "gene"
#> 
#> [[380]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[380]]$biotype
#> [1] "protein_coding"
#> 
#> [[380]]$description
#> [1] "chromosome 6 open reading frame 47 [Source:HGNC Symbol;Acc:19076]"
#> 
#> [[380]]$seq_region_name
#> [1] "6"
#> 
#> [[380]]$source
#> [1] "ensembl_havana"
#> 
#> 
#> [[381]]
#> [[381]]$external_name
#> [1] "GPANK1"
#> 
#> [[381]]$canonical_transcript
#> [1] "ENST00000375906.1"
#> 
#> [[381]]$end
#> [1] 31634060
#> 
#> [[381]]$id
#> [1] "ENSG00000204438"
#> 
#> [[381]]$start
#> [1] 31629006
#> 
#> [[381]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[381]]$version
#> [1] 6
#> 
#> [[381]]$strand
#> [1] -1
#> 
#> [[381]]$gene_id
#> [1] "ENSG00000204438"
#> 
#> [[381]]$source
#> [1] "ensembl_havana"
#> 
#> [[381]]$seq_region_name
#> [1] "6"
#> 
#> [[381]]$biotype
#> [1] "protein_coding"
#> 
#> [[381]]$feature_type
#> [1] "gene"
#> 
#> [[381]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[381]]$description
#> [1] "G patch domain and ankyrin repeats 1 [Source:HGNC Symbol;Acc:13920]"
#> 
#> 
#> [[382]]
#> [[382]]$gene_id
#> [1] "ENSG00000204435"
#> 
#> [[382]]$assembly_name
#> [1] "GRCh37"
#> 
#> [[382]]$strand
#> [1] 1
#> 
#> [[382]]$version
#> [1] 9
#> 
#> [[382]]$start
#> [1] 31633013
#> 
#> [[382]]$end
#> [1] 31638120
#> 
#> [[382]]$id
#> [1] "ENSG00000204435"
#> 
#> [[382]]$canonical_transcript
#> [1] "ENST00000375882.2"
#> 
#> [[382]]$external_name
#> [1] "CSNK2B"
#> 
#> [[382]]$description
#> [1] "casein kinase 2, beta polypeptide [Source:HGNC Symbol;Acc:2460]"
#> 
#> [[382]]$feature_type
#> [1] "gene"
#> 
#> [[382]]$logic_name
#> [1] "ensembl_havana_gene_homo_sapiens_37"
#> 
#> [[382]]$biotype
#> [1] "protein_coding"
#> 
#> [[382]]$seq_region_name
#> [1] "6"
#> 
#> [[382]]$source
#> [1] "ensembl_havana"
#> MHC genes detected: 382
#> MHC genes with unique HGNC ID detected: 374
get_gene_coordinates(mhc_genes %>% pull(external_name))
#> 
#> 
#> Dowloading MAGMA hg19 reference file ...
#> # A tibble: 19,261 × 5
#>    chr    start    end entrez    hgnc        
#>    <chr>  <int>  <int> <chr>     <chr>       
#>  1 chr1   69091  70008 79501     OR4F5       
#>  2 chr1  142447 174392 100996442 LOC100996442
#>  3 chr1  367659 368597 729759    OR4F29      
#>  4 chr1  621096 622034 81399     OR4F16      
#>  5 chr1  859993 879961 148398    SAMD11      
#>  6 chr1  879583 894679 26155     NOC2L       
#>  7 chr1  895967 901099 339451    KLHL17      
#>  8 chr1  901872 910488 84069     PLEKHN1     
#>  9 chr1  910579 917473 84808     PERM1       
#> 10 chr1  934342 936608 57801     HES4        
#> # ℹ 19,251 more rows
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
