# MsRDB-Manuscript-Code
This is the repository archiving code and data for MsRDB manuscript -- Multi-scale Adaptive Differential Abundance Analysis in
Microbial Compositional Data

## File Introduction

The "MsRDB" directory contains files for algorithms and results. Algorithm.R includes the code of new MsRDB test, ASVwise.R includes the code of RDb test and KNN-RDB test, and MultiATE.R includes the code of a generalized version of ATE function in ATE package. Simulation.R includes code of numerical experiments and visualization in Section 2.2, Vangay_Cell_2018.R includes code of numerical experiments and visualization in Section 2.3, and Bokulich_PNAS_2013.R includes code of numerical experiments and visualization in Section 2.4. <br />
<br />
The "Data" directory contains all data set used in the manuscript. <br />
<br />

## Results Replication

First, run all code in Algorithm.R, ASVwise.R and MultiATE.R. Then, the results can be replicated by running one of Simulation.R, Vangay_Cell_2018.R, and Bokulich_PNAS_2013.R.  <br />

## Session Information

R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 12.0.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] doRNG_1.8.2         rngtools_1.5        doParallel_1.0.16   iterators_1.0.13    foreach_1.5.1       microbiome_1.14.0   ALDEx2_1.24.0      
 [8] Rfast_2.0.3         RcppZiggurat_0.1.6  zCompositions_1.3.4 truncnorm_1.0-8     NADA_1.6-1.1        survival_3.2-12     MASS_7.3-54        
[15] ANCOMBC_1.2.2       phyloseq_1.36.0     forcats_0.5.1       stringr_1.4.0       dplyr_1.0.7         purrr_0.3.4         readr_2.0.1        
[22] tidyr_1.1.3         tibble_3.1.3        ggplot2_3.3.5       tidyverse_1.3.1     biomformat_1.20.0   igraph_1.2.6        DECIPHER_2.20.0    
[29] RSQLite_2.2.8       Biostrings_2.60.2   GenomeInfoDb_1.28.1 XVector_0.32.0      IRanges_2.26.0      S4Vectors_0.30.0    BiocGenerics_0.38.0
[36] kmer_1.1.2          dada2_1.20.0        Rcpp_1.0.7         

loaded via a namespace (and not attached):
 [1] readxl_1.3.1                backports_1.2.1             plyr_1.8.6                  splines_4.1.0               BiocParallel_1.26.2        
 [6] digest_0.6.27               fansi_0.5.0                 magrittr_2.0.1              memoise_2.0.0               cluster_2.1.2              
[11] tzdb_0.1.2                  modelr_0.1.8                RcppParallel_5.1.4          matrixStats_0.60.1          jpeg_0.1-9                 
[16] colorspace_2.0-2            blob_1.2.2                  rvest_1.0.1                 haven_2.4.3                 rbibutils_2.2.3            
[21] crayon_1.4.1                RCurl_1.98-1.4              jsonlite_1.7.2              ape_5.5                     glue_1.4.2                 
[26] gtable_0.3.0                phylogram_2.1.0             zlibbioc_1.38.0             DelayedArray_0.18.0         Rhdf5lib_1.14.2            
[31] scales_1.1.1                futile.options_1.0.1        DBI_1.1.1                   bit_4.0.4                   httr_1.4.2                 
[36] RColorBrewer_1.1-2          ellipsis_0.3.2              pkgconfig_2.0.3             farver_2.1.0                dbplyr_2.1.1               
[41] utf8_1.2.2                  tidyselect_1.1.1            labeling_0.4.2              rlang_0.4.11                reshape2_1.4.4             
[46] munsell_0.5.0               cellranger_1.1.0            tools_4.1.0                 cachem_1.0.6                cli_3.0.1                  
[51] generics_0.1.0              ade4_1.7-17                 broom_0.7.9                 fastmap_1.1.0               bit64_4.0.5                
[56] fs_1.5.0                    nlme_3.1-152                formatR_1.11                xml2_1.3.2                  compiler_4.1.0             
[61] rstudioapi_0.13             png_0.1-7                   reprex_2.0.1                stringi_1.7.3               futile.logger_1.4.3        
[66] lattice_0.20-44             Matrix_1.3-4                nloptr_1.2.2.2              vegan_2.5-7                 permute_0.9-5              
[71] multtest_2.48.0             vctrs_0.3.8                 pillar_1.6.2                lifecycle_1.0.0             rhdf5filters_1.4.0         
[76] Rdpack_2.1.2                data.table_1.14.0           bitops_1.0-7                GenomicRanges_1.44.0        R6_2.5.1                   
[81] latticeExtra_0.6-29         hwriter_1.3.2               ShortRead_1.50.0            codetools_0.2-18            lambda.r_1.2.4             
[86] assertthat_0.2.1            rhdf5_2.36.0                SummarizedExperiment_1.22.0 withr_2.4.2                 GenomicAlignments_1.28.0   
[91] Rsamtools_2.8.0             GenomeInfoDbData_1.2.6      mgcv_1.8-36                 hms_1.1.0                   grid_4.1.0                 
[96] MatrixGenerics_1.4.2        Rtsne_0.15                  Biobase_2.52.0              lubridate_1.7.10 