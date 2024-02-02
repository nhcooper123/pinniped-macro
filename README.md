# pinniped-macro
Code and data for pinniped macroevolution paper

Author(s): Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper

*This README is a work in progess*

This repository contains all the code and some data used in the [paper](XXX). 

To cite the paper: 
>  Park T, Burin G, Lazo-Cancino D, Rees J.P.G., Rule J.P, Slater G.J, Cooper N. Charting the Course of Pinniped Evolution: insights from molecular phylogeny and fossil record integration. Evolution.

To cite this repo: 
>  Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper. Code for the paper v1.

[![DOI](https://zenodo.org/badge/656148711.svg)](https://zenodo.org/doi/10.5281/zenodo.10611047)

![alt text](https://github.com/nhcooper123/pinniped-macro/raw/main/pinnped_tree_2.jpg)

This paper has three main components, a metatree of Pinnipedia, biogeographical analyses and diversification analyses.

-------
## Metatree

Code for the metatree is within the `metatree-scripts/` folder. Note that the data etc. are too large for GitHub so are available from [Zenodo](https://doi.org/10.5281/zenodo.8276115). The file structure of that repository is included below for completeness.

If you use the data please cite as follows: 
>  Travis Park and Graham Slater (2023). Pinniped metatree data and scripts (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.8276115

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8276115.svg)](https://doi.org/10.5281/zenodo.8276115)

The set of R scripts used to generate the metatree are:

1. **BuildBeastXML.R** A function used to assemble a nexus file, starting tree, topology constraints, and age priors for BEAST input. Is sourced by the file CreateBeastInput_Script.R
2. **Build_Metatree_Files_Combined.R** Script to read in input files from input_files folder and assemble MRP matrices for generation of metatree
2. **CreateBeastInput_Script.R** Script to read in a set of input files contained in the TimeTree folder and output a set of files that can be used to build a BEAST XML for Fossilized Birth Death tip dating analyses with a relaxed molecular clock
2. **makemrpandscandmpts.R** Utility script to read the output of TNT for each source study, produce separate MPT and Strict Cosensus files, and generate an MRP file that can be used for input into the metatree function
2. **process_beast_output.R** Script to read combined, post-burnin output from beast and find the median tree under the Kendall-Colijn distance metric. This tree is used in subsequent macroevolutionary analyses.
2. **runtntoutputcollator.R** Utility script to read the output from TNT analysis of the ensemble MRP matrix, find all unique shortest trees, and assemble a combined MPT file, consensus trees, and perform safe taxonomic reinsertion. 

### `Pinniped_metatree` folder ([Zenodo](https://doi.org/10.5281/zenodo.8276115))

The file structure of that repository is as follows: 

### InputData
Folder containing the following subfolders:

1. **Molecular_Data**. Folder containing two subdirectories (fulton and Lopes), each of which contains a further two subdirectories (mrp and xml) that contain the Matrix Representations with Parsimony (MRP) nexus files and XML files (see below), respectivly, for the two molecular studies used as backbone constraints
2. **MRP**. Matrix Representations with Parsimony (MRPs) obtained by re-analysing each morphological input data set under parsimony and retaining only unique biparitions until all such biparitions have been sampled.
3. **NEXUS**. The original character-taxon matrices for each morphological input data set.
4. **TNT**. TNT format input files used for reanalysis of each morphological input dataset. Each file contains the relevant commands to run the analyses.
5. **XML**. XML files that record important metadata about each input data set, particularly data set dependence and taxonomic reconciliation.

### metatree_files
Folder containing the metatree data organised in the following subfolders:

1. **Consensus_trees**. Folder containing Strict (SCC) and Majority Rule (MRC) trees before (STR_***) and after (no prefix) safely reinserting taxa
2. **Files**. All files output by running the metatree function to assemble a reconciled complete MRP matrix over all input datasets:
    - **Character weights**. Weights for each character
    - **DataSetWeights.txt**. Information of the contribution of each dataset to the overall analysis 
    - **FULL.nex**. MRP files in nexus format containing all taxa (not used in main analysis)
    - **FULL.tnt**. MRP files in tnt format containing all taxa (not used in main analysis)
    - **STR.nex**. Safely reduced input files in nexus format.
    - **STR.tnt**. Safely reduced input files in tnt format. The tnt format file is the one used in main analyses.
3. **STR.tnt**. Text file containing safely removed taxa and the rules used to reinsert them
    - **TaxonomyTree.tre** Newick string containing the paleobiology database taxonomy for PanPinnipeds in tree form. Used in the metatree analysis as a very downweighted input.
5. **MPTS** Folder containing all most parsimonious trees before and after Safe Taxonomic Reinsertion
6. **tnt** MPTs from each of 1000 independent TNT searches (.tnt) and associated search info from screen output (.txt) 

### Scripts

See `metatree-scripts/` description above.

### TimeTree
Files and Folders related to BEAST time-tree inference. Organized as follows:

1. **BEAST_Run**. Folder containing input xml file, log and tree files from the two independent MCMC analyses, combined post-burnin log and tree files, and the median tree from the combined posterior sample.
1. **partitionFinder**. PartitionFinder input and output files used to identify the appropriate partitioning scheme and evolutionary models for the augmented Fulton and Strobeck alignment. 
1. **strictclockru**. xml and output files from a short analysis of extant pinnipeds only under a strict molecular clock. The posterior mean root height divided by the approximate age of crown pinnipeds was used to define a mean for the lognormal prior on the relaxed molecular clock in the main analysis
1. **MRC.tre**. The majority Rule Consensus metatree used to derive topological constraints for the BEAST analysis
1. **pinnipeds_extant_beast.nex**. Augmented version of the Fulton and Strobeck alignment used in BEAST analysis
1. **Taxa_FAD.csv**. Stratigraphic uncertainty associated with the first appearances of fossil taxa. Used to define age ranges for fossils for the Fossilized Birth Death analysis.

### TimeTreeInference
Folder containing the R scripts and temporal data used to build the timetrees as well as the following subfolders:

1. **Dangerous**. Subfolder containing the results for the Dangerous timetree analyses.
1. **Risky**. Subfolder containing the results for the Risky timetree analyses.
1. **Safe**. Subfolder containing the results for the Safe timetree analyses.


-------
## Biogeography analyses

Code for the biogeographical analyses is within the `biogeography/` folder

### Data
All data are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq), but also within the  `data/`  and `raw-data/` folders.

If you use the data please cite as follows: 
>  Travis Park; Gustavo Burin; Graham J Slater; Natalie Cooper (2022). Data from the "Back to the water" project [Data set]. Natural History Museum. https://doi.org/10.5519/vmbrpkuq.

### Analyses

R-scripts 01-05 are modified versions of BioGeoBEARS R script(s) that can be found at [http://phylo.wikidot.com/biogeobears](http://phylo.wikidot.com/biogeobears). Copyright Nicholas J. Matzke. Please cite BioGeoBEARS if you use these scripts.

- **00-prepare-trees** removes species without geographical information from the trees.
- **01-pinniped-state-list-fix.R** modifies the states allowed for each analysis. THIS NEEDS TO BE RUN BEFORE SCRIPT 02. It is sourced in most other scripts to prevent you forgetting.
- **02A-fit-DEC-DECj** fits DEC and DECj models using BioGeoBEARS for all taxa removing impossible states.
- **02B-fit-DEC-DECj-unlikely** fits DEC and DECj models using BioGeoBEARS for all taxa removing impossible and unlikely states.
- **03A-fit-DEC-DECj_fossil** fits DEC and DECj models using BioGeoBEARS for fossil taxa removing impossible states.
- **03B-fit-DEC-DECj-unlikely_fossil** fits DEC and DECj models using BioGeoBEARS for fossil taxa removing impossible and unlikely states.
- **04A-fit-DEC-DECj_extant** fits DEC and DECj models using BioGeoBEARS for extant taxa removing impossible states.
- **04B-fit-DEC-DECj-unlikely_extant** fits DEC and DECj models using BioGeoBEARS for extant taxa removing impossible and unlikely states.
- **05-DEC-DECj-results.R** extracts results from DEC and DECj models for tables. 
- **06-DEC-DECj-plots.R** creates all BioGeoBEARS plots for the ESM.
- **07-Extract-BGB-results-for-plotting.R** extracts results from BioGeoBEARS results objects for making the main publication plot.
- **08-main-results-plot.R** code for the nice main publication plot.
- **09-worldmap-for-legend.R** quick snippet to get a world map for the legend.
- **messing-about-with-colours.R** attempts to set nice colour palettes for scripts 06, 07 and 08 which call it.
- **modifiedBGBplot.R** modified version of BioGeoBEARS plotting function to have more control over plots in 06.

The folder `biogeography/outputs` contains the main figures and tables. ESM figures/tables are in the higher level `/supplemental` folder.

-------

## Diversification rate analyses

Code for the diversification rate analyses is within the `diversification/` folder

### Data
All data are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq), but also within the  `data/` folder. These consist of the following:

1. **pinniped_median_full.tre** - full phylogeny.
1. **pinniped_median_extant.tre** - phylogeny with only extant taxa.  
1. **pinniped_median_fossil.tre** - phylogeny with only fossil taxa.
1. **pinniped_median_noanc.tre** - full phylogeny but with sampled ancestors removed.
1. **pinniped_median_fossil_noanc.tre** - phylogeny with only fossil taxa but with sampled ancestors removed.

If you use the data please cite as follows: 
>  Travis Park; Gustavo Burin; Graham J Slater; Natalie Cooper (2022). Data from the "Back to the water" project [Data set]. Natural History Museum. 

### FossilBAMM analyses
The fossilBAMM analyses are in the `analyses/` folder organised in folders according to which subset of the data were used (i.e. all taxa, extant only, fossil only, with or without sampled ancestors). Within each folder there is a control file for fossilBAMM (in `main_analysis/` this is called `divcontrol_pinnipedia.txt`) and the output files from fossilBAMM for use in further analyses. For `main_analysis/` these are:

1. `chain_swap_pinnipedia_noanc.txt`
1. `main_noanc.err`
2. `mcmc_out_pinnipedia_noanc.txt`
3. `run_info_pinnipedia_noanc.txt`
4. `event_data_pinnipedia_noanc.txt`

This folder also contains sensitivity analyses (withing the folders `sensitivity_analyses` and `sensitivity_analyses_with_sampled_ancestors`) where the analyses were repeated but using priors of 0.01, 0.1, 0.5, 2, or 10. Within the sensitivity analyses folders there are subfolders for each different prior and then within these the same files as described above for the main analyses.

#### BAMMTools analyses in R

These scripts process the outputs in the `analyses/` folder using BAMMTools in R to extract speciation, extinction and diversification rates, credicle rate shifts and plotting etc.

1. **fossilbamm_results_full.R** - code for all taxa with sampled ancestors.
2. **fossilbamm_results_noanc.R** - code for all taxa without sampled ancestors.
3. **fossilbamm_results_extant.R** - code for extant taxa.
1. **fossilbamm_results_fossil_only.R** - code for fossil taxa.
3. **fossilbamm_results_fossil_only_with_sampled_ancestors.R** - code for fossil taxa with sampled ancestors.
1. **fossilbamm_results_full_sensitivity.R** - code for sensitivity analyses on all taxa with sampled ancestors.
1. **fossilbamm_results_noanc_sensitivity.R** - code for sensitivity analyses on all taxa without sampled ancestors.

-------
## Other folders

* `/supplemental` contains the ESM files/figures and tables in LaTeX format.

-------
## Session Info
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication.

    ─ Session info ─────────────────────────────────────────────────────────────
    setting  value
    version  R version 4.2.2 (2022-10-31)
    os       macOS Monterey 12.5.1
    system   x86_64, darwin17.0
     ui       RStudio
     language (EN)
     collate  en_US.UTF-8
     ctype    en_US.UTF-8
     tz       Europe/London
     date     2024-02-02
     rstudio  2023.12.0+369 Ocean Storm (desktop)
     pandoc   3.1.1 @ /usr/local/bin/pandoc

    ─ Packages ─────────────────────────────────────────────────────────────────
     package           * version    date (UTC) lib source
     abind               1.4-5      2016-07-21 [1] CRAN (R 4.2.0)
     ade4                1.7-22     2023-02-06 [1] CRAN (R 4.2.0)
     ape               * 5.7        2023-02-16 [1] CRAN (R 4.2.0)
     aplot               0.1.9      2022-11-24 [1] CRAN (R 4.2.0)
     assertthat          0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
     backports           1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
     BAMMtools         * 2.1.10     2022-07-15 [1] CRAN (R 4.2.0)
     BioGeoBEARS       * 1.1.2      2023-02-16 [1] Github (nmatzke/BioGeoBEARS@3af1a12)
     bitops              1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
     broom               1.0.3      2023-01-25 [1] CRAN (R 4.2.2)
     cachem              1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
     callr               3.7.3      2022-11-02 [1] CRAN (R 4.2.0)
     caTools             1.18.2     2021-03-28 [1] CRAN (R 4.2.0)
     cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
     cladoRcpp         * 0.15.1     2023-02-16 [1] Github (nmatzke/cladoRcpp@3cdd71e)
     cli                 3.6.0      2023-01-09 [1] CRAN (R 4.2.0)
     cluster             2.1.4      2022-08-22 [2] CRAN (R 4.2.2)
     clusterGeneration   1.3.7      2020-12-15 [1] CRAN (R 4.2.0)
     coda                0.19-4     2020-09-30 [1] CRAN (R 4.2.0)
     codetools           0.2-19     2023-02-01 [2] CRAN (R 4.2.0)
     colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.2.0)
     combinat            0.0-8      2012-10-29 [1] CRAN (R 4.2.0)
     crayon              1.5.2      2022-09-29 [1] CRAN (R 4.2.0)
     DBI                 1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
     dbplyr              2.3.0      2023-01-16 [1] CRAN (R 4.2.0)
     deSolve             1.34       2022-10-22 [1] CRAN (R 4.2.0)
     devtools          * 2.4.5      2022-10-11 [1] CRAN (R 4.2.0)
     dichromat           2.0-0.1    2022-05-02 [1] CRAN (R 4.2.0)
     digest              0.6.31     2022-12-11 [1] CRAN (R 4.2.0)
     doParallel          1.0.17     2022-02-07 [1] CRAN (R 4.2.0)
     dotCall64           1.0-2      2022-10-03 [1] CRAN (R 4.2.0)
     dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.2.0)
     ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
     expm                0.999-7    2023-01-09 [1] CRAN (R 4.2.0)
     fansi               1.0.4      2023-01-22 [1] CRAN (R 4.2.0)
     fastmap             1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
     fastmatch           1.1-3      2021-07-23 [1] CRAN (R 4.2.0)
     FD                  1.0-12.1   2022-05-02 [1] CRAN (R 4.2.0)
     fdrtool             1.2.17     2021-11-13 [1] CRAN (R 4.2.0)
     forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.2.0)
     foreach             1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
     fs                  1.6.1      2023-02-06 [1] CRAN (R 4.2.0)
     gargle              1.3.0      2023-01-30 [1] CRAN (R 4.2.0)
     gdata               2.18.0.1   2022-05-10 [1] CRAN (R 4.2.0)
     geiger            * 2.0.10     2022-06-03 [1] CRAN (R 4.2.0)
     generics            0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
     GenSA               1.1.7      2018-01-17 [1] CRAN (R 4.2.0)
     geometry            0.4.7      2023-02-03 [1] CRAN (R 4.2.0)
     ggfun               0.0.9      2022-11-21 [1] CRAN (R 4.2.0)
     ggimage           * 0.3.1      2022-04-25 [1] CRAN (R 4.2.0)
     ggplot2           * 3.4.4      2023-10-12 [1] CRAN (R 4.2.0)
     ggplotify           0.1.0      2021-09-02 [1] CRAN (R 4.2.0)
     ggtree            * 3.6.2      2022-11-10 [1] Bioconductor
     glue                1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
     googledrive         2.0.0      2021-07-08 [1] CRAN (R 4.2.0)
     googlesheets4       1.0.1      2022-08-13 [1] CRAN (R 4.2.0)
     gplots              3.1.3      2022-04-25 [1] CRAN (R 4.2.0)
     gridGraphics        0.5-1      2020-12-13 [1] CRAN (R 4.2.0)
     gtable              0.3.1      2022-09-01 [1] CRAN (R 4.2.0)
     gtools              3.9.4      2022-11-27 [1] CRAN (R 4.2.0)
     haven               2.5.1      2022-08-22 [1] CRAN (R 4.2.0)
     here                1.0.1      2020-12-13 [1] CRAN (R 4.2.0)
     hms                 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
     htmltools           0.5.4      2022-12-07 [1] CRAN (R 4.2.0)
     htmlwidgets         1.6.1      2023-01-07 [1] CRAN (R 4.2.0)
     httpuv              1.6.9      2023-02-14 [1] CRAN (R 4.2.2)
     httr                1.4.4      2022-08-17 [1] CRAN (R 4.2.0)
     igraph              1.4.1      2023-02-24 [1] CRAN (R 4.2.0)
     iterators           1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
     jsonlite            1.8.4      2022-12-06 [1] CRAN (R 4.2.0)
     KernSmooth          2.23-20    2021-05-03 [2] CRAN (R 4.2.2)
     later               1.3.0      2021-08-18 [1] CRAN (R 4.2.0)
     lattice             0.20-45    2021-09-22 [2] CRAN (R 4.2.2)
     lazyeval            0.2.2      2019-03-15 [1] CRAN (R 4.2.0)
     lifecycle           1.0.3      2022-10-07 [1] CRAN (R 4.2.0)
     lubridate           1.9.1      2023-01-24 [1] CRAN (R 4.2.2)
     magic               1.6-1      2022-11-16 [1] CRAN (R 4.2.0)
     magick              2.7.3      2021-08-18 [1] CRAN (R 4.2.0)
     magrittr            2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
     mapproj             1.2.11     2023-01-12 [1] CRAN (R 4.2.0)
     maps              * 3.4.1      2022-10-30 [1] CRAN (R 4.2.0)
     MASS                7.3-58.2   2023-01-23 [2] CRAN (R 4.2.0)
     Matrix              1.5-3      2022-11-11 [2] CRAN (R 4.2.0)
     memoise             2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
     mgcv                1.8-41     2022-10-21 [2] CRAN (R 4.2.2)
     mime                0.12       2021-09-28 [1] CRAN (R 4.2.0)
     miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 4.2.0)
     mnormt              2.1.1      2022-09-26 [1] CRAN (R 4.2.0)
     modelr              0.1.10     2022-11-11 [1] CRAN (R 4.2.0)
     MultinomialCI       1.2        2021-05-11 [1] CRAN (R 4.2.0)
     munsell             0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
     mvtnorm             1.1-3      2021-10-08 [1] CRAN (R 4.2.0)
     nlme                3.1-162    2023-01-31 [2] CRAN (R 4.2.0)
     numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
     optimParallel       1.0-2      2021-02-11 [1] CRAN (R 4.2.0)
     optimx              2022-4.30  2022-05-10 [1] CRAN (R 4.2.0)
     pals              * 1.8        2023-08-23 [1] CRAN (R 4.2.0)
     patchwork         * 1.2.0      2024-01-08 [1] CRAN (R 4.2.2)
     permute             0.9-7      2022-01-27 [1] CRAN (R 4.2.0)
     phangorn            2.11.1     2023-01-23 [1] CRAN (R 4.2.0)
     phylobase           0.8.10     2020-03-01 [1] CRAN (R 4.2.0)
     phytools          * 1.5-1      2023-02-19 [1] CRAN (R 4.2.0)
     pillar              1.9.0      2023-03-22 [1] CRAN (R 4.2.0)
     pkgbuild            1.4.0      2022-11-27 [1] CRAN (R 4.2.0)
     pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
     pkgload             1.3.2      2022-11-16 [1] CRAN (R 4.2.0)
     plotrix             3.8-2      2021-09-08 [1] CRAN (R 4.2.0)
     plyr                1.8.8      2022-11-11 [1] CRAN (R 4.2.0)
     prettyunits         1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
     processx            3.8.0      2022-10-26 [1] CRAN (R 4.2.0)
     profvis             0.3.7      2020-11-02 [1] CRAN (R 4.2.0)
     progress            1.2.2      2019-05-16 [1] CRAN (R 4.2.0)
     promises            1.2.0.1    2021-02-11 [1] CRAN (R 4.2.0)
     ps                  1.7.2      2022-10-26 [1] CRAN (R 4.2.0)
     purrr             * 1.0.1      2023-01-10 [1] CRAN (R 4.2.2)
     quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.2.0)
     R6                  2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
     Rcpp                1.0.10     2023-01-22 [1] CRAN (R 4.2.0)
     readr             * 2.1.3      2022-10-01 [1] CRAN (R 4.2.0)
     readxl              1.4.1      2022-08-17 [1] CRAN (R 4.2.0)
     remotes             2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
     reprex              2.0.2      2022-08-17 [1] CRAN (R 4.2.0)
     reshape2            1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
     rexpokit            0.26.6.7   2020-07-04 [1] CRAN (R 4.2.0)
     rlang               1.1.3      2024-01-10 [1] CRAN (R 4.2.2)
     rncl                0.8.7      2023-01-08 [1] CRAN (R 4.2.0)
     RNeXML              2.4.11     2023-02-01 [1] CRAN (R 4.2.0)
     rprojroot           2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
     rstudioapi          0.14       2022-08-22 [1] CRAN (R 4.2.0)
     rvest               1.0.3      2022-08-19 [1] CRAN (R 4.2.0)
     scales              1.2.1      2022-08-20 [1] CRAN (R 4.2.0)
     scatterplot3d       0.3-42     2022-09-08 [1] CRAN (R 4.2.0)
     sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
     shiny               1.7.4      2022-12-15 [1] CRAN (R 4.2.0)
     spam                2.9-1      2022-08-07 [1] CRAN (R 4.2.0)
     SparseM             1.81       2021-02-18 [1] CRAN (R 4.2.0)
     statmod             1.5.0      2023-01-06 [1] CRAN (R 4.2.0)
     stringi             1.7.12     2023-01-11 [1] CRAN (R 4.2.2)
     stringr           * 1.5.0      2022-12-02 [1] CRAN (R 4.2.0)
     subplex             1.8        2022-04-12 [1] CRAN (R 4.2.0)
     tibble            * 3.2.1      2023-03-20 [1] CRAN (R 4.2.0)
     tidyr             * 1.3.0      2023-01-24 [1] CRAN (R 4.2.2)
     tidyselect          1.2.0      2022-10-10 [1] CRAN (R 4.2.0)
     tidytree            0.4.2      2022-12-18 [1] CRAN (R 4.2.0)
     tidyverse         * 1.3.2      2022-07-18 [1] CRAN (R 4.2.0)
     timechange          0.2.0      2023-01-11 [1] CRAN (R 4.2.2)
     treeio              1.22.0     2022-11-01 [1] Bioconductor
     tzdb                0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
     urlchecker          1.0.1      2021-11-30 [1] CRAN (R 4.2.0)
     usethis           * 2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
     utf8                1.2.3      2023-01-31 [1] CRAN (R 4.2.0)
     uuid                1.1-0      2022-04-19 [1] CRAN (R 4.2.0)
     vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.2.0)
     vegan               2.6-4      2022-10-11 [1] CRAN (R 4.2.0)
     withr               2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
     XML                 3.99-0.13  2022-12-04 [1] CRAN (R 4.2.0)
     xml2                1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
     xtable              1.8-4      2019-04-21 [1] CRAN (R 4.2.0)
     yulab.utils         0.0.6      2022-12-20 [1] CRAN (R 4.2.0)
   
## Checkpoint for reproducibility
To rerun all the code with packages as they existed on CRAN at time of our analyses we recommend using the `checkpoint` package, and running this code prior to the analysis:

```{r}
checkpoint("2024-02-02")
```
