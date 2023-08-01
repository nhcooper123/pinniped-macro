# pinniped-macro
Code and data for pinniped macroevolution paper

Author(s): Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper

*This README is a work in progess*

This repository contains all the code and some data used in the [paper](XXX). 

To cite the paper: 
>  Park T, Burin G, Lazo-Cancino D, Rees J.P.G., Rule J.P, Slater G.J, Cooper N. Macroevolutionary patterns of pinnipeds. in prep.

To cite this repo: 
>  Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper. Code for the paper v1.

![alt text](https://github.com/nhcooper123/pinniped-macro/raw/pinnped_tree_2.jpg)

This paper has three main components, a metatree of Pinnipedia, biogeographical analyses and diversification analyses.

-------
## Metatree

Code for the metatree is within the `metatree-scripts/` folder. Note that the data etc. are too large for GitHub so are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq) within the `Pinniped_metatree` folder. The file structure of that repository is included below for completeness.

The set of R scripts used to generate the metatree are:

1. **BuildBeastXML.R** A function used to assemble a nexus file, starting tree, topology constraints, and age priors for BEAST input. Is sourced by the file CreateBeastInput_Script.R
2. **Build_Metatree_Files_Combined.R** Script to read in input files from input_files folder and assemble MRP matrices for generation of metatree
2. **CreateBeastInput_Script.R** Script to read in a set of input files contained in the TimeTree folder and output a set of files that can be used to build a BEAST XML for Fossilized Birth Death tip dating analyses with a relaxed molecular clock
2. **makemrpandscandmpts.R** Utility script to read the output of TNT for each source study, produce separate MPT and Strict Cosensus files, and generate an MRP file that can be used for input into the metatree function
2. **process_beast_output.R** Script to read combined, post-burnin output from beast and find the median tree under the Kendall-Colijn distance metric. This tree is used in subsequent macroevolutionary analyses.
2. **runtntoutputcollator.R** Utility script to read the output from TNT analysis of the ensemble MRP matrix, find all unique shortest trees, and assemble a combined MPT file, consensus trees, and perform safe taxonomic reinsertion. 

### `Pinniped_metatree` folder ([NHM Data Portal](https://doi.org/10.5519/vmbrpkuq))

The file structure of that repository is as follows: 

### InputData
Folder containing the following subfolders:

1. **Molecular_Data** Folder containing two subdirectories (fulton and Lopes), each of which contains a further two subdirectories (mrp and xml) that contain the Matrix Representations with Parsimony (MRP) nexus files and XML files (see below), respectivly, for the two molecular studies used as backbone constraints
2. **MRP** Matrix Representations with Parsimony (MRPs) obtained by re-analysing each morphological input data set under parsimony and retaining only unique biparitions until all such biparitions have been sampled.
3. **NEXUS** The original character-taxon matrices for each morphological input data set.
4. **TNT** TNT format input files used for reanalysis of each morphological input dataset. Each file contains the relevant commands to run the analyses.
5. **XML** XML files that record important metadata about each input data set, particularly data set dependence and taxonomic reconciliation.

### metatree_files
Folder containing the metatree data organised in the following subfolders:

1. **Consensus_trees** Folder containing Strict (SCC) and Majority Rule (MRC) trees before (STR_***) and after (no prefix) safely reinserting taxa
2. **Files** All files output by running the metatree function to assemble a reconciled complete MRP matrix over all input datasets:
    - **Character weights** weights for each character
    - **DataSetWeights.txt** information of the contribution of each dataset to the overall analysis 
    - **FULL.nex** MRP files in nexus format containing all taxa (not used in main analysis)
    - **FULL.tnt** MRP files in tnt format containing all taxa (not used in main analysis)
    - **STR.nex** Safely reduced input files in nexus format.
    - **STR.tnt** Safely reduced input files in tnt format. The tnt format file is the one used in main analyses.
3. **STR.tnt** Text file containing safely removed taxa and the rules used to reinsert them
    - **TaxonomyTree.tre** Newick string containing the paleobiology database taxonomy for PanPinnipeds in tree form. Used in the metatree analysis as a very downweighted input.
5. **MPTS** Folder containing all most parsimonious trees before and after Safe Taxonomic Reinsertion
6. **tnt** MPTs from each of 1000 independent TNT searches (.tnt) and associated search info from screen output (.txt) 

### Scripts

See above

### TimeTree
Files and Folders related to BEAST time-tree inference. Organized as follows:

1. **BEAST_Run** Folder containing input xml file, log and tree files from the two independent MCMC analyses, combined post-burnin log and tree files, and the median tree from the combined posterior sample.
1. **partitionFinder** PartitionFinder input and output files used to identify the appropriate partitioning scheme and evolutionary models for the augmented Fulton and Strobeck alignment. 
1. **strictclockru** xml and output files from a short analysis of extant pinnipeds only under a strict molecular clock. The posterior mean root height divided by the approximate age of crown pinnipeds was used to define a mean for the lognormal prior on the relaxed molecular clock in the main analysis
1. **MRC.tre** The majority Rule Consensus metatree used to derive topological constraints for the BEAST analysis
1. **pinnipeds_extant_beast.nex** Augmented version of the Fulton and Strobeck alignment used in BEAST analysis
1. **Taxa_FAD.csv** Stratigraphic uncertainty associated with the first appearances of fossil taxa. Used to define age ranges for fossils for the Fossilized Birth Death analysis.

### TimeTreeInference
Folder containing the R scripts and temporal data used to build the timetrees as well as the following subfolders:

1. **Dangerous** Subfolder containing the results for the Dangerous timetree analyses.
1. **Risky** Subfolder containing the results for the Risky timetree analyses.
1. **Safe** Subfolder containing the results for the Safe timetree analyses.


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

TBC
   
## Checkpoint for reproducibility
To rerun all the code with packages as they existed on CRAN at time of our analyses we recommend using the `checkpoint` package, and running this code prior to the analysis:

```{r}
checkpoint("xxx")
```
