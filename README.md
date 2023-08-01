# pinniped-macro
Code and data for pinniped macroevolution paper

Author(s): Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper

*This README is a work in progess*

This repository contains all the code and some data used in the [paper](XXX). 

To cite the paper: 
>  Park T, Burin G, Lazo-Cancino D, Rees J.P.G., Rule J.P, Slater G.J, Cooper N. Macroevolutionary patterns of pinnipeds. in prep.

To cite this repo: 
>  Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper. Code for the paper v1.

![alt text](https://github.com/nhcooper123/pinniped-macro/raw/main/figures/example.png)

This paper has three main components, a metatree of Pinnipedia, biogeographical analyses and diversification analyses.

-------
## Metatree

Code for the metatree is within the `metatree-scripts/` folder. Note that the data etc. are too large for GitHub so are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq) within the `Pinniped_metatree` folder. The file structure of that repository is include below for completeness.


## InputData

Folder containing the following subfolders:

### Molecular_Data

Folder containing two subdirectories (fulton and Lopes), each of which contains a further two subdirectories (mrp and xml) that contain the Matrix Representations with Parsimony (MRP) nexus files and XML files (see below), respectivly, for the two molecular studies used as backbone constraints

### MRP

Matrix Representations with Parsimony (MRPs) obtained by re-analysing each morphological input data set under parsimony and retaining only unique biparitions until all such biparitions have been sampled.

### NEXUS

The original character-taxon matrices for each morphological input data set.

### TNT

TNT format input files used for reanalysis of each morphological input dataset. Each file contains the relevant commands to run the analyses.

### XML

XML files that record important metadata about each input data set, particularly data set dependence and taxonomic reconciliation.




## metatree_files

Folder containing the metatree data organised in the following subfolders:

### Consensus_trees

Folder containing Strict (SCC) and Majority Rule (MRC) trees before (STR_***) and after (no prefix) safely reinserting taxa

### Files

All files output by running the metatree function to assemble a reconciled complete MRP matrix over all input datasets:

#### Character weights 

weights for each character

#### DataSetWeights.txt 

information of the contribution of each dataset to the overall analysis 

#### FULL.nex / FULL.tnt

MRP files in nexus and tnt format containing all taxa (not used in main analysis)

#### STR.nex / STR.tnt

Safely reduced input files in nexus and tnt format. The tnt format file is the one used in main analyses.

### STR.tnt

text file containing safely removed taxa and the rules used to reinsert them

#### TaxonomyTreee.tre

newick string containing the paleobiology database taxonomy for PanPinnipeds in tree form. Used in the metatree analysis as a very downweighted input.


### MPTS

Folder containing all most parsimonious trees before and after Safe Taxonomic Reinsertion

### tnt

MPTs from each of 1000 independent TNT searches (.tnt) and associated search info from screen output (.txt) 


### README.md

This file.

### Scripts

The set of R scripts used to generate the metatree

#### BuildBeastXML.R

A function used to assemble a nexus file, starting tree, topology constraints, and age priors for BEAST input. Is sourced by the file CreateBeastInput_Script.R

#### Build_Metatree_Files_Combined.R

Script to read in input files from input_files folder and assemble MRP matrices for generation of metatree

#### CreateBeastInput_Script.R

Script to read in a set of input files contained in the TimeTree folder and output a set of files that can be used to build a BEAST XML for Fossilized Birth Death tip dating analyses with a relaxed molecular clock

#### makemrpandscandmpts.R

Utility script to read the output of TNT for each source study, produce separate MPT and Strict Cosensus files, and generate an MRP file that can be used for input into the metatree function

#### process_beast_output.R

Script to read combined, post-burnin output from beast and find the median tree under the Kendall-Colijn distance metric. This tree is used in subsequent macroevolutionary analyses.

#### runtntoutputcollator.R

Utility script to read the output from TNT analysis of the ensemble MRP matrix, find all unique shortest trees, and assemble a combined MPT file, consensus trees, and perform safe taxonomic reinsertion. 

### TimeTree

Files and Folders related to BEAST time-tree inference. Organized as follows:

#### BEAST_Run

Folder containing input xml file, log and tree files from the two independent MCMC analyses, combined post-burnin log and tree files, and the median tree from the combined posterior sample.

#### partitionFinder

PartitionFinder input and output files used to identify the appropriate partitioning scheme and evolutionary models for the augmented Fulton and Strobeck alignment. 

#### strictclockrun

xml and output files from a short analysis of extant pinnipeds only under a strict molecular clock. The posterior mean root height divided by the approximate age of crown pinnipeds was used to define a mean for the lognormal prior on the relaxed molecular clock in the main analysis

#### MRC.tre

The majority Rule Consensus metatree used to derive topological constraints for the BEAST analysis

#### pinnipeds_extant_beast.nex

augmented version of the Fulton and Strobeck alignment used in BEAST analysis

#### Taxa_FAD.csv

Stratigraphic uncertainty associated with the first appearances of fossil taxa. Used to define age ranges for fossils for the Fossilized Birth Death analysis.

### TimeTreeInference

Folder containing the R scripts and temporal data used to build the timetrees as well as the following subfolders:

#### Dangerous

Subfolder containing the results for the Dangerous timetree analyses.

#### Risky

Subfolder containing the results for the Risky timetree analyses.

#### Safe

Subfolder containing the results for the Safe timetree analyses.


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

* The folder `biogeography/outputs` contains the main figures and tables. ESM figures/tables are in the higher level `/supplemental` folder.

-------

## Diversification rate analyses

Code for the diversification rate analyses is within the `diversification/` folder

### Data
All data are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq), but also within the  `data/`  and `raw-data/` folders.

If you use the data please cite as follows: 
>  Travis Park; Gustavo Burin; Graham J Slater; Natalie Cooper (2022). Data from the "Back to the water" project [Data set]. Natural History Museum. 

### Analyses

-------
## Other folders

* `/figures` contains the figures.
* `/outputs` contains the tables.
* `/img` contains the silhouettes from from `PhyloPic.org` needed for plotting. Contributed by: XXX
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
