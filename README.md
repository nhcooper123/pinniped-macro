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

Code for the metatree is within the `metatree/` folder.

### Data
All data are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq), but are too large to upload to GitHub.

### Analyses


-------
## Biogeography analyses

Code for the biogeographical analyses is within the `biogeography/` folder

### Data
All data are available from the [NHM Data Portal](https://doi.org/10.5519/vmbrpkuq), but also within the  `data/`  and `raw-data/` folders.

If you use the data please cite as follows: 
>  Travis Park; Gustavo Burin; Graham J Slater; Natalie Cooper (2022). Data from the "Back to the water" project [Data set]. Natural History Museum. 

### Analyses

R-scripts 01-05 are modified versions of BioGeoBEARS R script(s) that can be found at [http://phylo.wikidot.com/biogeobears](http://phylo.wikidot.com/biogeobears). Copyright Nicholas J. Matzke. Please cite BioGeoBEARS if you use these scripts.

1. **00-prepare-trees** removes species without geographical information from the trees.
1. **01-pinniped-state-list-fix.R** modifies the states allowed for each analysis. THIS NEEDS TO BE RUN BEFORE SCRIPT 02. It is sourced in most other scripts to prevent you forgetting.
1. **02A-fit-DEC-DECj** fits DEC and DECj models using BioGeoBEARS for all taxa removing impossible states.
1. **02B-fit-DEC-DECj-unlikely** fits DEC and DECj models using BioGeoBEARS for all taxa removing impossible and unlikely states.
1. **03A-fit-DEC-DECj_fossil** fits DEC and DECj models using BioGeoBEARS for fossil taxa removing impossible states.
1. **03B-fit-DEC-DECj-unlikely_fossil** fits DEC and DECj models using BioGeoBEARS for fossil taxa removing impossible and unlikely states.
1. **04A-fit-DEC-DECj_extant** fits DEC and DECj models using BioGeoBEARS for extant taxa removing impossible states.
1. **04B-fit-DEC-DECj-unlikely_extant** fits DEC and DECj models using BioGeoBEARS for extant taxa removing impossible and unlikely states.
1. **05-DEC-DECj-results.R** extracts results from DEC and DECj models for tables. 
1. **06-DEC-DECj-plots.R** creates all BioGeoBEARS plots for the ESM.
1. **07-Extract-BGB-results-for-plotting.R** extracts results from BioGeoBEARS results objects for making the main publication plot.
1. **08-main-results-plot.R** code for the nice main publication plot.
1. **09-worldmap-for-legend.R** quick snippet to get a world map for the legend.
1. **messing-about-with-colours.R** attempts to set nice colour palettes for scripts 06, 07 and 08 which call it.
1. **modifiedBGBplot.R** modified version of BioGeoBEARS plotting function to have more control over plots in 06.

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
