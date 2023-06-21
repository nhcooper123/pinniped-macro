# pinniped-macro
Code and data for pinniped macroevolution paper

Author(s): Gustavo Burin, Travis Park, Graham Slater, Natalie Cooper

*This README is a work in progess*

This repository contains all the code and some data used in the [paper](XXX). 

To cite the paper: 
>  Park T, Burin G, Lazo-Cancino D, Rees J.P.G., Rule J.P, Slater G.J, Cooper N. Macroevolutionary patterns of pinnipeds. in prep.

To cite this repo: 
>  

![alt text](https://github.com/nhcooper123/pinniped-macro/raw/main/figures/example.png)

This paper has three main components, a metatree of Pinnipedia, biogeographical analyses and diversification analyses.

-------
## Metatree

Code for the metatree is within the `metatree/` folder

-------
## Biogeography analyses

Code for the biogeographical analyses is within the `biogeography/` folder

### Data
All cleaned data are available from the [NHM Data Portal](), but also within the  `data/`  and `raw-data` folders/.

If you use the cleaned data please cite as follows: 
>  

### Analyses

R-scripts 01 and 02 are modified versions of BioGeoBEARS R script(s) that can be found at [http://phylo.wikidot.com/biogeobears](http://phylo.wikidot.com/biogeobears). Copyright Nicholas J. Matzke. Please cite BioGeoBEARS if you use these scripts.

1. **00-prepare-trees** removes species without geographical information from the trees.
1. **01-pinniped-state-list-fix.R** modifies the states allowed for each analysis. THIS NEEDS TO BE RUN BEFORE SCRIPT 02.
1. **02-fit-DEC-DECj** fits DEC and DECj models using BioGeoBEARS.

-------

## Diversification rate analyses

Code for the diversification rate analyses is within the `diversification/` folder

### Data
All cleaned data are available from the [NHM Data Portal](), but also within the  `data/`  and `raw-data` folders/.

If you use the cleaned data please cite as follows: 
>  

### Analyses

-------
## Other folders

* `/figures` contains the figures.
* `/outputs` contains the tables.
* `/img` contains the silhouettes from from `PhyloPic.org` needed for plotting. Contributed by: XXX
* `/manuscript` contains the manuscript materials in LaTeX format.

-------
## Session Info
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication.

TBC
   
## Checkpoint for reproducibility
To rerun all the code with packages as they existed on CRAN at time of our analyses we recommend using the `checkpoint` package, and running this code prior to the analysis:

```{r}
checkpoint("xxx")
```
