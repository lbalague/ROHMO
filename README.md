# ROHMO

## Summary

`ROHMO` is an R package for detecting Runs of Homozygosity and Mosaic Chromosome Alteration regions from BAF and LRR values.

## Installation

`ROHMO` requires R version equal or newer than 3.3.0. 

The package can be installed using the R package `devtools`. `devtools` can be installed win the following code:

```r
install.packages("devtools")
```

Once `devtools` and the dependences are installed, the following code installs `omicRexposome` and the basic dependence `rexposome`:

```r
devtools::install_github("isglobal-brge/ROHMO")
```

## Basic Guide

ROH and mCA detection is done using the function `ROHMOseeker`. This function requires an argument `sub` being the name of the sample to analyze and the arguments `baf`and `lrr`, being data frames with a minimum of 4 columns: First to third column should contain probe names, chromosome and position. The fourth column should contain the BAF/LRR values for the sample indicated in `sub`. The function returns a list with the ROH and mCA regions detected. If the `plot` argument is set TRUE, the function will also return the parameters needed for plotting the results.


Function `ROHMOplotter` allows to plot the results of `ROHMOseeker`. The plot is saved in tiff format in the directory specified with the argument `dir`.
