# Penalized multistage models for ordered data

### Current Suggested Citation

Boonstra, Philip S. and Krauss, John C., "Inferring a consensus problem list using penalized multistage models for ordered data" (March 2019) *The University of Michigan Department of Biostatistics Working Paper Series*. Working Paper 126.
https://biostats.bepress.com/umichbiostat/paper126


[![DOI](https://zenodo.org/badge/175481636.svg)](https://zenodo.org/badge/latestdoi/175481636)


## Executive Summary
The function <samp>penRank_path</samp> contained in the file <samp>Functions.R</samp> represents the primary statistical contribution from this manuscript. This function will estimate the solution path for a Benter-Plackett-Luce model penalized with seamless <MATH>L<sub>0</sub> </MATH> penalties

## Further details

In more detail, there are six files included in this repository (in addition to this README and the authors' version of the manuscript): three CSV files (ending in <samp>.csv</samp>) and three <samp>R</samp> scripts (ending in  <samp>.R</samp>). The results reported in the manuscript were run using commit 1.

### CSV files

<samp>case</samp>X<samp>_20Dec2014.csv</samp>, where X = 23, 83, and 111 contain the problem list data as ranked lists. Case 23 = Case A; Case 111 = Case B; Case 83 = Case C

### <samp>R</samp> files

<samp>functions.R</samp> provides all of the necessary functions to fit the methods described in the paper. 

<samp>gather_data.R</samp> reads in the problem list data from the .csv files and recharacterizes from ranked lists to ordered lists. Case 23 == Case A; Case 111 == Case B; Case 83 == Case C

<samp>fit_model.R</samp> calls the previous two scripts and then calculates the solution paths and creates the tables and figures reported in the manuscript. 

## Acknowledgments 

This work was supported by the National Institutes of Health (UL1TR002240)