# Penalized multistage models for ordered data

### Current Suggested Citation

Boonstra, Philip S. and Krauss, John C., "Inferring a consensus problem list 
using penalized multistage models for ordered data" (March 2019) *The 
University of Michigan Department of Biostatistics Working Paper Series*.
Working Paper 126.
https://biostats.bepress.com/umichbiostat/paper126


[![DOI](https://zenodo.org/badge/175481636.svg)](https://zenodo.org/badge/latestdoi/175481636)


## Executive Summary
The function <samp>penRank_path</samp> contained in the file
<samp>functions_bpl.R</samp> represents the primary statistical contribution 
from this manuscript. This function will estimate the solution path for a 
Benter-Plackett-Luce model penalized with seamless <MATH>L<sub>0</sub> </MATH> 
penalties


## Further details

In more detail, there are eleven files included in this repository (in addition 
to this README and the authors' version of the manuscript): three CSV files 
(ending in <samp>.csv</samp>) and eight <samp>R</samp> scripts (ending in 
<samp>.R</samp>). The results reported in the manuscript were run using commit 13.

### CSV files

<samp>case</samp>X<samp>_20Dec2014.csv</samp>, where X = 23, 83, and 111 
contain the problem list data as ranked lists. Case 23 = Case A; 
Case 111 = Case B; Case 83 = Case C

### <samp>R</samp> files

<samp>functions_bpl.R</samp> provides all of the necessary functions to fit 
the BPL methods described in the paper

<samp>functions_ldrbo.R</samp> provides all of the necessary functions to 
calculate the consensus LDRBO reported in this paper and Krauss, et al. (2015)

<samp>gather_data.R</samp> reads in the problem list data from the .csv files 
and recharacterizes from ranked lists to ordered lists. Case 23 == Case A; 
Case 111 == Case B; Case 83 == Case C

<samp>fit_model_problists.R</samp> calls the previous two scripts and then 
calculates the solution paths and creates the tables and figures reported 
in the manuscript

<samp>run_rank_sims.R</samp> is the top-level script for conducting the 
simulation study. You provide the value of <samp>array_id</samp> on line 13
or the SLURM scheduler provides the value on line 15 to run the simulation 
scenario that you want to run. Choose any integer from 1 to 36 (actually you can
choose any positive integer and it will be mapped to the numbers 1 to 36 via
modular arithmetic). **If you run this script once for each of array_id = 
1, ..., 972, you will have run the entire simulation study reported in this
manuscript**

<samp>generate_params.R</samp> is called by <samp>run_rank_sims.R</samp> to 
create the individual data generating mechanisms used in the simulation study

<samp>function_simulations.R</samp> is also called by 
<samp>run_rank_sims.R</samp> and contains all of the code for doing the 
simulation study

<samp>process_results.R</samp> should be called after the simulations are
complete. Each time the script <samp>run_rank_sims.R</samp> is run, two 
files will be saved: simX_performance.csv and simX_bpl.csv, where X = the 
array_id used. Assuming these are all collected in a folder called 'out', running
<samp>process_results.R</samp> will collate all of these results together and 
create the tables and figures reported in the manuscript

## Acknowledgments 

This work was supported by the National Institutes of Health (UL1TR002240)