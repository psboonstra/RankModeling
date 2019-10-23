


library(tidyverse);
library(TopKLists);
library(pROC);
library(ExtMallows);
library(gtools);


# DESCRIPTION: 
# Flag for whether this is running on a local machine or on a cluster running SLURM
my_computer = T;

# 'which_run' indicates whether this is the first batch (= 1 ) 
# or the second batch ( = 2)
which_run = 1;

if(my_computer) {
  array_id = 1;
} else {
  #srun --time=1:00:00 --cpus-per-task=2 --mem=4G --pty /bin/bash
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
}

#Helpful to print warnings when they occur for debugging
options(warn = 1);

if(which_run == 1) {#first batch
  
  # 'jobs_per_scenario' is the number of parallel independent jobs to send for 
  # each scenario, and 'n_sim' is the number of independent datasets per job
  # to generate. The constraints are as follows: 
  # (i) jobs_per_scenario * length(arglist) < 1000 (because only 1000 jobs
  # can be submitted at once)
  # (ii) sum_{i=1}^{jobs_per_scenario} n_sim_i = total desired sims per scenario
  jobs_per_scenario = ifelse(my_computer, 
                             27, 
                             27);
  n_sim = ifelse(my_computer, 
                 1,# ifelse(array_id <= 36, 12, 38), #
                 ifelse(array_id <= 36, 12, 38));
  # 'array_id_offset' is added to the label of the saved object. Useful when a 
  # new batch of jobs is run and you want to continue the labeling scheme. 
  array_id_offset = 0;
  # Should the sim_id labels be randomly permuted across array_ids?
  permute_array_ids = ifelse(my_computer, T, T);
} else { 
  jobs_per_scenario = ifelse(my_computer, 
                             27, 
                             27);
  n_sim = ifelse(my_computer, 
                 4,# ifelse(array_id <= 36, 12, 38), #
                 ifelse(array_id <= 36, 12, 38));
  # 'array_id_offset' is added to the label of the saved object. Useful when a 
  # new batch of jobs is run and you want to continue the labeling scheme. 
  array_id_offset = 972;
  # Should the sim_id labels be randomly permuted across array_ids?
  permute_array_ids = ifelse(my_computer, T, T);
}


# Before calling the next line, the user should specify any parameters that 
# she/he wishes to change from their default values. Any specified values will
# overwrite their default values, which are set by sourcing generate_params.R 
source("generate_params.R");

rm(list=setdiff(ls(all=T),c("arglist","array_id","my_computer",
                            "permute_array_ids","jobs_per_scenario",
                            "array_id_offset")));

source("functions_bpl.R");
source("functions_simulations.R");
source("functions_ldrbo.R");


#This is the step that permutes the duplicated job ids (the first are left alone)
#If FALSE, then all jobs from the same scenario will occur in contiguous blocks. 
if(permute_array_ids) {
  permute_array_ids = seq(1, jobs_per_scenario * length(arglist), by = jobs_per_scenario);
  set.seed(2);#set.seed(1);
  permute_array_ids = 
    c(permute_array_ids, 
      setdiff(sample(jobs_per_scenario * length(arglist)), permute_array_ids));
  
} else {
  permute_array_ids = 1:(jobs_per_scenario * length(arglist));
}

curr_args = arglist[[ceiling(permute_array_ids[array_id]/jobs_per_scenario)]];
curr_args[["random_seed"]] = curr_args[["random_seed"]] + array_id_offset + array_id;

# Phil delete below ----
curr_args$true_param = curr_args$true_param[-1];
# Phil delete above ----

assign(paste0("sim",array_id_offset + array_id),do.call("simulator",args = curr_args));


do.call("save",list(paste0("sim",array_id_offset + array_id),
                    file = paste0("out/sim",array_id_offset + array_id,".RData")));

write_csv(get(paste0("sim",array_id_offset + array_id))$summarized_performance,
          path = paste0("out/sim",array_id_offset + array_id,"_performance.csv"));

write_csv(get(paste0("sim",array_id_offset + array_id))$summarized_bpl,
          path = paste0("out/sim",array_id_offset + array_id,"_bpl.csv"));
