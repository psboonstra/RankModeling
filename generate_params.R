require(tidyverse);

short_scenarios <- 
  crossing('0' = c(-1, 0.5),
           #
           nesting('1' = c(0.09, 1.5, 1.5), 
                   '2' = c(0.08, 1.4, 1.5), 
                   '3' = c(0.07, 1.3, 1.5),
                   '4' = c(0.06, 1.2, 1.5), 
                   '5' = c(0.05, 1.1, 1.5),
                   '6' = c(0, 1.0, 0),
                   '7' = c(0, 0.9, 0), 
                   '8' = c(0, 0.8, 0), 
                   '9' = c(0, 0.7, 0)),
           #
           nesting('10' = c(0)), 
           nesting('log_delta1' = c(log(0.9)), 
                   'log_delta2' = c(log(0.9)))) %>%
  filter((`1` > 0) |
           (`6` == 0)) %>%
  arrange(`2`,`0`);

long_scenarios <- 
  crossing('0' = c(-1, 0.5),
           #
           nesting('1' = c(0.09, 1.5, 1.5), 
                   '2' = c(0.08, 1.4, 1.5), 
                   '3' = c(0.07, 1.3, 1.5),
                   '4' = c(0.06, 1.2, 1.5), 
                   '5' = c(0.05, 1.1, 1.5),
                   '6' = c(0.04, 1.0, 1.5),
                   '7' = c(0.03, 0.9, 1.5), 
                   '8' = c(0.02, 0.8, 1.5), 
                   '9' = c(0.01, 0.7, 0),
                   '10' = c(0, 0.6, 0),
                   '11' = c(0, 0.5, 0),
                   '12' = c(0, 0.4, 0),
                   '13' = c(0, 0.3, 0),
                   '14' = c(0, 0.2, 0),
                   '15' = c(0, 0.1, 0)),
           #
           nesting('16' = c(0),
                   '17' = c(0),
                   '18' = c(0),
                   '19' = c(0),
                   '20' = c(0)), 
           nesting('log_delta1' = c(log(0.9)), 
                   'log_delta2' = c(log(0.9)))) %>%
  filter((`1` > 0) |
           (`6` == 0)) %>%
  arrange(`2`,`0`);

if(!"n_training_seq"%in%ls()){n_training_seq = c(30, 100, 500);}

if(!"array_id"%in%ls()){array_id = 1;}
if(!"array_id_offset"%in%ls()){array_id_offset = 0;}
if(!"n_sim"%in%ls()) {n_sim = 25;}

arglist = list();
scenario_id = 0;
for(j in 1:length(n_training_seq)) {
  for(i in 1:nrow(short_scenarios)) {
    scenario_id = scenario_id + 1;
    
    assign(paste("sim",array_id + array_id_offset,"_params",sep=""),
           list(array_id = array_id + array_id_offset,
                scenario_id = scenario_id,
                param_id = i,
                n_sim = n_sim,
                true_param = short_scenarios %>% slice(i) %>% unlist(use.names = T),
                n_training = n_training_seq[j],
                random_seed = sample(.Machine$integer.max - 1e4, 1), 
                data_seeds = NULL))
    arglist = c(arglist,list(get(paste("sim",array_id + array_id_offset,"_params",sep=""))));
  }
}
for(j in 1:length(n_training_seq)) {
  for(i in 1:nrow(long_scenarios)) {
    scenario_id = scenario_id + 1;
    
    assign(paste("sim",array_id + array_id_offset,"_params",sep=""),
           list(array_id = array_id + array_id_offset,
                scenario_id = scenario_id,
                param_id = nrow(short_scenarios) + i,
                n_sim = n_sim,
                true_param = long_scenarios %>% slice(i) %>% unlist(use.names = T),
                n_training = n_training_seq[j],
                random_seed = sample(.Machine$integer.max - 1e4, 1),
                data_seeds = NULL))
    arglist = c(arglist,list(get(paste("sim",array_id + array_id_offset,"_params",sep=""))));
  }
}
