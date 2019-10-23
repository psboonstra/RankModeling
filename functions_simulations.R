
comparison_matrix = function(x) {
  foo = matrix(x, nrow = length(x), ncol = length(x));
  t_foo = t(foo);
  (foo >= t_foo) -
    (foo < t_foo) - 
    diag(1, length(x))
}

validate_weights = function(est_weights, 
                            true_weights, 
                            true_weights_comparison_matrix, 
                            tie_breaker = 100) {
  
  n_weights = length(est_weights);
  if(missing(true_weights_comparison_matrix)) {
    true_weights_comparison_matrix = comparison_matrix(true_weights);
  } 
  
  #RMSE ordering
  if(!any(duplicated(est_weights)) && !any(duplicated(true_weights))) {
    tie_breaker = 0L;# no tie breaker needed if everything is unique
  }
  stopifnot((length(tie_breaker) == 1L && 
               tie_breaker %% 1 == 0L && 
               tie_breaker >= 0) ||
              (length(tie_breaker) == length(est_weights)));
  
  if(length(tie_breaker) == 1L && tie_breaker == 0) {
    # Take everything as is, ignoring any ties
    rmse_order =sqrt(mean((true_weights[order(true_weights,
                                              decreasing = TRUE)] -  
                             true_weights[order(est_weights,
                                                decreasing = TRUE)])^2));
  } else if(length(tie_breaker) == 1L && tie_breaker > 0) {
    # Randomly permute ties and take the average
    foo = 0;
    for(k in 1:tie_breaker) {
      foo = foo + 
        sqrt(mean((true_weights[order(true_weights,
                                      decreasing = TRUE)] -  
                     true_weights[order(est_weights,
                                        sample(n_weights), 
                                        decreasing = TRUE)])^2));
      
    }
    rmse_order = foo / tie_breaker;
  } else {
    # This corresponds to tie_breaker being a vector 
    rmse_order = sqrt(mean((true_weights[order(true_weights,
                                               decreasing = TRUE)] -  
                              true_weights[order(est_weights,
                                                 tie_breaker,
                                                 decreasing = TRUE)])^2));
    
  }
  
  # Emond and Mason's tau x
  tau_x = sum(comparison_matrix(est_weights) * 
                true_weights_comparison_matrix) / n_weights / (n_weights - 1);
  
  list(tau_x = tau_x, 
       rmse_order = rmse_order);
}


# function name and purpose: 'simulate_bpl' Simulates a set of lists according to the BPL model implied by 'param'. 
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 3/5/2019, 10:00am EST
#
# Function inputs:
#
### param: (vector with 'names' attribute, numeric) the parameter vector to evaluate with respect to the data. Must be named exactly 
# as follows: '0' (optional, if the stopping process is to be modeled), then a set of item labels, then 'log_delta1', then 'log_delta2'. 
# In contrast to other functions, the item labels do not have to be complete nor consecutive. You would leave out '7', for example, if you
# are simulating from stage 2 and item '7' was put first at the first stage, so that you want the stage 2 probabilities conditional upon
# stage 1 having ranked item 7. 
### stages: (vector, consecutive positive integers starting with 1) up to how many stages should be simulated? OR (single positive integer)
# which stage to simulate?
#
### n_to_sim: (integer, positive) how many lists to generate?
#
### Value
#
# If length(stages) == 1, a vector of items, corresponding to n_to_sim draws from the conditional multinomial distribution from that BPL model
# at that stage. If length(stages) > 1, a matrix with dimensions c(n_to_sim, length(stages)) corresponding to a n_to_sim draws from the sequence
# of conditional multinomial distributions from that BPL model. 

simulate_bpl = function(param, stages = 1, n_to_sim = 1) {
  if(length(stages) == 1) {
    result = sample(setdiff(names(param), c("log_delta1","log_delta2")),size = n_to_sim,  prob = bpl_probs(param, stage = stages), replace = TRUE)
  } else {
    stopifnot(all(stages == 1:length(stages)));
    result = matrix(NA, nrow = n_to_sim, ncol = length(stages));
    result[,1] = sample(as.numeric(setdiff(names(param), c("log_delta1","log_delta2"))),
                        size = n_to_sim,  
                        prob = bpl_probs(param, stage = 1), 
                        replace = TRUE);
    for(curr_sim in 1:n_to_sim) {
      for(curr_stage in 2:length(stages)) {
        curr_param = param[setdiff(names(param),result[curr_sim,1:curr_stage])];
        result[curr_sim,curr_stage] = 
          sample(as.numeric(setdiff(names(curr_param), c("log_delta1","log_delta2"))),
                 size = 1,  
                 prob = bpl_probs(curr_param, stage = curr_stage));
        if(result[curr_sim,curr_stage] == 0) {
          result[curr_sim,curr_stage] = NA;
          break;
        }
      }
    }
  }
  if(any((foo <- colSums(is.na(result))) == n_to_sim)) {
    if(min(which(foo == n_to_sim)) == 1) {
      results[,1, drop = F];
    } else {
      result[,1:(min(which(foo == n_to_sim)) - 1), drop = F]; 
    }
  } else {
    result;
  }
}


simulator = function(array_id = 1,
                     scenario_id = 1,
                     param_id = 1,
                     n_sim = 1,
                     #Set to NA to do a dry run of the simulator:
                     true_param = NA,
                     n_training = 30,
                     random_seed = sample(.Machine$integer.max,1),
                     data_seeds = NULL)#If non-null, should have length equal to n_sim
{ 
  require(tidyverse);
  require(TopKLists);
  require(pROC);
  require(ExtMallows);
  require(gtools);
  
  begin_all = Sys.time();
  set.seed(random_seed);
  if(!length(data_seeds)) {
    data_seeds = sample(.Machine$integer.max,n_sim);
  } else {
    if(length(data_seeds) != n_sim) {
      stop("'data_seeds' was provided but did not have length equal to 'n_sim'");
    }  
  }
  
  
  if(any(is.na(true_param))) {
    fit_methods = F;
    cat("One or more missing value(s) in 'true_param' so I can't run the simulation\n");
  } else {
    fit_methods = T;
    n_items = length(true_param) - 3;
    item_names = as.character(1:n_items);
    if(!(length(true_param) > 3) && !isTRUE(all.equal(names(true_param),c(0:n_items,"log_delta1","log_delta2"))) && !isTRUE(all.equal(names(true_param),c(1:n_items,"log_delta1","log_delta2")))) {
      stop("'true_param' must be named exactly as follows: '0', then one or more consecutive positive integer labels starting with '1', one for each unique item, then 'log_delta1', then 'log_delta2'");
    }
    if(!near(min(true_param[item_names]), 0)) {
      stop("The minimum item weight should be zero");
    }
    if("log_delta1"%in%names(true_param) && true_param["log_delta1"] > 0) {stop("'log_delta1' must be non-positive");}
    if("log_delta2"%in%names(true_param) && true_param["log_delta2"] > 0) {stop("'log_delta2' must be non-positive");}
    
    null_item_names = item_names[which(near(true_param[item_names], 0))];
    nonnull_item_names = item_names[which(!near(true_param[item_names], 0))];
    
  }
  
  method_names = 
    c("bpl",
      "bpl_aic",
      "bpl_bic",
      "ldrbo",
      "mc1",
      "mc2",
      "mc3",
      "cemc_s",
      "cemc_k",
      "EMM",
      "MM") %>% 
    factor() %>% 
    fct_inorder();
  
  summarized_performance = 
    crossing(crossing(sim_id = 1:n_sim,
                      array_id = array_id,
                      scenario_id = scenario_id, 
                      param_id = param_id, 
                      n_training = n_training,
                      n_obs_items = NA_real_,
                      n_obs_nonnull_items = NA_real_,
                      method_names, 
                      # tau_x
                      tau_x = NA_real_,
                      # auc 
                      auc = NA_real_,
                      rmse_order = NA_real_,
                      #running time to fit each method
                      run_time_secs = NA_real_));
  
  summarized_bpl = 
    crossing(crossing(sim_id = 1:n_sim,
                      array_id = array_id,
                      scenario_id = scenario_id, 
                      param_id = param_id, 
                      n_training = n_training,
                      n_obs_items = NA_real_,
                      n_obs_nonnull_items = NA_real_,
                      method_names = c("bpl","bpl_aic","bpl_bic") %>% 
                        factor() %>%
                        fct_inorder(),
                      neg_loglikelihood  = NA_real_,
                      aic = NA_real_,
                      lambda = NA_real_,
                      num_eff_params = NA_real_,
                      num_actual_params = NA_real_,
                      max_item_weight = NA_real_,
                      rmse = NA_real_,
                      # This will be different from rmse_order used in
                      # 'summarized_performance' because that randomly breaks
                      # ties whereas this will use additional information 
                      # from the solution path to determine how tied items
                      # are ranked higher and lower
                      rmse_order = NA_real_,
                      tpr = NA_real_, 
                      tnr = NA_real_));
  
  
  begin_sim = Sys.time();
  i=1;
  
  if(fit_methods) {
    
    # + Draw data ----
    for(i in 1:n_sim) {
      set.seed(data_seeds[i]);
      curr_true_param = true_param;
      
      cat("\n######################################\n");
      cat("# the generating BPL params are : \n");
      cat("# ")
      print(curr_true_param);
      cat("# the number of observations used for training is", n_training, "\n");
      cat("# the value of 'random_seed' is", random_seed,"\n");
      cat("# iteration:",i, "/", n_sim,"\n");
      cat("# the current data seed is",data_seeds[i],"\n");
      cat("######################################\n\n");
      
      set.seed(data_seeds[i]);
      dat_ordered <-
        simulate_bpl(param = curr_true_param, 
                     stages = seq_len(length(curr_true_param) - 3), 
                     n_to_sim = n_training);
      
      # It is not guaranteed that all possible items will be observed, e.g. 
      # if n_training is small and there is a large a fatigue parameter, then 
      # we may never witness some items with small relative weights.
      # We make the assumption here that we are unaware of any such items that 
      # were unobserved in the data by relabeling the observed items to be in
      # contiguous order
      
      observed_item_names = 
        as.character(sort(unique(as.numeric(dat_ordered))));
      which_observed_nonnull = 
        which(observed_item_names %in% nonnull_item_names);
      
      # This iteration will only evaluate methods with respect to items that 
      # were actually observed
      curr_true_param = 
        curr_true_param[c("0", observed_item_names, "log_delta1","log_delta2")];
      names(curr_true_param) = 
        c("0", as.character(seq_along(observed_item_names)), "log_delta1","log_delta2");
      
      # Now if necessary relabel items so that they are in contiguous order
      if(length(observed_item_names) < length(item_names)) {
        for(j in seq_along(observed_item_names)) {
          dat_ordered[dat_ordered == observed_item_names[j]] = j;
        }
      } 
      curr_item_names = 
        as.character(seq_along(observed_item_names));
      curr_nonnull_item_names =
        curr_item_names[which_observed_nonnull];
      rm(observed_item_names,which_observed_nonnull);
      
      curr_n_items = length(curr_item_names);
      curr_n_nonnull_items = length(curr_nonnull_item_names);
      
      # For TopKList calculations
      dat_ordered_as_list = matrix_to_list(dat_ordered);
      
      # For Emond and Mason's TauX statistic:
      curr_true_param_comparison_matrix = 
        comparison_matrix(curr_true_param[curr_item_names]);
      nonnull_curr_true_param_comparison_matrix = 
        comparison_matrix(curr_true_param[curr_nonnull_item_names]);
      
      # For AUC
      can_do_auc = 
        any(curr_true_param[curr_item_names] > 0) && 
        any(curr_true_param[curr_item_names] == 0);
      
      # + BPL ----
      cat("#", "BPL","\n");
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "bpl"));
      curr_row_bpl = 
        with(summarized_bpl, 
             which(sim_id == i & method_names == "bpl"));
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_bpl[curr_row_bpl, "n_obs_items"] = 
        curr_n_items;
      summarized_bpl[curr_row_bpl, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      
      bpl_soln_path <-
        penRank_path(dat = dat_ordered, 
                     num_inits = 3, 
                     lambda_seq = 0,
                     consider_truncating_lambda_seq = FALSE,
                     min_reps = 2, 
                     stable_reps = 2,
                     max_reps_keq1 = 100, 
                     max_reps_kgt1 = 25,
                     verbose = FALSE, 
                     safe_mode = FALSE);
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(bpl_soln_path$runtime);
      
      selected_iter = which.min(bpl_soln_path$best_neg_loglikelihood["neg_loglikelihood",]);
      curr_weights = 
        bpl_soln_path$best_fit_neg_loglikelihood[curr_item_names,selected_iter];
      curr_tie_breaker = 
        bpl_soln_path$all_params[[selected_iter]][curr_item_names,1]
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      summarized_bpl[curr_row_bpl,"neg_loglikelihood"] =
        bpl_soln_path$all_neg_loglikelihood[1,1];
      summarized_bpl[curr_row_bpl,"aic"] =
        bpl_soln_path$all_aic[1,1];
      summarized_bpl[curr_row_bpl,"lambda"] =
        bpl_soln_path$control$lambda_seq[1];
      summarized_bpl[curr_row_bpl,"num_eff_params"] =
        bpl_soln_path$num_eff_params[1,1];
      summarized_bpl[curr_row_bpl,"num_actual_params"] =
        bpl_soln_path$num_actual_params[1,1];
      summarized_bpl[curr_row_bpl,"max_item_weight"] =
        max(curr_weights);
      summarized_bpl[curr_row_bpl,"rmse"] =
        sqrt(mean((curr_weights - curr_true_param[curr_item_names])^2));
      summarized_bpl[curr_row_bpl,"tpr"] =
        sum((!near(curr_weights, 0)) * (!near(curr_true_param[curr_item_names], 0))) / 
        sum(!near(curr_true_param[curr_item_names], 0));
      summarized_bpl[curr_row_bpl,"tnr"] =
        sum(near(curr_weights, 0) * near(curr_true_param[curr_item_names], 0)) / 
        sum(near(curr_true_param[curr_item_names], 0));
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix, 
                             tie_breaker = curr_tie_breaker);
      summarized_bpl[curr_row_bpl, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      
      init_params = bpl_soln_path$best_fit_neg_loglikelihood;
      rm(bpl_soln_path, curr_weights, curr_tie_breaker, curr_row_bpl, curr_row_performance, selected_iter);
      
      # + BPL (AIC) ----
      cat("#", "BPL (AIC)","\n");
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "bpl_aic"));
      curr_row_bpl = 
        with(summarized_bpl, 
             which(sim_id == i & method_names == "bpl_aic"));
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_bpl[curr_row_bpl, "n_obs_items"] = 
        curr_n_items;
      summarized_bpl[curr_row_bpl, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      
      bpl_soln_path <-
        penRank_path(dat = dat_ordered, 
                     num_inits = ncol(init_params), 
                     num_lambda = 50,
                     min_reps = 2, 
                     stable_reps = 2,
                     max_reps_keq1 = 100, 
                     max_reps_kgt1 = 25,
                     verbose = FALSE, 
                     safe_mode = FALSE,
                     init_params = init_params);
      
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(bpl_soln_path$runtime);
      
      selected_iter = which.min(bpl_soln_path$best_aic["aic",]);
      curr_weights = 
        bpl_soln_path$best_fit_aic[curr_item_names,selected_iter];
      curr_tie_breaker = 
        bpl_soln_path$all_params[[selected_iter]][curr_item_names,1]
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      summarized_bpl[curr_row_bpl,c("aic","neg_loglikelihood","lambda","num_eff_params","num_actual_params")]  = 
        bpl_soln_path$best_aic[c("aic","neg_loglikelihood","lambda","num_eff_params","num_actual_params"),1];
      summarized_bpl[curr_row_bpl,"max_item_weight"] =
        max(curr_weights);
      summarized_bpl[curr_row_bpl,"rmse"] =
        sqrt(mean((curr_weights - curr_true_param[curr_item_names])^2));
      summarized_bpl[curr_row_bpl,"tpr"] =
        sum((!near(curr_weights, 0)) * (!near(curr_true_param[curr_item_names], 0))) / 
        sum(!near(curr_true_param[curr_item_names], 0));
      summarized_bpl[curr_row_bpl,"tnr"] =
        sum(near(curr_weights, 0) * near(curr_true_param[curr_item_names], 0)) / 
        sum(near(curr_true_param[curr_item_names], 0));
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix, 
                             tie_breaker = curr_tie_breaker);
      summarized_bpl[curr_row_bpl, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      rm(curr_weights, curr_tie_breaker, curr_row_bpl, curr_row_performance, selected_iter, init_params);
      
      # + BPL (BIC) ----
      cat("#", "BPL (BIC)","\n");
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "bpl_bic"));
      curr_row_bpl = 
        with(summarized_bpl, 
             which(sim_id == i & method_names == "bpl_bic"));
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_bpl[curr_row_bpl, "n_obs_items"] = 
        curr_n_items;
      summarized_bpl[curr_row_bpl, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(bpl_soln_path$runtime);
      
      selected_iter = which.min(bpl_soln_path$best_bic["bic",]);
      curr_weights = 
        bpl_soln_path$best_fit_bic[curr_item_names,selected_iter];
      curr_tie_breaker = 
        bpl_soln_path$all_params[[selected_iter]][curr_item_names,1]
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      summarized_bpl[curr_row_bpl,c("bic","neg_loglikelihood","lambda","num_eff_params","num_actual_params")]  = 
        bpl_soln_path$best_bic[c("bic","neg_loglikelihood","lambda","num_eff_params","num_actual_params"),1];
      summarized_bpl[curr_row_bpl,"max_item_weight"] =
        max(curr_weights);
      summarized_bpl[curr_row_bpl,"rmse"] =
        sqrt(mean((curr_weights - curr_true_param[curr_item_names])^2));
      summarized_bpl[curr_row_bpl,"tpr"] =
        sum((!near(curr_weights, 0)) * (!near(curr_true_param[curr_item_names], 0))) / 
        sum(!near(curr_true_param[curr_item_names], 0));
      summarized_bpl[curr_row_bpl,"tnr"] =
        sum(near(curr_weights, 0) * near(curr_true_param[curr_item_names], 0)) / 
        sum(near(curr_true_param[curr_item_names], 0));
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix, 
                             tie_breaker = curr_tie_breaker);
      summarized_bpl[curr_row_bpl, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      rm(curr_weights, bpl_soln_path, curr_tie_breaker, curr_row_bpl, curr_row_performance, selected_iter);
      
      
      # + LDRBO ----
      cat("#", "LDRBO","\n");
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "ldrbo"));
      
      begin_ldrbo = Sys.time();
      ldrbo_results = 
        consensus_ldrbo(dat = dat_ordered, 
                        psi = 1, 
                        max_size = curr_n_items,
                        min_size = curr_n_items, 
                        look_beyond_init = 5,
                        look_beyond_final = 5,
                        verbose = FALSE);
      
      ldrbo_results = 
        consensus_ldrbo(dat = dat_ordered, 
                        psi = 1, 
                        max_size = curr_n_items,
                        min_size = curr_n_items, 
                        look_beyond_init = 
                          ifelse(n_training <= 100, 500, 100),
                        look_beyond_final =
                          ifelse(n_training <= 100, 100, 10),
                        initial_order = 
                          ldrbo_results$consensus_list, 
                        window_seq = 
                          rep(6, curr_n_items),
                        verbose = T);
      
      end_ldrbo = Sys.time();
      
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(difftime(end_ldrbo, begin_ldrbo, units = "secs"));
      
      ldrbo_ranks = order(c(ldrbo_results$consensus_list, 
                            setdiff(1:n_items,ldrbo_results$consensus_list)));
      ldrbo_ranks[ldrbo_ranks > length(ldrbo_results$consensus_list)] = (n_items + 1);
      
      curr_weights = -(ldrbo_ranks);
      names(curr_weights) = curr_item_names;
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      rm(begin_ldrbo,ldrbo_results, end_ldrbo, curr_row_performance, curr_weights);
      
      
      # + MC ----
      cat("#", "MC 1--3","\n");
      
      begin_mc = Sys.time();
      mc_results = MC(input = dat_ordered_as_list);
      end_mc = Sys.time();
      
      for(curr_mc_name in c("mc1","mc2","mc3")) {
        curr_row_performance = 
          with(summarized_performance, 
               which(sim_id == i & method_names == curr_mc_name));
        curr_weights = 
          mc_results[[paste0(toupper(curr_mc_name),".Prob")]][
            order(mc_results[[paste0(toupper(curr_mc_name),".TopK")]])] %>%
          round(9);
        names(curr_weights) = curr_item_names;
        
        summarized_performance[curr_row_performance, "n_obs_items"] = 
          curr_n_items;
        summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
          curr_n_nonnull_items;
        summarized_performance[curr_row_performance, "run_time_secs"] =
          as.numeric(difftime(end_mc, begin_mc, units = "secs"))/3;
        
        foo = validate_weights(est_weights = curr_weights, 
                               true_weights = curr_true_param[curr_item_names], 
                               true_weights_comparison_matrix = curr_true_param_comparison_matrix);
        summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
        summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
        rm(foo);
        
        if(can_do_auc) {
          summarized_performance[curr_row_performance, "auc"] =
            as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
        }
        
        rm(curr_row_performance, curr_weights);
      }
      rm(begin_mc,mc_results, end_mc, curr_mc_name);
      
      # + CEMC ----
      cat("#", "CEMC (Spearman, Kendall)","\n");
      
      begin_cemc = Sys.time();
      cemc_results = CEMC(dat_ordered_as_list, dm = "s");
      end_cemc = Sys.time();
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "cemc_s"));
      curr_weights = -order(cemc_results$TopK);
      names(curr_weights) = curr_item_names;
      
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(difftime(end_cemc, begin_cemc, units = "secs"));
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      rm(begin_cemc, end_cemc, cemc_results, curr_row_performance, curr_weights);
      
      begin_cemc = Sys.time();
      cemc_results = CEMC(dat_ordered_as_list, dm = "k");
      end_cemc = Sys.time();
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "cemc_k"));
      curr_weights = -order(cemc_results$TopK);
      names(curr_weights) = curr_item_names;
      
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(difftime(end_cemc, begin_cemc, units = "secs"));
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      rm(begin_cemc, end_cemc, cemc_results, curr_row_performance, curr_weights);
      
      # + EMM ----
      cat("#", "EMM","\n");
      
      begin_emm = Sys.time();
      EMM_results = EMM(rankings = t(replace_na(dat_ordered, 0)));
      end_emm = Sys.time();
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "EMM"));
      curr_weights = -order(as.numeric(EMM_results$op.pi0));
      names(curr_weights) = curr_item_names;
      
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(difftime(end_emm, begin_emm, units = "secs"));
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      rm(begin_emm, end_emm, EMM_results, curr_row_performance, curr_weights);
      
      # + MM ----
      cat("#", "MM","\n");
      
      begin_mm = Sys.time();
      MM_results = MM(rankings = t(replace_na(dat_ordered, 0)));
      end_mm = Sys.time();
      
      curr_row_performance = 
        with(summarized_performance, 
             which(sim_id == i & method_names == "MM"));
      curr_weights = -order(as.numeric(MM_results$op.pi0));
      names(curr_weights) = curr_item_names;
      
      summarized_performance[curr_row_performance, "n_obs_items"] = 
        curr_n_items;
      summarized_performance[curr_row_performance, "n_obs_nonnull_items"] = 
        curr_n_nonnull_items;
      summarized_performance[curr_row_performance, "run_time_secs"] =
        as.numeric(difftime(end_mm, begin_mm, units = "secs"));
      
      foo = validate_weights(est_weights = curr_weights, 
                             true_weights = curr_true_param[curr_item_names], 
                             true_weights_comparison_matrix = curr_true_param_comparison_matrix);
      summarized_performance[curr_row_performance, "tau_x"] = foo$tau_x;
      summarized_performance[curr_row_performance, "rmse_order"] = foo$rmse_order;
      rm(foo);
      
      if(can_do_auc) {
        summarized_performance[curr_row_performance, "auc"] =
          as.numeric(auc(1 * (curr_true_param[curr_item_names] > 0), curr_weights));
      }
      
      rm(begin_mm, end_mm, MM_results, curr_row_performance, curr_weights);
      
      rm(curr_true_param, curr_item_names, curr_nonnull_item_names, curr_n_items, 
         curr_n_nonnull_items, curr_true_param_comparison_matrix, 
         nonnull_curr_true_param_comparison_matrix, can_do_auc);
    }
  }
  
  list(summarized_performance = summarized_performance,
       summarized_bpl = summarized_bpl,
       true_param = true_param, 
       n_training = n_training,
       random_seed = random_seed,
       data_seeds = data_seeds,
       total_run_time_secs = difftime(Sys.time(), begin_all, units = "secs"));
}
