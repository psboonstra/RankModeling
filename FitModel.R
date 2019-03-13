library(tidyverse);
library(RColorBrewer);
source("Functions.R");
source("gatherData.R");

fixed = NULL;
num_inits = 1;
conv_tol = 1e-3;
tau = conv_tol;
num_lambda = 250; 
lambda_min_ratio = 1e-6;
lambda_seq = NULL;
penalty_pow = 1;
proposal_seq = NULL;
proposal_seq_half_length = NULL;
max_proposal_seq = 1;
stable_reps = 5;
min_reps = 5;
multivariable_proposals = T;
delta_lb = 0.25;
max_reps_keq1 = 1e3;
max_reps_kgt1 = 250;
random_seed = sample(2^30.9, 1);
verbose = F;
safe_mode = F;
init_params = NULL;

for(case_name in c("23","111","83")) {
  cat(case_name, "\n\n");
  
  if(case_name == "23") {
    #Case 23 
    dat_ordered = case23_ordered;
    dat_ranked = case23_ranked;
    table_path = "../../case23_all_params.txt";
    figure_path = "../../case23_soln_path.pdf";
    
    case_color = "#800026"
    ldrbo = c(NA,1, 3, NA, 4, NA, 2, 6, 5, NA,NA, NA,NA, 7, rep(NA,13), 8,NA, NA, NA)
    unpenalized_soln_name = "case23_unpenalized_soln";
    penalized_soln_path_name = "case23_penalized_soln_path";
  } else if(case_name == "111") {
    #Case 111
    dat_ordered = case111_ordered;
    dat_ranked = case111_ranked;
    table_path = "../../case111_all_params.txt";
    figure_path = "../../case111_soln_path.pdf";
    
    case_color = "#00441B"
    ldrbo = c(NA, 2, 4, NA, NA, NA, 8, NA, 14, 12, NA, 11, NA, 1, 10, 3, 9, 7, 5, NA, 6, NA, NA, NA, 13, rep(NA, 19), 15, NA, NA, NA, NA, NA);
    unpenalized_soln_name = "case111_unpenalized_soln";
    penalized_soln_path_name = "case111_penalized_soln_path";
  } else if(case_name == "83") {
    #Case 83
    dat_ordered = case83_ordered;
    dat_ranked = case83_ranked;
    table_path = "../../case83_all_params.txt";
    figure_path = "../../case83_soln_path.pdf";
    
    case_color = "#00274C"
    ldrbo = c(NA, NA, 1, 6, NA, 2, NA, 3, NA, 5, 7, NA, 4, rep(NA,20 ));
    unpenalized_soln_name = "case83_unpenalized_soln";
    penalized_soln_path_name = "case83_penalized_soln_path";
  }
  
  ################
  obs_weights = rep(1, nrow(dat_ordered));
  
  #Rprof("penRank_run1.out",line.profiling = T)
  unpenalized_soln = penRank_path(dat = dat_ordered, obs_weights = obs_weights, fixed = fixed, num_inits = num_inits, 
                                  conv_tol = conv_tol, tau = 100, num_lambda = 2, lambda_min_ratio = 1e-10, lambda_seq = lambda_seq, consider_truncating_lambda_seq = FALSE,
                                  penalty_pow = penalty_pow,
                                  proposal_seq = proposal_seq, proposal_seq_half_length = proposal_seq_half_length, max_proposal_seq = max_proposal_seq,
                                  stable_reps = stable_reps, min_reps = min_reps, multivariable_proposals = multivariable_proposals,
                                  max_reps_keq1 = max_reps_keq1, max_reps_kgt1 = max_reps_kgt1,
                                  delta_lb = delta_lb,
                                  random_seed = random_seed, verbose = verbose, safe_mode = safe_mode);
  assign(unpenalized_soln_name,unpenalized_soln);
  penalized_soln_path = penRank_path(dat = dat_ordered, obs_weights = obs_weights, fixed = fixed, num_inits = num_inits, 
                                     conv_tol = conv_tol, tau = tau, num_lambda = num_lambda, lambda_min_ratio = lambda_min_ratio, lambda_seq = lambda_seq, consider_truncating_lambda_seq = TRUE,
                                     penalty_pow = penalty_pow,
                                     proposal_seq = proposal_seq, proposal_seq_half_length = proposal_seq_half_length, max_proposal_seq = max_proposal_seq,
                                     stable_reps = stable_reps, min_reps = min_reps, multivariable_proposals = multivariable_proposals,
                                     max_reps_keq1 = max_reps_keq1, max_reps_kgt1 = max_reps_kgt1,
                                     delta_lb = delta_lb,
                                     random_seed = random_seed, verbose = verbose, safe_mode = safe_mode);
  assign(penalized_soln_path_name,penalized_soln_path);
  
  # Table of results (Tables 1--3)----
  bpl_unpen = unpenalized_soln$all_params[[1]][,1];
  bpl_unpen[c("log_delta1","log_delta2")] = exp(bpl_unpen[c("log_delta1","log_delta2")]);
  bpl_pen_aic = penalized_soln_path$best_fit_aic[,1];
  bpl_pen_aic[c("log_delta1","log_delta2")] = exp(bpl_pen_aic[c("log_delta1","log_delta2")]);
  bpl_pen_bic = penalized_soln_path$best_fit_bic[,1];
  bpl_pen_bic[c("log_delta1","log_delta2")] = exp(bpl_pen_bic[c("log_delta1","log_delta2")]);
  lambda_values = data.frame(
    problem = "$\\lambda$",
    bpl_unpen = "0",
    bpl_pen_aic = formatC(penalized_soln_path$best_aic["lambda",1],format = "f", digits = 2),
    bpl_pen_bic = formatC(penalized_soln_path$best_bic["lambda",1],format = "f", digits = 2), 
    ldrbo = NA, 
    mean_rank = NA,
    median_rank = "\\\\",
    stringsAsFactors = F)
  
  param_counts = data.frame(
    problem = "$\\tilde p_\\lambda$",
    bpl_unpen = formatC(unpenalized_soln$num_actual_params[1,1],format = "f", digits = 0),
    bpl_pen_aic = formatC(penalized_soln_path$best_aic["num_actual_params",1],format = "f", digits = 0),
    bpl_pen_bic = formatC(penalized_soln_path$best_bic["num_actual_params",1],format = "f", digits = 0), 
    ldrbo = NA, 
    mean_rank = NA,
    median_rank = "\\\\",
    stringsAsFactors = F)
  
  results <- 
    bind_rows(
      data.frame(problem = c("$\\theta_0$",paste0("\\textsc{",tolower(colnames(dat_ranked)),"}"),"$\\delta_1$","$\\delta_2$"), 
                 bpl_unpen = formatC(bpl_unpen,format = "f", digits = 2), 
                 bpl_pen_aic = formatC(bpl_pen_aic,format = "f", digits = 2),
                 bpl_pen_bic = formatC(bpl_pen_bic,format = "f", digits = 2),
                 ldrbo = ldrbo,
                 mean_rank = c(NA, formatC(apply(dat_ranked,2,mean,na.rm=T),format = "f", digits = 1),NA, NA),
                 median_rank = paste0(c("",formatC(apply(dat_ranked,2,median,na.rm=T),format = "f", digits = 1),"", ""),"\\\\"),
                 stringsAsFactors = F) %>% 
        arrange(str_detect(problem, "delta"), desc(bpl_pen_aic),desc(bpl_unpen),desc(ldrbo),mean_rank),
      lambda_values,
      param_counts
    )
  results[nrow(results)-2, ncol(results)] = paste0(results[nrow(results)-1, ncol(results)], "\\hline")
  write_delim(results, path = table_path,delim = "& ",col_names = F, na = "");
  
  # Simulate from model to determine average expected length ----
  n_stages_to_sim = 30;
  n_to_sim = 3e3;
  foo = simulate_bpl(penalized_soln_path$best_fit_aic[,1], stages = 1:n_stages_to_sim, n_to_sim = n_to_sim)
  #mean
  sum((0:(n_stages_to_sim-1)) * colSums(foo == 0, na.rm = T))/n_to_sim
  #median
  max(which(cumsum(colSums(foo == 0, na.rm = T)) < (0.5 * n_to_sim)));
  #1 quartile
  max(which(cumsum(colSums(foo == 0, na.rm = T)) < (0.25 * n_to_sim)));
  #3 quartile
  max(which(cumsum(colSums(foo == 0, na.rm = T)) < (0.75 * n_to_sim)));
  
  summary(rowSums(!is.na(dat_ordered)));
  
  # Solution path ----
  
  foo = penalized_soln_path$all_params[[1]] %>%
    as.data.frame();
  colnames(foo) = as.character(penalized_soln_path$control$lambda_seq);
  foo <- 
    foo %>%
    rownames_to_column() %>%
    gather(key = lambda, value = estimate, -rowname) %>% 
    mutate(rowname = factor(rowname, levels = unique(rowname), ordered = T), 
           lambda = as.numeric(lambda),
           log_lambda = log(lambda)) %>%
    rename(param_name = rowname) %>%
    arrange(lambda,param_name);
  
  foo2 = penalized_soln_path$all_aic %>% 
    as.data.frame();
  colnames(foo2) = as.character(penalized_soln_path$control$lambda_seq);
  foo2 <- 
    foo2 %>% 
    gather(key = lambda, value = aic) %>%
    mutate(lambda = as.numeric(lambda));
  
  foo <- 
    full_join(foo, foo2);
  
  foo2 = penalized_soln_path$all_bic %>% 
    as.data.frame();
  colnames(foo2) = as.character(penalized_soln_path$control$lambda_seq);
  foo2 <- 
    foo2 %>% 
    gather(key = lambda, value = bic) %>%
    mutate(lambda = as.numeric(lambda));
  
  foo <- 
    full_join(foo, foo2);
  
  ggplot() +
    geom_line(data = filter(foo, 
                            log_lambda > log(foo[which.min(foo$bic),"lambda"]) - 1.75, 
                            log_lambda < log(foo[which.min(foo$aic),"lambda"]) + 1.75), 
              aes(x = lambda,
                  y = estimate,
                  group = param_name),
              alpha = 0.3) + 
    geom_vline(data = filter(foo, bic == min(bic), param_name == min(param_name)) %>% filter(lambda == min(lambda)),
               aes(xintercept = lambda, 
                   linetype = "BIC"),
               color = case_color,
               size = 1.25) + 
    geom_vline(data = filter(foo, aic == min(aic), param_name == min(param_name)) %>% filter(lambda == min(lambda)),
               aes(xintercept = lambda, 
                   linetype = "AIC"), 
               color = case_color,
               size = 1.25) + 
    scale_x_continuous(name = expression(lambda), 
                       expand = expand_scale(0), 
                       trans = "log10") +
    scale_y_continuous(name = "Estimate") +
    scale_linetype_discrete(name = "Minimum of") + 
    guides(group = NA) + 
    theme(text = element_text(size = 18),
          legend.position = "top",
          panel.grid.minor = element_blank());
  ggsave(filename = figure_path, height = 5, width = 8);
}

case23_params = case23_penalized_soln_path$best_fit_aic[,1];
foo <- bpl_probs(case23_params, stage = 1)
case23_highest_prob <- max(foo);
names(case23_highest_prob)[length(case23_highest_prob)] <- names(which.max(foo));
while(names(case23_highest_prob)[length(case23_highest_prob)]!="0") {
  case23_params <- case23_params[setdiff(names(case23_params),names(case23_highest_prob))];
  foo <- bpl_probs(case23_params, stage = length(case23_highest_prob) + 1)
  case23_highest_prob <- c(case23_highest_prob,max(foo));
  names(case23_highest_prob)[length(case23_highest_prob)] <- names(which.max(foo));
}
#

case111_params = case111_penalized_soln_path$best_fit_aic[,1];
foo <- bpl_probs(case111_params, stage = 1)
case111_highest_prob <- max(foo);
names(case111_highest_prob)[length(case111_highest_prob)] <- names(which.max(foo));
while(names(case111_highest_prob)[length(case111_highest_prob)]!="0") {
  case111_params <- case111_params[setdiff(names(case111_params),names(case111_highest_prob))];
  foo <- bpl_probs(case111_params, stage = length(case111_highest_prob) + 1)
  case111_highest_prob <- c(case111_highest_prob,max(foo));
  names(case111_highest_prob)[length(case111_highest_prob)] <- names(which.max(foo));
}
#

case83_params = case83_penalized_soln_path$best_fit_aic[,1];
foo <- bpl_probs(case83_params, stage = 1)
case83_highest_prob <- max(foo);
names(case83_highest_prob)[length(case83_highest_prob)] <- names(which.max(foo));
while(names(case83_highest_prob)[length(case83_highest_prob)]!="0") {
  case83_params <- case83_params[setdiff(names(case83_params),names(case83_highest_prob))];
  foo <- bpl_probs(case83_params, stage = length(case83_highest_prob) + 1)
  case83_highest_prob <- c(case83_highest_prob,max(foo));
  names(case83_highest_prob)[length(case83_highest_prob)] <- names(which.max(foo));
}


all_highest_prob <-
  rbind(data.frame(name = "A",
                   order = 1:length(case23_highest_prob),
                   item = c("STOP",colnames(case23_ranked))[1 + as.numeric(names(case23_highest_prob))],
                   prob = case23_highest_prob),
        data.frame(name = "B", 
                   order = 1:length(case111_highest_prob),
                   item = c("STOP",colnames(case111_ranked))[1 + as.numeric(names(case111_highest_prob))],
                   prob = case111_highest_prob),
        data.frame(name = "C", 
                   order = 1:length(case83_highest_prob),
                   item = c("STOP",colnames(case83_ranked))[1 + as.numeric(names(case83_highest_prob))],
                   prob =  case83_highest_prob))


ggplot(data = all_highest_prob,
       aes(x = order,
           y = prob,
           color = name,
           linetype = name,
           shape = name)) + 
  geom_point(size = 2) + 
  geom_path(size = 1.15) + 
  scale_x_continuous(name = "Stage",
                     breaks = 1:10) +
  scale_y_continuous(name = "Highest Probability at Stage",
                     limits = c(0, NA)) +
  scale_color_manual(name = "Case Name", 
                     values = c("#800026","#00441B","#00274C")) + 
  scale_linetype_discrete(name = "Case Name") + 
  scale_shape_discrete(name = "Case Name") + 
  theme(text = element_text(size = 18),
        legend.position = "top",
        legend.key.width = unit(1.5, "cm"));
ggsave(filename = "../../highest_prob_plot.pdf", height = 3.5, width = 8);
  

  
  
