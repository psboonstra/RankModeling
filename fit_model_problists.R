# Run this script to create three .txt files corresponding to Table 4--6 in the
# manuscript formatted as LaTeX tables, three .png files corresponding to Figures
# S1--S3 in the manuscript, and one .png file corresponding to Figure 3 in the 
# manuscript


library(tidyverse);
library(RColorBrewer);
library(TopKLists);
library(ExtMallows);
library(gtools);

source("functions_bpl.R");
source("functions_ldrbo.R");
fig_extension = "png";
source("gather_data.R");

fixed = NULL;
num_inits = 3;
conv_tol = 1e-3;
tau = conv_tol;
num_lambda = 100; 10;#
lambda_min_ratio = 1e-6;
lambda_seq = NULL;
penalty_pow = 1;
proposal_seq = NULL;
proposal_seq_half_length = NULL;
max_proposal_seq = 1;
stable_reps = 3;
min_reps = 3;
multivariable_proposals = T;
delta_lb = 0.25;
max_reps_keq1 = 1e3;
max_reps_kgt1 = 250;
random_seed = sample(2^30.9, 1);
verbose = FALSE;
safe_mode = FALSE;


# Case 23 = Case A
# Case 111 = Case B
# Case 83 = Case C


for(case_name in c("23","111","83")) {
  cat(case_name, "\n\n");
  
  if(case_name == "23") {
    #Case 23, i.e. Case A
    dat_ordered = case23_ordered;
    dat_ranked = case23_ranked;
    table_write_to_path = "case23_all_params.txt";
    figure_write_to_path = paste0("case23_soln_path.",fig_extension);
    
    case_color = "#800026"
    unpenalized_soln_name = "case23_unpenalized_soln";
    penalized_soln_path_name = "case23_penalized_soln_path";
    
    look_beyond_init = 1000;
    look_beyond_final = 100;
    
  } else if(case_name == "111") {
    #Case 111, i.e. Case B
    dat_ordered = case111_ordered;
    dat_ranked = case111_ranked;
    table_write_to_path = "case111_all_params.txt";
    figure_write_to_path = paste0("case111_soln_path.",fig_extension);
    
    case_color = "#00441B"
    unpenalized_soln_name = "case111_unpenalized_soln";
    penalized_soln_path_name = "case111_penalized_soln_path";
    
    look_beyond_init = 500;
    look_beyond_final = 250;
    
  } else if(case_name == "83") {
    #Case 83, i.e. Case C
    dat_ordered = case83_ordered;
    dat_ranked = case83_ranked;
    table_write_to_path = "case83_all_params.txt";
    figure_write_to_path = paste0("case83_soln_path.",fig_extension);
    
    case_color = "#00274C"
    unpenalized_soln_name = "case83_unpenalized_soln";
    penalized_soln_path_name = "case83_penalized_soln_path";
    
    look_beyond_init = 1000;
    look_beyond_final = 100;
  }
  
  ################
  dat_ordered_as_list = matrix_to_list(dat_ordered);
  obs_weights = rep(1, nrow(dat_ordered));
  n_items = length(unique(sort(dat_ordered)));
  
  #Rprof("penRank_run1.out",line.profiling = T)
  unpenalized_soln = penRank_path(dat = dat_ordered, 
                                  obs_weights = obs_weights, 
                                  fixed = fixed, 
                                  num_inits = num_inits, 
                                  conv_tol = conv_tol, 
                                  tau = tau, 
                                  num_lambda = num_lambda, 
                                  lambda_min_ratio = lambda_min_ratio, 
                                  lambda_seq = 0, 
                                  consider_truncating_lambda_seq = FALSE,
                                  penalty_pow = penalty_pow,
                                  proposal_seq = proposal_seq, 
                                  proposal_seq_half_length = proposal_seq_half_length, 
                                  max_proposal_seq = max_proposal_seq,
                                  stable_reps = stable_reps, 
                                  min_reps = min_reps,
                                  multivariable_proposals = multivariable_proposals,
                                  max_reps_keq1 = max_reps_keq1,
                                  max_reps_kgt1 = max_reps_kgt1,
                                  delta_lb = delta_lb,
                                  random_seed = random_seed, 
                                  verbose = verbose, 
                                  safe_mode = safe_mode,
                                  init_params = NULL);
  assign(unpenalized_soln_name,unpenalized_soln);
  
  selected_iter_neg_loglikelihood = which.min(unpenalized_soln$best_neg_loglikelihood["neg_loglikelihood",]);
  
  penalized_soln_path = penRank_path(dat = dat_ordered, 
                                     obs_weights = obs_weights, 
                                     fixed = fixed, 
                                     num_inits = num_inits, 
                                     conv_tol = conv_tol, 
                                     tau = tau, 
                                     num_lambda = num_lambda, 
                                     lambda_min_ratio = lambda_min_ratio, 
                                     lambda_seq = lambda_seq, 
                                     consider_truncating_lambda_seq = TRUE,
                                     penalty_pow = penalty_pow,
                                     proposal_seq = proposal_seq, 
                                     proposal_seq_half_length = proposal_seq_half_length, 
                                     max_proposal_seq = max_proposal_seq,
                                     stable_reps = stable_reps, 
                                     min_reps = min_reps,
                                     multivariable_proposals = multivariable_proposals,
                                     max_reps_keq1 = max_reps_keq1,
                                     max_reps_kgt1 = max_reps_kgt1,
                                     delta_lb = delta_lb,
                                     random_seed = random_seed, 
                                     verbose = verbose, 
                                     safe_mode = safe_mode,
                                     init_params = unpenalized_soln$best_fit_neg_loglikelihood);
  assign(penalized_soln_path_name,penalized_soln_path);
  selected_iter_aic = which.min(penalized_soln_path$best_aic["aic",]);
  selected_iter_bic = which.min(penalized_soln_path$best_bic["bic",]);
  
  
  # Table of results (Tables 4--6)----
  
  # + Bottom part: tuning parameter values ----
  lambda_values = data.frame(
    problem = "$\\lambda$",
    bpl_unpen = "0",
    bpl_pen_aic = formatC(penalized_soln_path$best_aic["lambda",selected_iter_aic],format = "f", digits = 2),
    bpl_pen_bic = formatC(penalized_soln_path$best_bic["lambda",selected_iter_bic],format = "f", digits = 2), 
    ldrbo = NA, 
    mean_rank = NA,
    median_rank = NA, 
    mc1 = NA, 
    mc2 = NA, 
    mc3 = NA, 
    cemc_spearman = NA,
    cemc_kendall = NA, 
    EMM = NA, 
    MM = "\\\\",
    stringsAsFactors = F);
  
  # + Middle part: parameter counts ----
  param_counts = data.frame(
    problem = "$\\tilde p_\\lambda$",
    bpl_unpen = formatC(unpenalized_soln$best_neg_loglikelihood["num_actual_params",selected_iter_neg_loglikelihood],format = "f", digits = 0),
    bpl_pen_aic = formatC(penalized_soln_path$best_aic["num_actual_params",selected_iter_aic],format = "f", digits = 0),
    bpl_pen_bic = formatC(penalized_soln_path$best_bic["num_actual_params",selected_iter_bic],format = "f", digits = 0), 
    ldrbo = NA, 
    mean_rank = NA,
    median_rank = NA, 
    mc1 = NA, 
    mc2 = NA, 
    mc3 = NA, 
    cemc_spearman = NA,
    cemc_kendall = NA, 
    EMM = NA, 
    MM = "\\\\",
    stringsAsFactors = F);
  
  # + Top and most important part: rankings ----
  
  bpl_unpen = unpenalized_soln$best_fit_neg_loglikelihood[,selected_iter_neg_loglikelihood];
  bpl_unpen[c("log_delta1","log_delta2")] = exp(bpl_unpen[c("log_delta1","log_delta2")]);
  bpl_pen_aic = penalized_soln_path$best_fit_aic[,selected_iter_aic];
  bpl_pen_aic[c("log_delta1","log_delta2")] = exp(bpl_pen_aic[c("log_delta1","log_delta2")]);
  bpl_pen_bic = penalized_soln_path$best_fit_bic[,selected_iter_bic];
  bpl_pen_bic[c("log_delta1","log_delta2")] = exp(bpl_pen_bic[c("log_delta1","log_delta2")]);
  
  num_nonzero_aic = sum(!near(bpl_pen_aic[as.character(1:n_items)], 0))
  
  ldrbo_results = consensus_ldrbo(dat = dat_ordered, 
                                  psi = 1,
                                  min_size = 1, 
                                  look_beyond_init = look_beyond_init,
                                  look_beyond_final = look_beyond_final,
                                  verbose = verbose);
  ldrbo_ranks = order(c(ldrbo_results$consensus_list, 
                        setdiff(1:n_items,ldrbo_results$consensus_list)));
  ldrbo_ranks[ldrbo_ranks > length(ldrbo_results$consensus_list)] = "";
  
  mc_results = MC(dat_ordered_as_list);
  mc1A = order(mc_results$MC1.TopK);
  mc1B = paste0("(",trimws(format(100*mc_results$MC1.Prob[mc1A],format = "f", digits = 2)),")");
  
  mc2A = order(mc_results$MC2.TopK);
  mc2B = paste0("(",trimws(format(100*mc_results$MC2.Prob[mc2A],format = "f", digits = 2)),")");
  
  mc3A = order(mc_results$MC3.TopK);
  mc3B = paste0("(",trimws(format(100*mc_results$MC3.Prob[mc3A],format = "f", digits = 2)),")");
  
  cemc_spearman_results = 
    order(CEMC(dat_ordered_as_list, dm = "s")$TopK);
  
  cemc_kendall_results =  
    order(CEMC(dat_ordered_as_list, dm = "k")$TopK);
  
  EMM_results = 
    order(as.numeric(EMM(rankings = t(replace_na(dat_ordered, 0)))$op.pi0));
  
  MM_results = 
    order(as.numeric(MM(rankings = t(replace_na(dat_ordered, 0)))$op.pi0));
  
  rankings <- 
    data.frame(
      problem = c("$\\theta_0$",paste0("\\textsc{",tolower(colnames(dat_ranked)),"}"),"$\\delta_1$","$\\delta_2$"), 
      bpl_unpen = bpl_unpen,
      bpl_pen_aic = bpl_pen_aic,
      bpl_pen_bic = bpl_pen_bic,
      ldrbo = c("", ldrbo_ranks , "", ""),
      mean_rank = c(NA, apply(dat_ranked,2,mean,na.rm=T),NA, NA),
      median_rank = c(NA, apply(dat_ranked,2,median,na.rm=T),NA, NA),
      mc1A = c("", mc1A, "", ""),
      mc1B = c("", mc1B, "", ""),
      mc2A = c("", mc2A, "", ""),
      mc2B = c("", mc2B, "", ""),
      mc3A = c("", mc3A, "", ""),
      mc3B = c("", mc3B, "", ""),
      cemc_spearman = c("", cemc_spearman_results, "", ""),
      cemc_kendall = c("",cemc_kendall_results,"", ""),
      EMM = c("",EMM_results,"", ""),
      MM = c("",MM_results,"", ""),
      eol = "\\\\",
      stringsAsFactors = F) %>% 
    arrange(!str_detect(problem, fixed("theta_0")), 
            str_detect(problem, "delta"),
            desc(bpl_pen_aic),
            desc(bpl_unpen),
            desc(ldrbo),
            mean_rank) %>%
    mutate(bpl_unpen = formatC(bpl_unpen,format = "f", digits = 2), 
           bpl_pen_aic = formatC(bpl_pen_aic,format = "f", digits = 2),
           bpl_pen_bic = formatC(bpl_pen_bic,format = "f", digits = 2),
           rank_counter = cumsum(str_detect(problem, fixed("\\textsc"))),
           ldrbo_pre = replace_na(ifelse(abs(as.numeric(ldrbo) - rank_counter) <= 1 | 
                                           ((as.numeric(ldrbo) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           ldrbo_post = replace_na(ifelse(abs(as.numeric(ldrbo) - rank_counter) <= 1 |
                                            ((as.numeric(ldrbo) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           mean_rank_pretty = formatC(mean_rank,format = "f", digits = 1),
           mean_rank_pre = replace_na(ifelse(abs(mean_rank - rank_counter) <= 1 | 
                                               ((mean_rank > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           mean_rank_post = replace_na(ifelse(abs(mean_rank - rank_counter) <= 1 |
                                                ((mean_rank > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           median_rank_pretty = formatC(median_rank,format = "f", digits = 1),
           median_rank_pre = replace_na(ifelse(abs(median_rank - rank_counter) <= 1 | 
                                                 ((median_rank > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           median_rank_post = replace_na(ifelse(abs(median_rank - rank_counter) <= 1 | 
                                                  ((median_rank > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           mc1_pre = replace_na(ifelse(abs(as.numeric(mc1A) - rank_counter) <= 1 | 
                                         ((as.numeric(mc1A) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           mc1_post = replace_na(ifelse(abs(as.numeric(mc1A) - rank_counter) <= 1 | 
                                          ((as.numeric(mc1A) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           mc2_pre = replace_na(ifelse(abs(as.numeric(mc2A) - rank_counter) <= 1 | 
                                         ((as.numeric(mc2A) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           mc2_post = replace_na(ifelse(abs(as.numeric(mc2A) - rank_counter) <= 1 |
                                          ((as.numeric(mc2A) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           # 
           mc3_pre = replace_na(ifelse(abs(as.numeric(mc3A) - rank_counter) <= 1 | 
                                         ((as.numeric(mc3A) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           mc3_post = replace_na(ifelse(abs(as.numeric(mc3A) - rank_counter) <= 1 | 
                                          ((as.numeric(mc3A) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           cemc_spearman_pre = replace_na(ifelse(abs(as.numeric(cemc_spearman) - rank_counter) <= 1 | 
                                                   ((as.numeric(cemc_spearman) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           cemc_spearman_post = replace_na(ifelse(abs(as.numeric(cemc_spearman) - rank_counter) <= 1 | 
                                                    ((as.numeric(cemc_spearman) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           cemc_kendall_pre = replace_na(ifelse(abs(as.numeric(cemc_kendall) - rank_counter) <= 1 | 
                                                  ((as.numeric(cemc_kendall) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           cemc_kendall_post = replace_na(ifelse(abs(as.numeric(cemc_kendall) - rank_counter) <= 1 | 
                                                   ((as.numeric(cemc_kendall) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           EMM_pre = replace_na(ifelse(abs(as.numeric(EMM) - rank_counter) <= 1 | 
                                         ((as.numeric(EMM) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           EMM_post = replace_na(ifelse(abs(as.numeric(EMM) - rank_counter) <= 1 | 
                                          ((as.numeric(EMM) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),""),
           #
           MM_pre = replace_na(ifelse(abs(as.numeric(MM) - rank_counter) <= 1 | 
                                        ((as.numeric(MM) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "\\textbf{",""),""),
           MM_post = replace_na(ifelse(abs(as.numeric(MM) - rank_counter) <= 1 | 
                                         ((as.numeric(MM) > num_nonzero_aic) & (rank_counter > num_nonzero_aic)), "}",""),"")) %>%
    mutate(ldrbo = paste0(ldrbo_pre, ldrbo, ldrbo_post),
           mean_rank = paste0(mean_rank_pre, mean_rank_pretty, mean_rank_post),
           median_rank = paste0(median_rank_pre, median_rank_pretty, median_rank_post),
           mc1 = paste0(mc1_pre, mc1A,  mc1B, mc1_post),
           mc2 = paste0(mc2_pre, mc2A, mc2B, mc2_post),
           mc3 = paste0(mc3_pre, mc3A, mc3B, mc3_post),
           cemc_spearman = paste0(cemc_spearman_pre, cemc_spearman, cemc_spearman_post),
           cemc_kendall = paste0(cemc_kendall_pre, cemc_kendall, cemc_kendall_post),
           EMM = paste0(EMM_pre, EMM, EMM_post),
           MM = paste0(MM_pre, MM, MM_post, eol)) %>%
    select(problem:median_rank, mc1:mc3, cemc_spearman, cemc_kendall, EMM, MM) %>%
    mutate_all(funs(ifelse(. == "", NA, .))) %>%
    mutate_all(funs(ifelse(. == "NA", NA, .)));
  
  results = 
    bind_rows(
      rankings, 
      lambda_values,
      param_counts
    );
  
  results[nrow(results)-2, ncol(results)] = paste0(results[nrow(results)-1, ncol(results)], "\\hline")
  write_delim(results, path = table_write_to_path,delim = "& ",col_names = F, na = "");
  
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
  ggsave(filename = figure_write_to_path, height = 5, width = 8);
}


# Figure 3 in manuscript ----

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
ggsave(filename = "highest_prob_plot.pdf", height = 3.5, width = 8);




