
library(tidyverse);library(knitr);library(kableExtra);
source("functions_bpl.R");
source("functions_ldrbo.R");
source("functions_simulations.R");

# Table 1: Data generating models ----
n_training_seq = c(1);
source("generate_params.R");


all_scenarios_nonnull_label =
  all_scenarios_null_label = 
  expected_length = 
  theta0 = 
  n_items = rep(NA, length(arglist));
for(i in seq_len(length(arglist))) {
  theta0[i] = arglist[[i]]$true_param["0"];
  n_items[i] = length(arglist[[i]]$true_param) - 3;
  expected_length[i] = 
    paste0("\\{", 
           lapply(simulate_bpl(arglist[[i]]$true_param,
                               stages = seq_len(n_items[i]), 
                               n_to_sim = 3e3) %>%
                    is.na() %>%
                    apply(1, which), 
                  min,n_items[i]) %>%
             unlist() %>%
             "-"(1) %>%
             quantile(p = c(0.25, 0.5, 0.75)) %>%
             formatC(format = "f", digits = 0) %>%
             paste0(collapse = ","),
           "\\}");
  
  foo <- arglist[[i]]$true_param[as.character(seq_len(n_items[i]))]
  if(any(foo > 0)) {
    foo2 = foo[foo > 0];
    if(all(near(foo2,foo2[1]))) {
      all_scenarios_nonnull_label[i] = 
        paste0("$\\{",foo2[1],"\\}_{k=1}^{",length(foo2), "}$");
    } else if(all(near(diff(foo2), -0.1))) {
      all_scenarios_nonnull_label[i] = 
        paste0("$\\{",foo2[1],"-0.1(k-1)\\}_{k=1}^{",length(foo2), "}$");
    } else if(all(near(diff(foo2), -0.01))) {
      all_scenarios_nonnull_label[i] = 
        paste0("$\\{",foo2[1],"-0.01(k-1)\\}_{k=1}^{",length(foo2), "}$");
    } else {
      all_scenarios_nonnull_label[i] = 
        paste0("$\\{",paste0(foo[foo > 0], collapse = ","),"\\}$");
    }
    
  } 
  all_scenarios_null_label[i] = sum(near(foo, 0));
}

all_scenarios <-
  bind_cols(label = seq_len(length(arglist)),
            theta0 = theta0,
            nonnull = all_scenarios_nonnull_label,
            number_null = all_scenarios_null_label, 
            n_items = n_items,
            length_q1q2q3 = formatC(expected_length,format="f",digits=1));

null_scenario_labels = which(is.na(all_scenarios$nonnull));

all_scenarios <-
  all_scenarios %>%
  mutate(nonnull = ifelse(!is.na(nonnull), nonnull, "$\\{\\}$"));

all_scenarios_colnames = 
  c("Label",
    "$\\theta_0$",
    "$\\{\\theta_k:\\theta_k > 0\\}$", 
    "$\\#\\{\\theta_k:\\theta_k = 0\\}$", 
    "$\\#\\{\\theta_k\\}$",
    "$\\ell_i:\\{Q_1,Q_2,Q_3\\}$")
linesep_index <- rep("", nrow(all_scenarios));


all_scenarios %>%
  knitr::kable(format = "latex",
               col.names = all_scenarios_colnames,
               booktabs = F,
               longtable = F,
               escape = F,#) %>% 
               linesep = linesep_index) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11);

# Read in results from all method comparisons ----

col_types <- cols(
  sim_id = col_integer(),
  array_id = col_integer(), 
  scenario_id = col_integer(),
  param_id = col_integer(), 
  n_training = col_integer(), 
  n_obs_items = col_integer(),
  n_obs_nonnull_items = col_integer(),
  method_names = col_character(),    
  tau_x= col_double(),
  auc = col_double(),
  rmse_order = col_double(),
  run_time_secs= col_double()
);

raw_all_performance = tibble();
for(i in 1:972) {
  foo <- try(read_csv(paste0("out/sim",i,"_performance.csv"), col_types = col_types));
  if(!"try-error" %in% class(foo)) {
    raw_all_performance = 
      bind_rows(raw_all_performance, foo);
  } else {
    cat("sim ", i, ", not found\n");
  }
}


#number sims
raw_all_performance %>%
  group_by(n_training, param_id, method_names) %>%
  count() %>%
  group_by(n_training, param_id) %>%
  summarize(number_max_sims = max(n),
            number_min_sims = min(n)) %>%
  arrange(n_training, param_id) %>%
  as.data.frame();

# Read in results from BPL comparisons ----

col_types <- cols(
  sim_id = col_integer(),
  array_id = col_integer(),
  scenario_id = col_integer(),
  param_id = col_integer(),
  n_training = col_integer(),
  n_obs_items = col_integer(),
  n_obs_nonnull_items = col_integer(),
  method_names = col_character(),
  neg_loglikelihood = col_double(),
  aic = col_double(),
  lambda = col_double(),
  num_eff_params = col_double(),
  num_actual_params = col_double(),
  max_item_weight = col_double(),
  rmse = col_double(),
  tpr = col_double(),
  tnr = col_double()
);

raw_all_bpl = tibble();
for(i in 1:972) {
  foo <- try(read_csv(paste0("out/sim",i,"_bpl.csv"), col_types = col_types));
  if(!"try-error" %in% class(foo)) {
    raw_all_bpl = 
      bind_rows(raw_all_bpl, foo);
  } else {
    cat("sim ", i, ", not found\n");
  }
}

raw_all_bpl <- 
  left_join(raw_all_bpl, 
            select(raw_all_performance, sim_id, array_id, scenario_id, param_id, n_training, method_names, run_time_secs))

#number sims
raw_all_bpl %>%
  group_by(n_training, param_id, method_names) %>%
  count() %>%
  group_by(n_training, param_id) %>%
  summarize(number_max_sims = max(n),
            number_min_sims = min(n)) %>%
  arrange(n_training, param_id) %>%
  as.data.frame();


## Table 2 (Ordered RMSE) ----

methods_to_include = c("bpl",
                       "bpl_aic",
                       "bpl_bic",
                       "ldrbo",
                       "mc1",
                       "mc2",
                       "mc3",
                       "cemc_s",
                       "cemc_k", 
                       "EMM",
                       "MM");
pretty_names_methods_to_include = 
  c(bpl = "BPL($\\lambda=0$)",
    bpl_aic = "BPL($\\lambda_\\textsc{aic}$)",
    bpl_bic = "BPL($\\lambda_\\textsc{bic}$)",
    ldrbo = "LDRBO",
    mc1 = "MC1",
    mc2 = "MC2",
    mc3 = "MC3",
    cemc_s = "CEMC$_s$",
    cemc_k = "CEMC$_k$",
    EMM = "EMM",
    MM = "MM") %>%
  factor() %>% 
  fct_inorder();

all_performance <- 
  raw_all_performance %>%
  mutate(method_names = 
           factor(method_names, 
                  levels = methods_to_include,
                  labels = pretty_names_methods_to_include[methods_to_include], 
                  ordered = T),
         scenario_id = factor(scenario_id),
         param_id = factor(param_id));


dim(all_performance);

table4_colnames = c("Model", 
                    "$n$", 
                    levels(all_performance$method_names));

table4 <- 
  all_performance %>%
  group_by(param_id, n_training, method_names) %>%
  summarize(mean_auc = mean(1e2 * auc),
            mean_tau_x = mean(1e2 * tau_x),
            mean_rmse_order = mean(1000 * rmse_order), 
            sd_rmse_order = sd(1000 * rmse_order), 
            mean_run_time = formatC(mean(run_time_secs), format = "f", digits = 0)) %>%
  group_by(param_id,  n_training) %>%
  mutate(best1 = ifelse(mean_rmse_order <= 1.05 * min(mean_rmse_order),"\\textbf{", ""),
         best2 = ifelse(mean_rmse_order <= 1.05 * min(mean_rmse_order),"}", "")) %>%
  ungroup() %>%
  arrange(param_id, n_training, method_names) %>%
  mutate(
    mean_auc = formatC(mean_auc, format = "f", digits = 0),
    mean_tau_x= formatC(mean_tau_x, format = "f", digits = 0),
    mean_rmse_order = formatC(mean_rmse_order, format = "f", digits = 0),
    sd_rmse_order = formatC(sd_rmse_order, format = "f", digits = 0),
    combined_metric = paste0(best1, mean_rmse_order, best2)) %>%
  select(param_id, 
         n_training, 
         method_names,
         combined_metric) %>%
  spread(key = method_names, 
         value = combined_metric);

table4_aux <- 
  raw_all_bpl %>% 
  mutate(
    scenario_id = factor(scenario_id),
    param_id = factor(param_id)) %>%
  filter(method_names %in% c("bpl_aic","bpl_bic")) %>% 
  select(n_training, param_id, method_names, rmse_order) %>%
  group_by(n_training, param_id, method_names) %>%
  summarize(bpl_pen_addon = 
              paste0("(",formatC(mean(1000 * rmse_order),
                                 format = "f", 
                                 digits = 0),")")) %>%
  arrange(param_id, n_training)  %>%
  spread(key = method_names, 
         value = bpl_pen_addon);

table4 <-
  left_join(table4, table4_aux) %>%
  mutate(`BPL($\\lambda_\\textsc{aic}$)` = 
           paste0(`BPL($\\lambda_\\textsc{aic}$)`, bpl_aic),
         `BPL($\\lambda_\\textsc{bic}$)` = 
           paste0(`BPL($\\lambda_\\textsc{bic}$)`, bpl_bic)) %>%
  select(-bpl_aic, -bpl_bic);
  

linesep_index <- rep("", nrow(table4));

table4 %>%
  knitr::kable(format = "latex",
               col.names = table4_colnames,
               booktabs = F,
               longtable = F,
               escape = F,#) %>% 
               linesep = linesep_index) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11);

## Table 3 (BPM performance) ----

methods_to_include = c("bpl",
                       "bpl_aic",
                       "bpl_bic");
pretty_names_methods_to_include = 
  c(bpl = "BPL($\\lambda=0$)",
    bpl_aic = "BPL($\\lambda_\\textsc{aic}$)",
    bpl_bic = "BPL($\\lambda_\\textsc{bic}$)") %>%
  factor() %>% 
  fct_inorder();

all_bpl <- 
  raw_all_bpl %>%
  mutate(method_names = 
           factor(method_names, 
                  levels = methods_to_include,
                  labels = pretty_names_methods_to_include[methods_to_include], 
                  ordered = T),
         scenario_id = factor(scenario_id),
         param_id = factor(param_id));



# Table 5: RMSE, running time 

table5_header = 
  c(" " = 2, 
    "RMSE" = 2, 
    "TPR" = 2, 
    "TNR" = 2, 
    "Youden" = 2,
    "Run Time" = 2);

table5_colnames = c("Model", 
                    "$n$", 
                    rep(levels(pretty_names_methods_to_include), times = 5));

table5_ingredients <- 
  all_bpl %>%
  group_by(param_id, n_training, method_names) %>%
  summarize(mean_rmse = mean(1e3 * rmse),
            mean_tpr = mean(1e2 * tpr),
            mean_tnr = mean(1e2 * tnr),
            mean_youden = mean(1e2 * (tpr + tnr - 1)), 
            mean_run_time = mean(run_time_secs)) %>%
  arrange(param_id, n_training, method_names);

table5A <-
  table5_ingredients %>%
  select(param_id, n_training, method_names, mean_rmse) %>%
  group_by(param_id,  n_training) %>%
  mutate(best1 = ifelse(mean_rmse <= 1.05 * min(mean_rmse),"\\textbf{", ""),
         best2 = ifelse(mean_rmse <= 1.05 * min(mean_rmse),"}", "")) %>%
  ungroup() %>%
  mutate(combined_metric = paste0(best1, formatC(mean_rmse,format = "f", digits = 0), best2)) %>%
  select(param_id, n_training, method_names, combined_metric) %>%
  spread(key = method_names, value = combined_metric);

table5B <-
  table5_ingredients %>%
  select(param_id, n_training, method_names, mean_tpr) %>%
  group_by(param_id,  n_training) %>%
  mutate(best1 = ifelse(mean_tpr >= 0.95 * max(mean_tpr),"\\textbf{", ""),
         best2 = ifelse(mean_tpr >= 0.95 * max(mean_tpr),"}", "")) %>%
  ungroup() %>%
  mutate(combined_metric = 
           ifelse(is.nan(mean_tpr), 
                  yes = "-", 
                  no = paste0(best1,  formatC(mean_tpr,format = "f", digits = 0), best2))) %>%
  select(param_id, n_training, method_names, combined_metric) %>%
  spread(key = method_names, value = combined_metric);

table5C <-
  table5_ingredients %>%
  select(param_id, n_training, method_names, mean_tnr) %>%
  group_by(param_id,  n_training) %>%
  mutate(best1 = ifelse(mean_tnr >= 0.95 * max(mean_tnr),"\\textbf{", ""),
         best2 = ifelse(mean_tnr >= 0.95 * max(mean_tnr),"}", "")) %>%
  ungroup() %>%
  mutate(combined_metric = paste0(best1,  formatC(mean_tnr,format = "f", digits = 0), best2)) %>%
  select(param_id, n_training, method_names, combined_metric) %>%
  spread(key = method_names, value = combined_metric);

table5D <-
  table5_ingredients %>%
  select(param_id, n_training, method_names, mean_youden) %>%
  group_by(param_id,  n_training) %>%
  mutate(best1 = ifelse(mean_youden >= 0.95 * max(mean_youden),"\\textbf{", ""),
         best2 = ifelse(mean_youden >= 0.95 * max(mean_youden),"}", "")) %>%
  ungroup() %>%
  mutate(combined_metric = paste0(best1,  formatC(mean_youden,format = "f", digits = 0), best2)) %>%
  select(param_id, n_training, method_names, combined_metric) %>%
  spread(key = method_names, value = combined_metric);

table5E <-
  table5_ingredients %>%
  select(param_id, n_training, method_names, mean_run_time) %>%
  group_by(param_id,  n_training) %>%
  mutate(best1 = ifelse(mean_run_time <= 1.05 * min(mean_run_time),"\\textbf{", ""),
         best2 = ifelse(mean_run_time <= 1.05 * min(mean_run_time),"}", "")) %>%
  ungroup() %>%
  mutate(combined_metric = paste0(best1,  formatC(mean_run_time,format = "f", digits = 0), best2)) %>%
  select(param_id, n_training, method_names, combined_metric) %>%
  spread(key = method_names, value = combined_metric);

table5 <-
  full_join(
    full_join(
      full_join(table5A, table5B, by = c("param_id", "n_training")), 
      full_join(table5C, table5D, by = c("param_id", "n_training")),
      by = c("param_id", "n_training")),
    table5E,
    by = c("param_id", "n_training"));

linesep_index <- rep("", nrow(table5));


table5 %>%
  knitr::kable(format = "latex",
               col.names = table5_colnames,
               booktabs = F,
               longtable = F,
               escape = F,#) %>% 
               linesep = linesep_index) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(table5_header);


