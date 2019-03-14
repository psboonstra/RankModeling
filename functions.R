# function name and purpose: 'BPL_neg_loglik' calculates the negative-log-likelihood of a BPL-type model when the dampening function
# is as parameterized in Section 3.1.
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 3/5/2019, 10:00am EST
#
# Function inputs:
#
### param: (vector with 'names' attribute, numeric) the parameter vector to evaluate with respect to the data. Must be named exactly 
#as follows: '0' (optional, if the stopping process is to be modeled), then consecutive positive integer labels, one for each 
#unique item, then 'log_delta1', then 'log_delta2'
#
### dat: (data.frame, integers, positive) A data.frame of ordered lists of positive integers. Each row is a different list, each column 
# is a rank, each cell is the integer label of an item. So for example the [1,1] element gives the label of the item that the first 
# ranker ordered first,  i.e. ranked highest. The constraints are as follows: If any cell is non-missing, then every cell in that same
# row to the left should also be non-missing. Second, for any number appearing anywhere in this data.frame, every single number that is 
# less than this number should appear at least once somewhere (not necessarily in the same list)
#
### obs_weights: (vector, integers, positive) A vector with length equal to the number of rows of 'dat' giving the number of observations that
# each list represents. This is gives a simple way to considerably speed up the calculation if there are a large number of duplicated lists. Just
# provide a vector of 1's if each list represents a single ranker
#
### fixed: (vector with 'names' attribute, numeric) a numeric vector containing a subset of elements in 'param' the values of which 
# will be over-written by the values in 'fixed'. This is not useful in isolation but helpful for the coordinate ascent algorithm
#
### safe_mode: (logical) if TRUE, then perform error checks to make sure that all inputs are as expected. Set to FALSE to speed things up a bit
#
### num_uniq: (integer, positive) number of unique items in 'dat'. Only required if safe_mode == FALSE, otherwise will be equal to sum(!is.na(unique(as.numeric(dat))))
#
### num_obs: (integer, positive) number of observations in 'dat'. Only required if safe_mode == FALSE, otherwise will be equal to nrow(dat)
#
### num_ranks: (integer, positive) 1 + length of longest observed list. Only required if safe_mode == FALSE, otherwise will be equal to 1 + ncol(dat)
#
### max_columns: (integer, positive) number of columns need to loop over. Only required if safe_mode == FALSE, otherwise will be equal to  min(num_uniq, num_ranks)
#
### list_lengths: (vector, positive) vector of individual list lengths. Only required if safe_mode == FALSE, otherwise will be equal to unlist(lapply(apply(apply(dat,1,is.na),2,which),"[",1)) 
#
### Value
#The weighted sum of negative of the BPL log-likelihood function (Equation 3.2) over the n observations, weighted by obs_weights

BPL_neg_loglik = function(param, 
                          dat, 
                          obs_weights, 
                          fixed = NULL,
                          safe_mode = TRUE, 
                          num_uniq = NULL, 
                          num_obs = NULL, 
                          num_ranks = NULL, 
                          max_columns = NULL,
                          list_lengths = NULL) {
  
  #augment the data with a column of NAs so that there is always room to add a zero after the end of an actual list
  aug_dat = cbind(dat,NA);
  #conduct the checks that all inputs are as expected
  if(safe_mode) {
    num_obs = nrow(aug_dat);
    num_ranks = ncol(aug_dat);
    num_uniq = sum(!is.na(unique(as.numeric(dat))));
    max_columns = min(num_uniq, num_ranks);
    list_lengths = unlist(lapply(apply(apply(aug_dat,1,is.na),2,which),"[",1)) - 1;
    if(!isTRUE(all(apply(apply(is.na(aug_dat),1,cumsum),2,diff) | (apply(is.na(dat),1,cumsum) == 0))) | 
       isTRUE(any(apply(dat,1,duplicated) & !apply(dat,1,is.na)))) {
      stop("'dat' must contain uninterrupted left-to-right sets of unique integers, i.e. no NAs between items");
    }
    if(!isTRUE(all.equal(names(param),c(0:num_uniq,"log_delta1","log_delta2"))) && !isTRUE(all.equal(names(param),c(1:num_uniq,"log_delta1","log_delta2")))) {
      stop("'param' must be named exactly as follows: '0' (optional, if the stopping process is to be modeled), then consecutive positive integer labels, one for each unique item, then 'log_delta1', then 'log_delta2'");
    }
    if(!all(names(fixed)%in%names(param))) {stop("'fixed' must be NULL or a subset of named elements from param");}
    if(num_uniq != max(dat,na.rm = TRUE)) {stop("Item labels in 'dat' must be consecutively labeled integers starting with 1");}
    if(num_uniq != length(param) - 2 - ("0"%in%names(param))) {stop("Length of 'param' must be equal to the number of unique items ranked plus 2 (for 'log_delta1' and 'log_delta2') plus 1 (if '0' is modeled");}
    if("log_delta1"%in%names(param) && param["log_delta1"] > 0) {stop("'log_delta1' must be non-positive");}
    if("log_delta2"%in%names(param) && param["log_delta2"] > 0) {stop("'log_delta2' must be non-positive");}
  }
  #Overwrite with fixed values
  param[names(fixed)] = fixed;
  log_item_weights = param[as.character(1:num_uniq)];
  item_weights = exp(log_item_weights);
  delta1 = exp(param["log_delta1"]);
  delta2 = exp(param["log_delta2"]);
  #Insert zero after end of every list to indicate 'stop'
  aug_dat[cbind(1:num_obs,list_lengths + 1)] = 0;
  #Check if fatigue is to be monitored
  if("0"%in%names(param)) {
    log_fatigue_weight = param["0"];
    fatigue_weight = exp(log_fatigue_weight);
  } else {
    log_fatigue_weight = fatigue_weight = 0;
  }
  #Contribution from Pr(xis=k|ui,Ois,xis!=0)
  log_num = log_denom = matrix(0,nrow=num_obs,ncol=num_ranks);
  log_num[,1] = log_item_weights[aug_dat[,1]];
  #Denominator is same for everyone at rank 1 (and does not include the constant '1' since no list can stop at rank 1)
  log_denom[,1] = log(sum(item_weights));
  index_to_loop = 1:max_columns;
  #all dampening effects at once
  dampening = (delta2 * delta1^(index_to_loop-1) + (1-delta2)^(2 * index_to_loop - 1))
  for(s in index_to_loop[-1]) {
    #which lists ranked an item at this stage?
    which_ranked_item = which(aug_dat[,s] > 0);
    #which lists ended precisely at this stage?
    which_ended = which(aug_dat[,s] == 0);
    if(length(which_ranked_item)) {
      log_num[which_ranked_item,s] = log_item_weights[aug_dat[which_ranked_item,s]] * dampening[s];
    }
    if(length(which_ended)) {
      log_num[which_ended,s] = log_fatigue_weight;
    }
    which_either = which(aug_dat[,s] >= 0);
    log_denom[which_either,s] = log(fatigue_weight + sum(item_weights^dampening[s]) - rowSums(matrix(item_weights[aug_dat[which_either,1:(s-1),drop = FALSE]]^dampening[s], ncol = s-1),na.rm = TRUE));
  }  
  sum(obs_weights * log_denom) - sum(obs_weights * log_num);
}

# function name and purpose: 'BPL_neg_loglik_CA' calculates a grid of the negative-log-likelihood value sof a BPL-type model 
# when the dampening function is as parameterized in Section 3.1. However, one element of 'param' can be missing (NA), and 
# in this case the function calculates and returns multiple values of the negative log-likelihood by plugging in 
# each element of 'CA_param' into the missing value 
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 3/5/2019, 10:00am EST
#
# Function inputs:
#
### param: (vector with 'names' attribute, numeric) the parameter vector to evaluate with respect to the data. Must be named exactly 
#as follows: '0' (optional, if the stopping process is to be modeled), then consecutive positive integer labels, one for each 
#unique item, then 'log_delta1', then 'log_delta2'. Up to one element of 'param' can be missing (NA), in which case CA_param
# should be non-missing
#
### CA_param (vector, numeric) a grid of values to plug in for the missing element of 'param'. The function will return a vector
# of negative log-likelihoods as long as this grid
#
### dat: (data.frame, integers, positive) A data.frame of ordered lists of positive integers. Each row is a different list, each column 
# is a rank, each cell is the integer label of an item. So for example the [1,1] element gives the label of the item that the first 
# ranker ordered first,  i.e. ranked highest. The constraints are as follows: If any cell is non-missing, then every cell in that same
# row to the left should also be non-missing. Second, for any number appearing anywhere in this data.frame, every single number that is 
# less than this number should appear at least once somewhere (not necessarily in the same list)
#
### obs_weights: (vector, integers, positive) A vector with length equal to the number of rows of 'dat' giving the number of observations that
# each list represents. This is gives a simple way to considerably speed up the calculation if there are a large number of duplicated lists. Just
# provide a vector of 1's if each list represents a single ranker
#
### safe_mode: (logical) if TRUE, then perform error checks to make sure that all inputs are as expected. Set to FALSE to speed things up a bit
#
### num_uniq: (integer, positive) number of unique items in 'dat'. Only required if safe_mode == FALSE, otherwise will be equal to sum(!is.na(unique(as.numeric(dat))))
#
### num_obs: (integer, positive) number of observations in 'dat'. Only required if safe_mode == FALSE, otherwise will be equal to nrow(dat)
#
### num_ranks: (integer, positive) 1 + length of longest observed list. Only required if safe_mode == FALSE, otherwise will be equal to 1 + ncol(dat)
#
### max_columns: (integer, positive) number of columns need to loop over. Only required if safe_mode == FALSE, otherwise will be equal to  min(num_uniq, num_ranks)
#
### list_lengths: (vector, positive) vector of individual list lengths. Only required if safe_mode == FALSE, otherwise will be equal to unlist(lapply(apply(apply(dat,1,is.na),2,which),"[",1)) 
#
### CA_param_counts_by_rank: (vector, integer, positive) vector of counts that the missing element of 'param' is ranked at each stage. Only required if 
# 'param' contains a missing element, that missing element is an item, and safe_mode == FALSE. Otherwise, will be equal to 
# c(colSums(obs_weights * (dat == as.numeric(CA_flag)), na.rm = TRUE),0)[1:max_columns]
#
### Value
#A vector as long as 'CA_param', with each element giving the weighted sum of negative of the BPL log-likelihood function (Equation 3.2) over the n observations, 
# weighted by obs_weights, after plugging in each element of 'CA_param' into the missing value of 'param'. If 'param' has no missing elements, 
# then this just returns the negative BPL log-likelihood


BPL_neg_loglik_CA = function(param,
                             CA_param,
                             dat,
                             obs_weights, 
                             safe_mode = TRUE, 
                             num_uniq = NULL, 
                             num_obs = NULL, 
                             num_ranks = NULL, 
                             max_columns = NULL,
                             list_lengths = NULL, 
                             CA_param_counts_by_rank = NULL) {
  
  #name of parameter to grid-search over
  CA_flag = names(param)[which(is.na(param))];
  #length of grid
  length_CA_param = length(CA_param);
  #augment the data with a column of NAs so that there is always room to add a zero after the end of an actual list
  aug_dat = cbind(dat,NA);
  #conduct the checks that all inputs are as expected
  if(safe_mode) {
    num_obs = nrow(aug_dat);
    num_ranks = ncol(aug_dat);
    num_uniq = sum(!is.na(unique(as.numeric(dat))));
    max_columns = min(num_uniq, num_ranks);
    list_lengths = unlist(lapply(apply(apply(aug_dat,1,is.na),2,which),"[",1)) - 1;
    num_uniq = sum(!is.na(unique(as.numeric(dat))));
    if(length(CA_flag) > 1) {stop("coordinate ascent cannot maximize with respect to more than one parameter");}
    if(!isTRUE(all.equal(names(param),c(0:num_uniq,"log_delta1","log_delta2"))) && !isTRUE(all.equal(names(param),c(1:num_uniq,"log_delta1","log_delta2")))) {
      stop("'param' must be named exactly as follows: '0' (optional, if the stopping process is to be modeled), then consecutive positive integer labels, one for each unique item, then 'log_delta1', then 'log_delta2'");
    }
    if(num_uniq != max(dat,na.rm = TRUE)) {stop("Item labels in 'dat' must be consecutively labeled integers starting with 1");}
    if(num_uniq!=length(param) - 2 - ("0"%in%names(param))) {stop("Length of 'param' must be equal to the number of unique items ranked plus 2 (for 'log_delta1' and 'log_delta2') plus 1 (if '0' is modeled");}
    if(length(CA_flag) && !(CA_flag %in% c("0","log_delta1","log_delta2"))) {
      CA_param_counts_by_rank = c(colSums(obs_weights * (dat == as.numeric(CA_flag)), na.rm = TRUE),0)[1:max_columns];#Last rank indicates termination, so no items are ever ranked there
      #unlist(lapply(apply(is.na(log_num),2,which),"length"))
    }
    if("log_delta1"%in%names(param) && CA_flag != "log_delta1" && param["log_delta1"] > 0) {stop("'log_delta1' must be non-positive");}
    if("log_delta2"%in%names(param) && CA_flag != "log_delta2" && param["log_delta2"] > 0) {stop("'log_delta2' must be non-positive");}
    if(length(CA_flag) && CA_flag == "log_delta1" && any(CA_param > 0)) {stop("'log_delta1' must be non-positive");}
    if(length(CA_flag) && CA_flag == "log_delta2" && any(CA_param > 0)) {stop("'log_delta2' must be non-positive");}
  }
  log_item_weights = param[as.character(1:num_uniq)];
  item_weights = exp(log_item_weights);
  delta1 = exp(param["log_delta1"]);
  delta2 = exp(param["log_delta2"]);
  #Insert zero after end of every list to indicate 'stop'
  aug_dat[cbind(1:num_obs,list_lengths + 1)] = 0;
  #Check if fatigue is to be monitored
  if("0"%in%names(param)) {
    log_fatigue_weight = param["0"];
    fatigue_weight = exp(log_fatigue_weight);
  } else {
    log_fatigue_weight = fatigue_weight = 0;
  }
  #Do this if a grid of dampening parameters is provided
  if(length(CA_flag) && (CA_flag %in% c("log_delta1","log_delta2"))) {
    log_num = 
      log_denom = array(0, dim = c(num_obs,num_ranks,length_CA_param));
    log_num[,1,] = log_item_weights[aug_dat[,1]];
    #Denominator is same for everyone at rank 1 (and does not include the constant '1' since no list can stop at rank 1)
    log_denom[,1,] = log(sum(item_weights));
    index_to_loop = 1:max_columns;
    expanded_index_to_loop = matrix(index_to_loop,nrow = length_CA_param, ncol = max_columns, byrow = TRUE);
    if(CA_flag %in% "log_delta1") {
      dampening_matrix = delta2 * matrix(exp(CA_param),nrow = length_CA_param, ncol = max_columns)^(expanded_index_to_loop - 1) + (1-delta2)^(2 * expanded_index_to_loop - 1);
    } else {
      matrix_delta2 = matrix(exp(CA_param), nrow = length_CA_param, ncol = max_columns);
      dampening_matrix = matrix_delta2 * delta1^(expanded_index_to_loop - 1) + (1 - matrix_delta2)^(2 * expanded_index_to_loop - 1);
    }
    dampened_item_weights = array(rep(item_weights, each = prod(dim(dampening_matrix))), dim = c(length_CA_param,num_ranks,num_uniq)) ^ array(dampening_matrix,dim = c(length_CA_param,num_ranks,num_uniq));
    sum_dampened_item_weights = rowSums(dampened_item_weights,dims = 2);
    for(s in index_to_loop[-1]) {
      which_ranked_item = which(aug_dat[,s] > 0);
      which_ended = which(aug_dat[,s] == 0);
      if(length(which_ranked_item)) {
        log_num[which_ranked_item,s,] = tcrossprod(log_item_weights[aug_dat[which_ranked_item,s]],dampening_matrix[,s]);
      }
      if(length(which_ended)) {
        log_num[which_ended,s,] = log_fatigue_weight;
      }
      which_either = which(aug_dat[,s] >= 0);
      num_either = sum(aug_dat[,s] >= 0, na.rm = TRUE);
      dampened_weight_loss = array(dampened_item_weights[,s,][cbind(rep(1:length_CA_param, each = num_either*(s-1)),  rep(as.numeric(aug_dat[which_either,1:(s-1),drop = FALSE]), times = length_CA_param))], dim = c(num_either, s - 1, length_CA_param));
      log_denom[which_either,s,] = log(fatigue_weight +
                                         matrix(sum_dampened_item_weights[,s], nrow = num_either, ncol = length_CA_param, byrow = TRUE) - 
                                         rowSums(aperm(dampened_weight_loss,c(1,3,2)), dims = 2))
    }  
    colSums(obs_weights * log_denom - obs_weights * log_num, dims = 2);
    #Do this if a grid of fatigue parameters is provided
  } else if(length(CA_flag) && (CA_flag == "0")) {
    log_num = 
      log_denom = array(0, dim = c(num_obs,num_ranks,length_CA_param));
    log_num[,1,] = log_item_weights[aug_dat[,1]];
    #Denominator is same for everyone at rank 1 (and does not include the constant '1' since no list can stop at rank 1)
    log_denom[,1,] = log(sum(item_weights));
    index_to_loop = 1:max_columns;
    dampening = (delta2 * delta1^(index_to_loop-1) + (1-delta2)^(2 * index_to_loop - 1));
    for(s in index_to_loop[-1]) {
      which_ranked_item = which(aug_dat[,s] > 0);
      which_ended = which(aug_dat[,s] == 0);
      if(length(which_ranked_item)) {
        log_num[which_ranked_item,s,] = log_item_weights[aug_dat[which_ranked_item,s]] * dampening[s];
      }
      if(length(which_ended)) {
        log_num[which_ended,s,] = rep(CA_param, each = length(which_ended));
      }
      which_either = which(aug_dat[,s] >= 0);
      num_either = sum(aug_dat[,s] >= 0, na.rm = TRUE);
      log_denom[which_either,s,] = log(rep(exp(CA_param), each = num_either) + rep(sum(item_weights^dampening[s]) - rowSums(matrix(item_weights[aug_dat[which_either,1:(s-1),drop = FALSE]]^dampening[s],ncol=s-1),na.rm = TRUE), times = length_CA_param));
    }  
    colSums(obs_weights * log_denom - obs_weights * log_num, dims = 2);
    #Do this if a grid of item weights is provided, or if 'param' is complete and no grid is provided
  } else {
    ##Insert zero after end of every list to indicate 'stop'
    index_to_loop = 1:max_columns;
    log_num = 
      denom_flag = matrix(0,nrow=num_obs,ncol=max_columns);
    denom = matrix(1,nrow=num_obs,ncol=max_columns);
    log_num[,1] = log_item_weights[aug_dat[,1]];
    denom[,1] = sum(item_weights,na.rm = TRUE);
    denom_flag[,1] = 1;
    dampening = (delta2 * delta1^(index_to_loop-1) + (1-delta2)^(2 * index_to_loop - 1));
    as_numeric_CA_flag = as.numeric(CA_flag);
    for(s in index_to_loop[-1]) {
      which_ranked_item = which(aug_dat[,s] > 0);
      which_ended = which(aug_dat[,s] == 0);
      if(length(which_ranked_item)) {
        log_num[which_ranked_item,s] = log_item_weights[aug_dat[which_ranked_item,s]] * dampening[s];
      }
      if(length(which_ended)) {
        log_num[which_ended,s] = log_fatigue_weight;
      }
      which_either = which(aug_dat[,s] >= 0);
      aug_dat_reduced = aug_dat[which_either,1:(s-1),drop = FALSE];
      dampened_item_weights = item_weights ^ dampening[s];
      denom[which_either,s] = 
        fatigue_weight + 
        sum(dampened_item_weights, na.rm = TRUE) - 
        rowSums(matrix(dampened_item_weights[aug_dat_reduced], ncol = s-1), na.rm = TRUE);
      if(length(CA_flag)) {
        denom_flag[which_either,s] = 1 - rowSums(aug_dat_reduced == as_numeric_CA_flag);
      }
    }  
    dampening_denom_flag = denom_flag * matrix(dampening, nrow = num_obs, ncol = max_columns, byrow = TRUE) + (1 - denom_flag);
    #delta1_denom_flag = delta1^(denom_flag * matrix(0:(ncol(denom_flag)-1),nrow=nrow(denom_flag),ncol=ncol(denom_flag),byrow = TRUE));
    if(length(CA_flag)) {
      rowSums(matrix(obs_weights, nrow = length_CA_param, ncol = prod(dim(denom)), byrow = TRUE) * log(tcrossprod(rep(1,length_CA_param),as.numeric(denom)) + tcrossprod(exp(CA_param),as.numeric(denom_flag))^tcrossprod(rep(1,length_CA_param),as.numeric(dampening_denom_flag)))) - 
        sum(obs_weights * log_num,na.rm = TRUE) - 
        sum(CA_param_counts_by_rank * dampening) * CA_param;
    } else {
      sum(obs_weights * log(denom)) - sum(obs_weights * log_num);
    }
  }
}

# function name and purpose: 'over_identified_cost' is a self-contained check to see if the parameter violates any of 
# the constraints needed to identify a BPL-type model. It is designed to accompany 'BPL_neg_loglik_CA', thus the user may provide
# a vector 'CA_param' to separately determine whether any of the elements of 'CA_param' result in a constraint violation. 
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 3/5/2019, 10:00am EST
#
# Function inputs:
#
### param: (vector with 'names' attribute, numeric) the parameter vector to evaluate with respect to the data. Must be named exactly 
#as follows: '0' (optional, if the stopping process is to be modeled), then consecutive positive integer labels, one for each 
#unique item, then 'log_delta1', then 'log_delta2'. Up to one element of 'param' can be missing (NA), in which case CA_param
# should be non-missing
#
### CA_param (vector, numeric) a grid of values to plug in for the missing element of 'param'. The function will return a vector
# of negative log-likelihoods as long as this grid
#
### safe_mode: (logical) if TRUE, then perform error checks to make sure that all inputs are as expected. Set to FALSE to speed things up a bit
#
### tiny_positive: (numeric, positive) self-explanatory
#
### Value
# A vector as long as 'CA_param', if it is provided and non-missing, otherwise a scalar if 'param' is fully observed. Each element will 
# be 0 if the function determines that the BPL is identified under this configuration, or Inf otherwise. 

over_identified_cost = function(param,
                                CA_param, 
                                safe_mode = TRUE,
                                tiny_positive = .Machine$double.eps^0.5) {
  #name of parameter to grid-search over
  CA_flag = names(param)[which(is.na(param))];
  item_names = intersect(1:length(param), names(param));
  if(safe_mode) {
    if(tiny_positive <= 0) {warning("In 'over_identified_cost', 'tiny_positive' seems not to be a tiny positive...");}
    if("log_delta1"%in%names(param) && CA_flag != "log_delta1" && param["log_delta1"] > 0) {stop("'log_delta1' must be non-positive");}
    if("log_delta2"%in%names(param) && CA_flag != "log_delta2" && param["log_delta2"] > 0) {stop("'log_delta2' must be non-positive");}
    if(length(CA_flag) && CA_flag == "log_delta1" && any(CA_param > 0)) {stop("'log_delta1' must be non-positive");}
    if(length(CA_flag) && CA_flag == "log_delta2" && any(CA_param > 0)) {stop("'log_delta2' must be non-positive");}
  }
  #This block handles the case that CA_flag is an actual item
  if(length(CA_flag) && !(CA_flag %in% c("0","log_delta1","log_delta2"))) {
    theta_min = pmin(min(param[item_names],na.rm = TRUE), CA_param);
  } else if(length(CA_flag) && (CA_flag %in% c("0","log_delta1","log_delta2"))) {
    theta_min = rep(min(param[item_names]),length(CA_param));
  } else {
    #Otherwise we assume param is completely provided
    theta_min = min(param[item_names]);
  }
  (abs(theta_min) > tiny_positive | theta_min < 0) / tiny_positive;
}

# function name and purpose: 'seamless_L0' counts the effective number of non-zero parameters in a given BPL-type model
# for a given value of tau. 
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 3/5/2019, 10:00am EST
#
# Function inputs:
#
### param: (vector with 'names' attribute, numeric) the parameter vector to evaluate with respect to the data. Must be named exactly 
#as follows: '0' (optional, if the stopping process is to be modeled), then consecutive positive integer labels, one for each 
#unique item, then 'log_delta1', then 'log_delta2'. Up to one element of 'param' can be missing (NA), in which case CA_param
# should be non-missing
#
### CA_param (vector, numeric) a grid of values to plug in for the missing element of 'param'. The function will return a vector
# of negative log-likelihoods as long as this grid
#
### tau (numeric, non-negative) When positive, this makes the count smoothly varying between 0 and 1. 
#
### penalty_pow (numeric, positive) defaults to 1. This was not changed in the paper, but it is the power to which
# both 'tau' and each parameter are raised, which tends to make the shape of the seamless penalty more ridge-like close to zero. 
#
### safe_mode: (logical) if TRUE, then perform error checks to make sure that all inputs are as expected. Set to FALSE to speed things up a bit
#
### num_uniq: (integer, positive) number of unique items in 'dat'. Only required if safe_mode == FALSE, otherwise will be equal to sum(!is.na(unique(as.numeric(dat))))
#
### Value
# A vector as long as 'CA_param', if it is provided and non-missing, otherwise a scalar if 'param' is fully observed. Each element will 
# be 0 if the function determines that the BPL is identified under this configuration, or Inf otherwise. 

seamless_L0 = function(param,
                       CA_param,
                       tau, 
                       penalty_pow = 1,
                       safe_mode = TRUE,
                       num_uniq = NULL) {
  #name of parameter to grid-search over
  CA_flag = names(param)[which(is.na(param))];
  item_names = intersect(1:length(param), names(param));
  fatigue_param_exists = ("0"%in%names(param));
  if(safe_mode) {
    num_uniq = length(item_names);
    if("log_delta1"%in%names(param) && CA_flag != "log_delta1" && param["log_delta1"] > 0) {stop("'log_delta1' must be non-positive");}
    if("log_delta2"%in%names(param) && CA_flag != "log_delta2" && param["log_delta2"] > 0) {stop("'log_delta2' must be non-positive");}
    if(length(CA_flag) && CA_flag == "log_delta1" && any(CA_param > 0)) {stop("'log_delta1' must be non-positive");}
    if(length(CA_flag) && CA_flag == "log_delta2" && any(CA_param > 0)) {stop("'log_delta2' must be non-positive");}
    if(penalty_pow <= 0) {stop("'penalty_pow' must be strictly positive");}
    if(tau < 0) {stop("'tau' must be strictly positive");}
  }
  #This block handles the case that CA_flag is an actual item
  if(length(CA_flag) && (length_CA_param <- length(CA_param)) > 1) {
    if(!(CA_flag %in% c("0","log_delta1","log_delta2"))) {
      theta_min = pmin(min(param[item_names],na.rm = TRUE), CA_param);
      delta1 = exp(param["log_delta1"]);
      delta2 = exp(param["log_delta2"]);
      #foo is a three-dimensional array: 
      #first dimension (rows) has length equal to length of CA_param
      #second dimension (columns) has length equal to length of item_names
      #third dimension also has length equal to number of items to rank and corresponds to rank
      foo = array(rep(sapply(param[item_names],"-",theta_min), each = num_uniq), 
                  dim = c(length_CA_param, num_uniq, num_uniq),
                  dimnames = list(NULL,item_names,item_names));
      foo[,,CA_flag] = CA_param - theta_min; 
      foo3 = array(rep(0:(num_uniq-1),each = length_CA_param) ,dim = dim(foo));
      foo4 = foo * (delta2 * (delta1 ^ foo3) + (1 - delta2)^(2 * foo3 + 1));
      #Final count:
      result = fatigue_param_exists +
        log2(1/(1 + abs(tau/param["log_delta1"])^penalty_pow) + 1) + 
        log2(1/(1 + abs(tau/param["log_delta2"])^penalty_pow) + 1) + 
        rowMeans(rowSums(log2(1/(1+(tau/foo4)^penalty_pow) + 1), dims = 2));
      #This block handles the case that CA_flag is log_delta1
    } else if(CA_flag == "log_delta1") {
      theta_min = rep(min(param[item_names]), length(CA_param));
      delta1 = exp(CA_param);
      delta2 = exp(param["log_delta2"]);
      foo = array(rep(sapply(param[item_names], "-", theta_min), each = num_uniq), 
                  dim = c(length_CA_param, num_uniq, num_uniq),
                  dimnames = list(NULL, item_names, item_names));
      #foo is a 3-dimensional array: 
      #first dimension has length equal to length of CA_param
      #second dimension has length equal to number of items to rank and corresponds to rank (this is the dimension that CA is exploring)
      #third dimension also has length equal to number of items to rank 
      foo2 = matrix(0:(num_uniq-1),nrow = length_CA_param, ncol = num_uniq, byrow = TRUE);
      foo3 = matrix(delta1,nrow = length_CA_param, ncol = num_uniq);
      foo4 = foo * array(delta2 * foo3 ^ foo2 + (1 - delta2)^(2 * foo2 + 1), dim = dim(foo));
      #Final count:
      result = fatigue_param_exists + 
        log2(1/(1 + abs(tau/CA_param)^penalty_pow) + 1) + 
        log2(1/(1 + abs(tau/param["log_delta2"])^penalty_pow) + 1) + 
        rowSums(rowMeans(log2(1/(1+(tau/foo4)^penalty_pow)+1), dims = 2));
      #This block handles the case that CA_flag is log_delta2
    } else if(CA_flag == "log_delta2") {
      theta_min = rep(min(param[item_names]), length(CA_param));
      delta1 = exp(param["log_delta1"]);
      delta2 = exp(CA_param);
      foo = array(rep(sapply(param[item_names], "-", theta_min), each = num_uniq), 
                  dim = c(length_CA_param, num_uniq, num_uniq),
                  dimnames = list(NULL, item_names, item_names));
      #foo is a 3-dimensional array: 
      #first dimension has length equal to length of CA_param
      #second dimension has length equal to number of items to rank and corresponds to rank (this is the dimension that CA is exploring)
      #third dimension also has length equal to number of items to rank 
      foo2 = matrix(0:(num_uniq-1),nrow = length_CA_param, ncol = num_uniq, byrow = TRUE);
      foo3 = matrix(delta2,nrow = length_CA_param, ncol = num_uniq);
      foo4 = foo * array(foo3 * delta1 ^ foo2 + (1 - foo3)^(2 * foo2 + 1), dim = dim(foo));
      #Final count:
      result = fatigue_param_exists + 
        log2(1/(1 + abs(tau/param["log_delta1"])^penalty_pow) + 1) + 
        log2(1/(1 + abs(tau/CA_param)^penalty_pow) + 1) + 
        rowSums(rowMeans(log2(1/(1+(tau/foo4)^penalty_pow)+1),dims=2));
    } else if(CA_flag == "0") {
      delta1 = exp(param["log_delta1"]);
      delta2 = exp(param["log_delta2"]);
      #
      foo = matrix(param[item_names] - min(param[item_names]), nrow = num_uniq, ncol = num_uniq, byrow = TRUE);
      foo2 = 0:(num_uniq-1);
      foo3 = foo * matrix((delta2 * delta1 ^ (foo2) + (1 - delta2) ^ (2 * foo2 + 1)), nrow = num_uniq, ncol = num_uniq);
      #Final count:
      result = fatigue_param_exists + 
        log2(1/(1 + abs(tau/param["log_delta1"])^penalty_pow) + 1) + 
        log2(1/(1 + abs(tau/param["log_delta2"])^penalty_pow) + 1) + 
        mean(rowSums(log2(1/(1+(tau/foo3)^penalty_pow) + 1)));
    } 
  } else {
    #Otherwise we assume param is completely provided
    if(length(CA_flag)) {
      param[CA_flag] = CA_param;
    }
    theta_min = min(param[item_names]);
    delta1 = exp(param["log_delta1"]);
    delta2 = exp(param["log_delta2"]);
    #
    foo = matrix(param[item_names] - min(param[item_names]), nrow = num_uniq, ncol = num_uniq, byrow = TRUE);
    foo2 = 0:(num_uniq-1);
    foo3 = foo * matrix((delta2 * delta1 ^ (foo2) + (1 - delta2) ^ (2 * foo2 + 1)), nrow = num_uniq, ncol = num_uniq);
    result = fatigue_param_exists +
      log2(1/(1 + abs(tau/param["log_delta1"])^penalty_pow) + 1) + 
      log2(1/(1 + abs(tau/param["log_delta2"])^penalty_pow) + 1) + 
      mean(rowSums(log2(1/(1+(tau/foo3)^penalty_pow) + 1)));
  }
  result;
} 

# function name and purpose: 'penRank_path' calculates the full solution path of a penalized BPL model as well as estimating
# the optimal penalization as given by that minimizing the small-sample AIC and the BIC.
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 3/5/2019, 10:00am EST
#
# Function inputs:
#
### dat: (data.frame, integers, positive) A data.frame of ordered lists of positive integers. Each row is a different list, each column 
# is a rank, each cell is the integer label of an item. So for example the [1,1] element gives the label of the item that the first 
# ranker ordered first,  i.e. ranked highest. The constraints are as follows: If any cell is non-missing, then every cell in that same
# row to the left should also be non-missing. Second, for any number appearing anywhere in this data.frame, every single number that is 
# less than this number should appear at least once somewhere (not necessarily in the same list)
#
### obs_weights: (vector, integers, positive) A vector with length equal to the number of rows of 'dat' giving the number of observations that
# each list represents. This is gives a simple way to considerably speed up the calculation if there are a large number of duplicated lists. Just
# provide a vector of 1's if each list represents a single ranker
#
### fixed: (either NULL or matrix with dimnames[[1]] being a subset of c('0', positive integers that are item labels, 'log_delta1', 
# and 'log_delta2')). Each column is a set of fixed parameter values for the respective rowname. Leave any element of a column missing
# that you don't want to fix. So, for example, fixed = matrix(c(0, NA, 0, NA, 3, 3), nrow = 2,  byrow = T, dimnames = list(c("0","2"),NULL)),
# means that you want to run the algorithm when the fatigue parameter (theta0) is fixed at 0 and everything else is estimated, and you 
# want to run the algorithm when the item 2 parameter (theta2) is fixed at 3 and everything else is estimated, and you want to run the
# algorithm when both theta0 and theta2 are fixed at 0 and 3, respectively, with everything else to be estimated. This is all duplicated
# for each set of initial values provided. So, if num_init=2, this will imply six separate runs of the algorithm. 
#
### num_inits: (integer, positive) the number of initial starting values
#
### conv_tol: (numeric, positive) covergence tolerance; looking back over the previous 'stable_reps' (default 5) iterations 
# of the coordinate ascent algorithm, if if the maximum absolute difference between any current parameter value and the previous
# 'stable_reps' values of the same parameter is no more than 'conv_tol', then numerical convergence is satisifed (as long as 
# at least 'min_reps' [default 5] iterations have been run) 
#
### tau (numeric, non-negative) When positive, this makes the count smoothly varying between 0 and 1. 
#
### num_lambda (integer, positive) the maximum length of the grid of lambda values to evaluate over. The algorithm will move 
# increasing order, starting with the smallest value of lambda, but not all lambdas will necessarily be evaluated, however; 
# the algorithm will break after a minimum AIC is reached and it encounters a lambda that yields an AIC increase of 10% relative 
# to the minimum. See 'lambda_min_ratio'
#
### lambda_min_ratio (numeric, positive) if 'lambda_seq' is not provided, the algorithm will calculate its own. First, it
# estimates lambda_max that results in the most parsimonious model, then calculates lambda_min = lambda_max * lambda_min_ratio * tau^penalty_pow, 
# and then take the equal spaced (on the log-scale) sequence of values between these two extremes. 
#
### lambda_seq (numeric, positive) optional user-provided set of tuning parameters lambda. If not provided, the algorithm will
# calculate its own set.
#
### consider_truncating_lambda_seq (logical) defaults to TRUE. The algorithm starts with the smallest value of lambda and then 
# calculates the solution path for increasing lambda. It is expected that the AIC will go down but start to go up as lambda gets
# very large. Let lambda_best be the value in lambda_seq that gives the smallest (best) AIC. If consider_truncating_lambda_seq == TRUE, 
# the algorithm will truncate the all values in lambda_seq greater than the smallest value of lambda_star such that 
# (i) lambda_star > lambda_min and (ii) AIC(lambda_start) > 1.1 * AIC(lambda_min). 
#
### penalty_pow (numeric, positive) defaults to 1. This was not changed in the paper, but it is the power to which
# both 'tau' and each parameter are raised, which tends to make the shape of the seamless penalty more ridge-like close to zero. 
#
### proposal_seq (vector) the proposal sequence (Gamma in the manuscript) to use. The default behavior is as described in the manuscript. 
# If provided, at a minimum 'proposal_seq' should be a set of at least two numbers, with at least one element >= conv_tol and one element 
# <= -conv_tol. 
#
### proposal_seq_half_length (integer, positive) an alternative to specifying 'proposal_seq'. This tells the algorithm to use as
# the set of positive numbers in the proposal sequence the set of equal log-spaced values between 'conv_tol' and 'max_proposal_seq'. 
#
### max_proposal_seq (numeric, positive) the largest element of 'proposal_seq' to use
#
### stable_reps (integer, positive) the number of iterations to keep track of and look backwards over to determine when 
# convergence has been achieved. See also 'conv_tol'
#
### min_reps (integer, positive) the minimum number of iterations to conduct, regardless of convergence
#
### multivariable_proposals (logical) should the multivariable proposals described in the manuscript be conducted? Some people are
# concerned about introducing a stochastic element into the CA algorithm, and setting this to FALSE removes the stochastic element of the
# optimization. Empirically, leaving it as TRUE speeds up convergence without changing any results. 
#
### delta_lb (numeric, positive) this is a parameter that provides a hard lower bound for the values of delta1 and delta2. In theory this
# should not be needed, i.e. it should just be left at its natural lower bound of 0, but the algorithm seemed to misbehave sometimes without 
# a positive lower bound. 
#
### max_reps_kgt1 (integer, positive) the maximum allowed number of iterations for all but the smallest value of lambda. 
#
### max_reps_keq1 (integer, positive) the maximum allowed number of iterations for the smallest value of lambda. 
#
### random_seed (integer, positive) the random seed to use for initializing the starting parameter values
#
### verbose (logical) if TRUE, the algorithm's progress is communicated. Tends to slow down everything
#
### safe_mode: (logical) if TRUE, then perform error checks to make sure that all inputs are as expected. Set to FALSE to speed things up a bit
#
### init_params: (data.frame with 'row.names' attribute) the parameter vector to evaluate with respect to the data. The number of 
# columns implicitly gives the number of starting values to look over. Must have rownames exactly as follows: '0' (optional, if the stopping process 
# is to be modeled), then consecutive positive integer labels, one for each unique item, then 'log_delta1', then 'log_delta2'
#
### Value
#
###A list of many results from the fitting process:
## control = a list of the parameters used by the algorithm 
## monitor_joint_shifts = a matrix with dimension c(length(lambda_seq), 4), giving the counts of the number of each type of multivariable
# proposal that was accepted. This is just to see if the multivariable proposals are actually offering sensible proposals
## hit_max_reps = matrix with num_inits * ncol(fixed) rows and num_lambda columns, indicating when the algorithm quit due to hitting 
# the maximum allowed number of iterations (either max_reps_keq1 for the first lambda or max_reps_kgt1 for other lambdas) instead
# of converging?
## actual_reps = matrix with the same dimension as hit_max_reps indicating the number of iterations till converging or quitting
## num_eff_params = matrix with the same dimension as hit_max_reps giving the effective number of parameters, equivalently the L0 penalty
# before being scaled by lambda
## num_actual_params = atrix with the same dimension as hit_max_reps giving the actual number of parameters in the model
## all_neg_loglikelihood = matrix with the same dimension as hit_max_reps giving the negative log-likelihood values
## all_bic = matrix with the same dimension as hit_max_reps giving the bic values
## all_aic = matrix with the same dimension as hit_max_reps giving the aic values
## all_params = list as long as num_inits * ncol(fixed), with each element being a matrix with number of rows equal to the largest possible
# number of parameters and number of columsn equal to the number of lambdas evaluated. 
## best_fit_neg_loglikelihood = matrix with number of rows equal to the largest possible number of parameters and number of columns equal to
# num_inits * ncol(fixed), giving the best solution according to maximized log-likelihood
## best_neg_loglikelihood = matrix with four rows and number of columns equal to num_inits * ncol(fixed), giving the value of the negative
# log-likelihood, the value of lambda, the number of effective parameters, and the number of actual parameters corresponding to the fit with
# the largest log-likelihood
## best_fit_bic = similar to best_fit_neg_loglikelihood but using BIC as the objective function
## best_bic = similar to best_neg_loglikelihood but using BIC as the objective function
## best_fit_aic = imilar to best_fit_neg_loglikelihood but using AIC as the objective function 
## best_aic =  similar to best_neg_loglikelihood but using AIC as the objective function
## runtime = the total runtime of the function

penRank_path = function(dat, 
                        obs_weights = rep(1, nrow(dat)), 
                        fixed = NULL, 
                        num_inits = 4, 
                        conv_tol = 1e-3, 
                        tau = NULL,
                        num_lambda = 200, 
                        lambda_min_ratio = 1e-5, 
                        lambda_seq = NULL, 
                        consider_truncating_lambda_seq = TRUE,
                        penalty_pow = 1,
                        proposal_seq = NULL, 
                        proposal_seq_half_length = NULL,
                        max_proposal_seq = 1,
                        stable_reps = 5, 
                        min_reps = 5,
                        multivariable_proposals = TRUE,
                        delta_lb = 0.25, 
                        max_reps_kgt1 = 1e3, 
                        max_reps_keq1 = 250, 
                        random_seed = sample(.Machine$integer.max,1),
                        verbose = FALSE, 
                        safe_mode = FALSE, 
                        init_params = NULL) {
  
  # + Initialize ----
  begin = Sys.time();
  min_reps = max(min_reps, stable_reps);
  #All parameter estimates are always rounded to the specified convergence tolerance
  sig_digits = ceiling(abs(log10(conv_tol)));
  #Translate the lower bound of the dampening function's parameters appropriately;
  log_delta_lb = floor(log(delta_lb) * (10^sig_digits))/(10^sig_digits);
  #Numerical tolerance for zero;
  tiny_positive = .Machine$double.eps^0.5;
  num_obs = nrow(dat);
  num_ranks = ncol(dat) + 1;
  num_uniq = sum(!is.na(unique(as.numeric(dat))));
  max_columns = min(num_uniq, num_ranks);
  item_names = as.character(1:num_uniq);
  if(any(is.na(dat))) {
    param_names = c("0",item_names,"log_delta1","log_delta2");
  } else {
    param_names = c(item_names,"log_delta1","log_delta2");
  }
  
  #Check for invalid / problematic arguments
  if(tiny_positive > conv_tol || conv_tol > 1) {stop(c("'conv_tol' should be between 'tiny_positive', i.e. ",formatC(tiny_positive,format="e",digits = 2),", and 1"));}
  if(num_obs != length(obs_weights)) {stop("'obs_weights' must have length equal to 'nrow(dat)'");}
  if(num_uniq != max(dat,na.rm = TRUE)) {stop("Item labels in 'dat' must be consecutively labeled integers starting with 1");}
  if(!isTRUE(all(apply(apply(is.na(cbind(dat,NA)),1,cumsum),2,diff) | (apply(is.na(dat),1,cumsum) == 0))) | 
     isTRUE(any(apply(dat,1,duplicated) & !apply(dat,1,is.na)))) {
    stop("'dat' must contain uninterrupted left-to-right sets of unique integers, i.e. no NAs between items");
  }
  if(!is.null(init_params) && !all(rownames(init_params) == param_names)) {stop("If provided, 'init_params' must have rownames exactly as follows: '0' (but only if any lists stop early), then consecutive positive integer labels, one for each unique item, then 'log_delta1', then 'log_delta2'");}
  if(!is.null(init_params) && sum(init_params[item_names,] < 0)) {stop("There are negative item weights in the user-provided matrix of 'init_params'; all item weights should be non-negative")}
  if((is.null(rownames(fixed)) || !all(rownames(fixed)%in%param_names)) && !is.null(fixed)) {stop("'fixed' must be NULL or a matrix with rownames being a subset of the parameter names");}
  if(is.null(proposal_seq) && !is.null(max_proposal_seq) && max_proposal_seq < conv_tol) {stop("'max_proposal_seq' must be greater than 'conv_tol'");}
  if(is.null(proposal_seq) && !is.null(proposal_seq_half_length) && (proposal_seq_half_length < 1 || proposal_seq_half_length %% 1 != 0)) {stop("'proposal_seq_half_length' must be a positive integer");}
  if(!is.null(proposal_seq) && !(any(proposal_seq <= -conv_tol) && any(proposal_seq >= conv_tol))) {stop(paste0("At a minimum, 'proposal_seq' should be a set of at least two numbers, with at least one element >= ", formatC(conv_tol,format="e", digits = 2), " and one element <= ",  formatC(-conv_tol,format="e", digits = 2)));}
  if(penalty_pow <= 0) {stop("'penalty_pow' must be strictly positive");}
  if(!is.null(tau) && tau <= 0) {stop("'tau', if provided, must be strictly positive");}
  
  #Additional setup
  if(is.null(tau)) {
    tau = conv_tol;
  }
  if(!is.null(lambda_seq)) {
    num_lambda = length(lambda_seq);
  } 
  length_params = length(param_names);
  list_lengths = unlist(lapply(apply(apply(cbind(dat,NA),1,is.na),2,which),"[",1)) - 1;
  sum_obs_weights = sum(obs_weights);
  
  #If fixed is null, then all parameters are to be estimated
  if(is.null(fixed)) {
    fixed = matrix(numeric(0));
    #If init_params is provided but does not match num_inits, force it to be so;
    if(!is.null(init_params) && ncol(init_params) != num_inits) {
      num_inits = ncol(init_params);
      cat("Forcing 'num_inits = ncol(init_params)'\n");
    }
    #If init_params is provided but does not match num_inits, force it to be so and note that any provided initial values that are also in fixed will be overwritten
  } else if(!is.null(init_params) && ncol(init_params) != num_inits) {
    num_inits = ncol(init_params);
    cat("Forcing 'num_inits = ncol(init_params)' and overwriting all values in rows of 'init_param' that are provided in matched rows of 'fixed'\n");
  }
  if(!(length(multivariable_proposals) %in% c(1, num_inits))) {stop("'multivariable_proposals' should either be a single logical or a vector of logicals as long as the number of initial values ('num_inits')");}
  multivariable_proposals = rep(multivariable_proposals, length = num_inits * ncol(fixed));
  
  #Set up 'init_params'
  set.seed(random_seed);  
  if(is.null(init_params)) {
    init_params = matrix(0, nrow = length_params,ncol = num_inits * ncol(fixed), dimnames = list(param_names,NULL));
    if("0"%in%param_names) {
      init_params["0", ] = rnorm(num_inits, sd = 2.5);
    }
    init_params[item_names, ] = rexp(num_uniq * num_inits, rate = 0.5);
    init_params["log_delta1",] = rep(log(runif(num_inits, delta_lb, 1)), each = ncol(init_params) / num_inits);
    init_params["log_delta2",] = rep(log(runif(num_inits, delta_lb, 1)), each = ncol(init_params) / num_inits);
  } else {
    #If 'init_params' is provided and fixed is also provided, 'init_params' will be replicated, once for each column of fixed
    init_params = matrix(rep(init_params,ncol(fixed)),nrow = length_params);
  }
  #'fixed_indicators' matches the dimensions of 'init_params' and indicates which parameters in which column are to be left alone, i.e. fixed
  fixed_indicators = matrix(0, nrow = nrow(init_params), ncol = ncol(init_params), dimnames = list(param_names,NULL));
  if(length(fixed)) {
    for(m in rownames(fixed)) {
      foo = rep(fixed[m,], each = num_inits);
      #This is the overwriting step; any fixed parameters are initialized as such and left alone
      init_params[m,] = ifelse(is.na(foo),init_params[m,], foo);
      fixed_indicators[m,] = ifelse(is.na(foo),fixed_indicators[m,], 1);
    }
    rm(foo);
  } 
  #Start off everything so that the smallest initial item weight is always exactly zero, while ensuring that fixed parameters are not overwritten
  if(!any(item_names %in% rownames(fixed))) {
    init_params[item_names,] = init_params[item_names,,drop = FALSE] - matrix(apply(init_params[item_names,,drop = FALSE],2,min),nrow = num_uniq, ncol = ncol(init_params),byrow = TRUE);
  } else {
    for(i in 1:ncol(init_params)) {
      item_names_to_estimate = setdiff(item_names,param_names[which(fixed_indicators[,i] == 1)]);
      init_params[item_names_to_estimate,i] = init_params[item_names_to_estimate,i,drop = FALSE] - min(init_params[item_names,i,drop = FALSE]);
    }
    rm(item_names_to_estimate);
  }
  #Check that the dampening parameters are in their proper support; note: the random initial values are constructed to always pass this check
  if(any(init_params["log_delta1",] > 0) || any(init_params["log_delta1",] < log_delta_lb)) {stop("All initial values of 'log_delta1' must be in ['log_delta_lb', 0]");} 
  if(any(init_params["log_delta2",] > 0) || any(init_params["log_delta2",] < log_delta_lb)) {stop("All initial values of 'log_delta2' must be in ['log_delta_lb', 0]");} 
  init_params = round(init_params, sig_digits);
  
  #Calculate this now to avoid re-calculating at each call to the lower-level functions
  all_item_counts_by_rank = matrix(NA, nrow = num_uniq, ncol = max_columns, dimnames = list(item_names, NULL));
  for(m in as.numeric(item_names)) {  
    all_item_counts_by_rank[m,] = c(colSums(obs_weights * (dat == m), na.rm = TRUE), 0)[1:max_columns]
  }
  
  store_bic = store_aic = matrix(Inf, nrow = num_inits * ncol(fixed), ncol = num_lambda);
  hit_max_reps = actual_reps = matrix(0, nrow = num_inits * ncol(fixed), ncol = num_lambda);
  best_fit_bic = best_fit_aic = best_fit_neg_loglikelihood = 
    init_params * 0;
  store_params = vector("list",num_inits * ncol(fixed));
  
  if(is.null(proposal_seq)) {
    if(is.null(proposal_seq_half_length)) {
      proposal_seq_half_length = max(2, sig_digits);
    }
    proposal_seq = round(exp(seq(log(conv_tol), log(max_proposal_seq), length = proposal_seq_half_length)), sig_digits);
    proposal_seq = sort(c(-proposal_seq, proposal_seq));
  } 
  
  positive_proposal_seq = proposal_seq[proposal_seq > tiny_positive];
  
  #The initial 'lambda_seq' is is equally spaced (on the log-scale) with distance 1/2, and then after an initial 
  #lambda_max is identified, it is further refined by looking at equally spaced increments of 1/8 above and below
  if(is.null(lambda_seq)) {
    #Remember all of the true values of these tuning parameters; we will reset them after 'lambda_seq' has been constructed
    identifying_max_lambda = 
      first_pass = TRUE;
    num_lambda_requested = num_lambda;
    stable_reps_actual = stable_reps;
    min_reps_actual = min_reps;
    tau_actual = tau;
    stable_reps = min_reps = 3; 
    tau = min(1e-2, tau);
    construction_message = "\n Constructing 'lambda_seq' (first pass) \n";
    lambda_seq = exp(seq(log(length_params * sum_obs_weights),log(length_params * sum_obs_weights) + 2, length = 5)); 
  } else {
    identifying_max_lambda = FALSE;
  }
  #There are potentially multiple passes through this top loop. All but the last pass are for constructing the value of 'lambda_seq'.
  #Thus if it's provided, there will only be one pass through this top loop
  while(T) {
    num_lambda = length(lambda_seq);
    store_num_eff_params = 
      store_num_actual_params = 
      store_neg_loglikelihood = matrix(Inf, nrow = num_inits * ncol(fixed), ncol = num_lambda);
    #the values of j_start and j_end depend upon whether this 'lambda_seq' is still being constructed or if this is the actual
    #solution path. If the former, then we don't care about the multiple initial starts and only focus on the first. 
    j_start = ifelse(identifying_max_lambda, which.min(colSums(fixed_indicators)), 1);
    j_end = ifelse(identifying_max_lambda, j_start, num_inits * ncol(fixed));
    j = j_start;
    #Count the number of times that a multivariate proposal was accepted
    monitor_joint_shifts = matrix(0,  nrow = 4, ncol = num_lambda, dimnames = list(c("shift up","shift down", "permute estimates","swap delta"),lambda_seq));
    
    #Loop through the initial starting values of the parameters
    for(j in j_start:j_end) {
      min_possible_num_param = seamless_L0(init_params[,j] * fixed_indicators[,j], tau = tau, penalty_pow = penalty_pow, safe_mode = TRUE, num_uniq = num_uniq);
      which_to_estimate = setdiff(param_names,param_names[which(fixed_indicators[,j] == 1)]);
      which_to_estimate_except_dampening = setdiff(which_to_estimate,c("log_delta1","log_delta2"));
      which_dampening_to_estimate = intersect(which_to_estimate,c("log_delta1","log_delta2"));
      which_items_to_estimate = setdiff(which_to_estimate,c("0","log_delta1","log_delta2"));
      num_items_to_estimate = length(which_items_to_estimate);
      
      store_params[[j]] = matrix(NA, nrow = length_params, ncol = num_lambda, dimnames = list(param_names,formatC(lambda_seq,format="f",digits=3)));
      if(!identifying_max_lambda) {
        cat("\n Initial Value ", j, "\n");
      } else {
        cat(construction_message); 
      }
      param_for_compare = matrix(NA,nrow = length_params, ncol = stable_reps + 1, dimnames = list(param_names,NULL));
      param_for_compare[,1] = init_params[,j];
      k = 1;
      while(k <= num_lambda) {
        if(verbose && !identifying_max_lambda) {cat("\n====",j,"--",k,"====\n");}
        lambda = lambda_seq[k];
        i = curr_column = 1;
        previous_column =  stable_reps + 1;
        #curr_best_cost keeps track of the smallest achieved penalized log-likelihood that we've achieved so far with this lambda;
        #the corresponding parameter values that yield this best cost will be kept in param_for_compare[,curr_column];
        curr_best_cost = BPL_neg_loglik_CA(param_for_compare[,curr_column],NA,dat,obs_weights,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths,all_item_counts_by_rank[m,]) + 
          lambda * seamless_L0(param_for_compare[,curr_column],NA, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq) + 
          over_identified_cost(param_for_compare[,curr_column],NA, safe_mode = safe_mode, tiny_positive = tiny_positive);
        while((i <= min_reps) || max(abs(param_for_compare[which_to_estimate,-curr_column] - param_for_compare[which_to_estimate,curr_column])) >= conv_tol - tiny_positive) {
          #First consider coordinate-wise maximization (skipping if the previous full iteration didn't result in any changes)
          if((i <= min_reps) || max(abs(param_for_compare[which_to_estimate,ifelse(previous_column == 1, stable_reps + 1, previous_column - 1)] - 
                                        param_for_compare[which_to_estimate,curr_column])) > tiny_positive) {
            curr_column = 1 + (curr_column)%%(stable_reps + 1);
            previous_column = 1 + (previous_column)%%(stable_reps + 1);
            curr_try = param_for_compare[,previous_column];
            for(m in sample(which_to_estimate)) {
              #start with the current values, and for exactly one parameter, propose a new sequence of values to consider, subject to the natural bounds, 
              #and always ensuring that the 'null' value is proposed  
              if(m %in% c("log_delta1","log_delta2")) {
                CA_curr_try = unique(pmax(log_delta_lb, pmin(0, curr_try[m] + proposal_seq)));
              } else if(m == "0") {
                CA_curr_try = curr_try[m] + proposal_seq;
              } else {
                CA_curr_try = unique(pmax(0, c(0, curr_try[m] + proposal_seq)));
              }
              curr_try[m] = NA;
              curr_try_cost = BPL_neg_loglik_CA(curr_try,CA_curr_try,dat,obs_weights,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths,all_item_counts_by_rank[m,]) + 
                lambda * seamless_L0(curr_try,CA_curr_try, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq) + 
                over_identified_cost(curr_try,CA_curr_try, safe_mode = safe_mode, tiny_positive = tiny_positive);
              if(curr_try_cost[which_min_curr_try_cost <- which.min(curr_try_cost)] < curr_best_cost) {
                curr_try[m] = CA_curr_try[which_min_curr_try_cost];
                curr_best_cost = curr_try_cost[which_min_curr_try_cost];
              } else {
                curr_try[m] = param_for_compare[m,previous_column]
              }
            }
            param_for_compare[,curr_column] = curr_try
          } else {
            if(verbose) {cat("skipping univariate CA\n\n");}
            curr_column = 1 + (curr_column)%%(stable_reps + 1);
            previous_column = 1 + (previous_column)%%(stable_reps + 1);
            param_for_compare[,curr_column] = param_for_compare[,previous_column];
          }
          if(multivariable_proposals[j]) {
            #Now consider some multivariable adjustments to discourage getting stuck at local maxima
            shift_up = shift_down = permute_estimates = swap_delta =
              param_for_compare[,curr_column];
            #shift_down: Bring all item weights down by a constant amount (but item parameters should be non-negative); make dampening parameter closer to zero
            shift_down_amount = sample(positive_proposal_seq,size = 1)
            shift_down[which_items_to_estimate] = pmax(0, shift_down[which_items_to_estimate] - shift_down_amount);
            #shift_up: Shift up all positive weights by a constant amount and one random zero weight by that same amount (but only if there is more than one); make dampening parameter closer to zero
            shift_up_amount = sample(positive_proposal_seq,size = 1);
            zero_items = item_names[which(param_for_compare[item_names,curr_column] <= tiny_positive)];
            if(length(zero_items) > 1 && length(estimated_zero_items <- intersect(which_items_to_estimate, zero_items))) {
              items_to_shift_up = c(sample(estimated_zero_items, 1), setdiff(which_items_to_estimate, zero_items));
            } else {
              items_to_shift_up = setdiff(which_items_to_estimate, zero_items);
            }
            shift_up[items_to_shift_up] = shift_up[items_to_shift_up] + shift_up_amount;
            #permute_estimates: randomly permute estimated item weights up or down one place; make dampening parameter closer to zero
            old_order = order(permute_estimates[which_items_to_estimate]);
            new_order = old_order[order((1:num_items_to_estimate) + sample(c(-1,0,1), num_items_to_estimate, replace = TRUE,prob = c(1,2,1)))];

            permute_estimates[which_items_to_estimate[old_order]] = permute_estimates[which_items_to_estimate[new_order]];
            rm(new_order,old_order);
            
            if("0" %in% which_to_estimate) {
              shift_up["0"] = shift_up["0"] + sample(proposal_seq,size = 1);
              shift_down["0"] = shift_down["0"] + sample(proposal_seq,size = 1);
              permute_estimates["0"] = permute_estimates["0"] + sample(proposal_seq,size = 1);
            }
            
            if(length(which_dampening_to_estimate)) {
              shift_up[which_dampening_to_estimate] =  
                pmax(log_delta_lb, pmin(0, shift_up[which_dampening_to_estimate] + sample(proposal_seq,size = 1)));
              shift_down[which_dampening_to_estimate] =  
                pmin(0, shift_down[which_dampening_to_estimate] + sample(positive_proposal_seq,size = 1));
              permute_estimates[which_dampening_to_estimate] =  
                pmin(0, permute_estimates[which_dampening_to_estimate] + sample(positive_proposal_seq,size = 1));
            }
            
            #Calculate, store the penalized log-likelihoods of these multivariable proposals.
            shift_down_cost = BPL_neg_loglik(shift_down,dat,obs_weights,fixed = NULL,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths) + 
              lambda * seamless_L0(shift_down, NA, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq) + 
              over_identified_cost(shift_down, NA, safe_mode = safe_mode, tiny_positive = tiny_positive);
            shift_up_cost = BPL_neg_loglik(shift_up,dat,obs_weights,fixed = NULL,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths) + 
              lambda * seamless_L0(shift_up, NA, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq) + 
              over_identified_cost(shift_up, NA, safe_mode = safe_mode, tiny_positive = tiny_positive);
            permute_estimates_cost = BPL_neg_loglik(permute_estimates, dat, obs_weights, fixed = NULL,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths) + 
              lambda * seamless_L0(permute_estimates, NA, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq) + 
              over_identified_cost(permute_estimates, NA, safe_mode = safe_mode, tiny_positive = tiny_positive);
            
            #swap_delta: if both dampening parameters are to be estimated, consider a simpler dampening model that sets log_delta2 = 0 while providing an equivalent amount
            #of dampening in the first two stages
            if(length(which_dampening_to_estimate) == 2) {
              stage1_dampening = exp(sum(swap_delta[c("log_delta1","log_delta2")])) + (1-exp(swap_delta["log_delta2"]))^3;
              swap_delta[c("log_delta1","log_delta2")] =  
                pmax(log_delta_lb, pmin(0,c(round(log(stage1_dampening),sig_digits),0)));
              swap_delta_cost = BPL_neg_loglik(swap_delta,dat,obs_weights,fixed = NULL,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths) + 
                lambda * seamless_L0(swap_delta, NA, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq) + 
                over_identified_cost(swap_delta, NA, safe_mode = safe_mode, tiny_positive = tiny_positive);
            } else {
              swap_delta_cost = Inf;
            }
            
            if(permute_estimates_cost < min(swap_delta_cost, shift_down_cost, shift_up_cost, curr_best_cost)) {
              if(verbose) {cat("permuting estimates\n\n");}
              curr_best_cost = permute_estimates_cost;
              param_for_compare[,curr_column] = permute_estimates;
              monitor_joint_shifts["permute estimates", k] = monitor_joint_shifts["permute estimates", k] + 1;
            } else if(swap_delta_cost < min(shift_down_cost, shift_up_cost, curr_best_cost)) {
              if(verbose) {cat("swapping delta\n\n");}
              curr_best_cost = swap_delta_cost;
              param_for_compare[,curr_column] = swap_delta;
              monitor_joint_shifts["swap delta", k] = monitor_joint_shifts["swap delta", k] + 1;
            } else if (shift_down_cost < min(shift_up_cost, curr_best_cost)) {
              if(verbose) {cat("shifting down\n\n");}
              curr_best_cost = shift_down_cost;
              param_for_compare[,curr_column] = shift_down;
              monitor_joint_shifts["shift down", k] = monitor_joint_shifts["shift down", k] + 1;
            } else if(shift_up_cost < curr_best_cost) {
              if(verbose) {cat("shifting up\n\n");}
              curr_best_cost = shift_up_cost;
              param_for_compare[,curr_column] = shift_up;
              monitor_joint_shifts["shift up", k] = monitor_joint_shifts["shift up", k] + 1;
            }
          }
          if(i >= ifelse(k == 1, max_reps_keq1, max_reps_kgt1)) {break;}
          i = i + 1;
        }
        if(verbose && !identifying_max_lambda) {
          if(k == 1) { 
            cat(i - 1,formatC((param_for_compare[which_to_estimate,curr_column] - init_params[which_to_estimate,j]),format="f"),"\n");
          } else {
            cat(i - 1,formatC((param_for_compare[which_to_estimate,curr_column] - store_params[[j]][which_to_estimate,k-1]),format="f"),"\n");
          }
        }
        store_num_actual_params[j,k] = seamless_L0(param_for_compare[,curr_column],CA_param = NULL, tau = tiny_positive, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq);
        store_num_eff_params[j,k] = seamless_L0(param_for_compare[,curr_column],CA_param = NULL, tau = tau, penalty_pow = penalty_pow, safe_mode = safe_mode, num_uniq = num_uniq);
        #If we are actually calculating the solution path (as opposed to constructing 'lambda_seq'), then store the results
        if(!identifying_max_lambda) {
          store_neg_loglikelihood[j,k] = BPL_neg_loglik(param_for_compare[,curr_column],dat,obs_weights,fixed = NULL,safe_mode = safe_mode, num_uniq, num_obs, num_ranks, max_columns, list_lengths)
          store_bic[j,k] = 2 * store_neg_loglikelihood[j,k] + 
            log(sum_obs_weights) * store_num_actual_params[j,k];
          store_aic[j,k] = 2 * store_neg_loglikelihood[j,k] + 
            2 * store_num_actual_params[j,k] * (sum_obs_weights / pmax(tiny_positive, sum_obs_weights -  store_num_eff_params[j,k] - 1));
          store_params[[j]][,k] =  param_for_compare[,curr_column]; 
          actual_reps[j,k] = i - 1;
          if(i >= ifelse(k == 1, max_reps_keq1, max_reps_kgt1)) {
            hit_max_reps[j,k] = 1;
          }
          if(consider_truncating_lambda_seq && store_aic[j,k] > (1.1 * min(store_aic[j,]))) {
            k = num_lambda;
          } 
        }
        k = k + 1;
      }
      if(!identifying_max_lambda) {
        store_params[[j]] = round(store_params[[j]], sig_digits);
        best_fit_bic[,j] = store_params[[j]][,which.min(store_bic[j,])];
        best_fit_aic[,j] = store_params[[j]][,which.min(store_aic[j,])];
        best_fit_neg_loglikelihood[,j] = store_params[[j]][,which.min(store_neg_loglikelihood[j,])];
      }
    }
    #If we are actually calculating the solution path, then we're done!
    if(!identifying_max_lambda) {
      break;
    } else if(first_pass) {
      #Grow the largest value of lambda until maximum parsimony is achieved
      which_lambda_max = min(c(which(store_num_eff_params[j,1:num_lambda] <= min_possible_num_param + tiny_positive)),num_lambda + 1);
      if(which_lambda_max == 1) {
        lambda_seq = c(exp(seq(log(lambda_seq[1]) - 2.5, log(lambda_seq[1]) - 0.5, length =  5)),lambda_seq);
      } else if(which_lambda_max == num_lambda + 1) {
        lambda_seq = c(lambda_seq, exp(seq(log(lambda_seq[which_lambda_max - 1]) + 0.5, log(lambda_seq[which_lambda_max - 1]) + 2.5, length =  5)));
      } else {
        lambda_seq = exp(seq(log(lambda_seq[1]), log(lambda_seq[which_lambda_max]), by =  1/4));
        first_pass = FALSE;
        construction_message = "\n Constructing 'lambda_seq' (second pass) \n";
      }
      rm(which_lambda_max);
    } else {
      #Refine the value of lambda_max to get a more precise estimate of the smallest value that still yields maximum parsimony
      which_lambda_max = min(c(which(store_num_eff_params[j,1:num_lambda] <= min_possible_num_param + tiny_positive)),num_lambda + 1);
      if(which_lambda_max == 1) {
        log_lambda_max = log(lambda_seq[1]) - 1/4;
      } else if(which_lambda_max == num_lambda + 1) {
        log_lambda_max = log(lambda_seq[num_lambda]) + 1/4;
      } else {
        log_lambda_max = log(lambda_seq[which_lambda_max]);
      }
      lambda_seq = exp(seq(log_lambda_max + log(lambda_min_ratio) + penalty_pow * log(tau_actual), log_lambda_max, length = num_lambda_requested));
      identifying_max_lambda = FALSE;
      stable_reps = stable_reps_actual;
      min_reps = min_reps_actual;
      tau = tau_actual;
      store_params = vector("list",num_inits * ncol(fixed));
      rm(log_lambda_max, which_lambda_max, stable_reps_actual, min_reps_actual, num_lambda_requested);
    }
  }
  #If any of the prepopulated results weren't filled in because 'consider_truncating_lambda_seq == TRUE', then drop these
  if(any(colSums(store_neg_loglikelihood == Inf) == nrow(store_neg_loglikelihood))) {
    keep_columns = which(colSums(store_neg_loglikelihood == Inf) < nrow(store_neg_loglikelihood));
    lambda_seq = lambda_seq[keep_columns, drop = FALSE];
    store_bic = store_bic[, keep_columns, drop = FALSE];
    store_aic = store_aic[, keep_columns, drop = FALSE];
    store_neg_loglikelihood = store_neg_loglikelihood[, keep_columns, drop = FALSE];
    store_num_eff_params = store_num_eff_params[, keep_columns, drop = FALSE];
    store_num_actual_params = store_num_actual_params[, keep_columns, drop = FALSE];
    monitor_joint_shifts = monitor_joint_shifts[, keep_columns, drop = FALSE];
    hit_max_reps = hit_max_reps[, keep_columns, drop = FALSE];
    actual_reps = actual_reps[, keep_columns, drop = FALSE];
    for(j in 1:length(store_params)) {
      store_params[[j]] = store_params[[j]][, keep_columns, drop = FALSE];
    }
  }
  
  colnames(store_aic) = 
    colnames(store_bic) = 
    colnames(monitor_joint_shifts) = 
    colnames(store_neg_loglikelihood) = 
    colnames(store_num_eff_params) = 
    colnames(hit_max_reps) = 
    colnames(actual_reps) = 
    formatC(lambda_seq, format = "f", digits = sig_digits);
  
  
  best_neg_loglikelihood = rbind(apply(store_neg_loglikelihood,1,min), 
                                 lambda_seq[apply(store_neg_loglikelihood,1,which.min)],
                                 store_num_eff_params[cbind(1:nrow(store_num_eff_params),apply(store_neg_loglikelihood,1,which.min))],
                                 store_num_actual_params[cbind(1:nrow(store_num_actual_params),apply(store_neg_loglikelihood,1,which.min))])
  best_bic = rbind(apply(store_bic,1,min), 
                   store_neg_loglikelihood[cbind(1:nrow(store_num_eff_params),apply(store_bic,1,which.min))],
                   lambda_seq[apply(store_bic,1,which.min)],
                   store_num_eff_params[cbind(1:nrow(store_num_eff_params),apply(store_bic,1,which.min))],
                   store_num_actual_params[cbind(1:nrow(store_num_actual_params),apply(store_bic,1,which.min))]);
  best_aic = rbind(apply(store_aic,1,min), 
                   store_neg_loglikelihood[cbind(1:nrow(store_num_eff_params),apply(store_aic,1,which.min))],
                   lambda_seq[apply(store_aic,1,which.min)],
                   store_num_eff_params[cbind(1:nrow(store_num_eff_params),apply(store_aic,1,which.min))],
                   store_num_actual_params[cbind(1:nrow(store_num_actual_params),apply(store_aic,1,which.min))]);
  rownames(best_neg_loglikelihood) = 
    c("neg_loglikelihood", "lambda","num_eff_params", "num_actual_params");
  rownames(best_bic) = 
    c("bic","neg_loglikelihood", "lambda","num_eff_params", "num_actual_params");
  rownames(best_aic) =
    c("aic", "neg_loglikelihood","lambda","num_eff_params", "num_actual_params");
  
  if(num_lambda > 1 && (any(apply(store_bic,1,which.min) == 1) || any(apply(store_aic,1,which.min) == 1))) {
    cat("The smallest AIC and/or BIC occurred at the minimum value of lambda; if this is unexpected, consider decreasing the minimum value of lambda");
  }
  return(list(control = list(conv_tol = conv_tol,
                             tau = tau, 
                             penalty_pow = penalty_pow,
                             lambda_seq = lambda_seq, 
                             consider_truncating_lambda_seq = consider_truncating_lambda_seq,
                             proposal_seq = proposal_seq, 
                             stable_reps = stable_reps, 
                             min_reps = min_reps, 
                             max_reps_keq1 = max_reps_keq1, 
                             max_reps_kgt1 = max_reps_kgt1,
                             random_seed = random_seed,
                             delta_lb = delta_lb,
                             init_params = init_params, 
                             fixed = fixed),
              monitor_joint_shifts = monitor_joint_shifts,
              hit_max_reps = hit_max_reps, 
              actual_reps = actual_reps, 
              num_eff_params = store_num_eff_params, 
              num_actual_params = store_num_actual_params, 
              all_neg_loglikelihood = store_neg_loglikelihood,
              all_bic = store_bic, 
              all_aic = store_aic, 
              all_params = store_params, 
              best_fit_neg_loglikelihood = best_fit_neg_loglikelihood,
              best_neg_loglikelihood = best_neg_loglikelihood, 
              best_fit_bic = best_fit_bic, 
              best_bic = best_bic, 
              best_fit_aic = best_fit_aic, 
              best_aic = best_aic, 
              runtime = Sys.time() - begin));
}


# function name and purpose: 'bpl_probs' calculates the model-based probabilities of based upon the provided parameter values and the
# given stage. Useful for simulating from fitted models. 
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
### stage: (positive integer) what stage should the probabilities be simulated at? Remember that the dampening function implies the 
# same parameter values will give different probabilities at different stages. 
#
### Value
#
# A named vector of probabilities with length equal to length(param) - 2 containing one probability for each element of param except for
# 'log_delta1' and 'log_delta2'. The sum of the vector equals 1, and these are the model-estimated probabilities of choosing the corresponding
# item (including item 0) at the next stage.

bpl_probs = function(param, stage = 1) {
  param_names = names(param);
  item_names = setdiff(param_names, c("0","log_delta1","log_delta2"));
  fatigue_weight = ifelse("0"%in%param_names && stage > 1, exp(param["0"]), 0);
  fitted_probs = numeric(length(param) - 2);
  names(fitted_probs) = setdiff(param_names, c("log_delta1","log_delta2"))
  dampening = exp(param["log_delta2"] + (stage - 1) * param["log_delta1"]) + 
    (1-exp(param["log_delta2"]))^(2 * stage - 1) ;
  fitted_probs[item_names] = exp(dampening * param[item_names]) / (fatigue_weight + sum(exp(dampening * param[item_names])))
  if("0"%in%param_names) {fitted_probs["0"] = 1 - sum(fitted_probs);}
  return(fitted_probs);
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
    result[,1] = sample(setdiff(names(param), c("log_delta1","log_delta2")),size = n_to_sim,  prob = bpl_probs(param, stage = 1), replace = TRUE);
    for(curr_sim in 1:n_to_sim) {
      for(curr_stage in 2:length(stages)) {
        curr_param = param[setdiff(names(param),result[curr_sim,1:curr_stage])];
        result[curr_sim,curr_stage] = sample(setdiff(names(curr_param), c("log_delta1","log_delta2")),size = 1,  prob = bpl_probs(curr_param, stage = curr_stage), replace = TRUE);
        if(result[curr_sim,curr_stage]=="0") {
          break;
        }
      }
    }
  }
  apply(result, 2, as.numeric);
}

ordered_barplot = function(data_ranked, 
                           col_palette_range, 
                           file_path = "foo.pdf", 
                           title_text = "",
                           skip_lonely = T) {
  require(RColorBrewer);
  if(skip_lonely && any(colSums(!is.na(data_ranked)) == 1)) {
    data_ranked_lonely = data_ranked[,which(colSums(!is.na(data_ranked)) == 1)];
    lonely_column_order_to_plot = order(-as.numeric(colSums(!is.na(data_ranked_lonely))),as.numeric(apply(data_ranked_lonely,2,mean,na.rm=T)));
    data_ranked = data_ranked[,which(colSums(!is.na(data_ranked)) > 1)];
  } else {
    data_ranked_lonely = data.frame(NULL);
  }
  
  column_order_to_plot = order(-as.numeric(colSums(!is.na(data_ranked))),as.numeric(apply(data_ranked,2,mean,na.rm=T)));
  
  xmax = max(max(colSums(!is.na(data_ranked))), sum(!is.na(data_ranked_lonely)));
  ymax = ncol(data_ranked) + as.logical(ncol(data_ranked_lonely));
  
  col_func = colorRampPalette(col_palette_range);
  plot_colors = col_func(5);
  plot_pch = c(15, 16, 17, 18, 15);
  plot_cex = c(1.15, 1.25, 1.05, 1.10, 0.85);
  
  pdf(file = file_path, 
      width = 8, 
      height = 1.4 + 2.6 * (ymax/28), 
      family ="serif");
  par(mar=c(3.5,1.5 + 8 * (max(nchar(colnames(data_ranked)[column_order_to_plot]))/24),1.0,0.6), 
      oma = c(0, 0, 1, 0));
  plot.new();
  plot.window(xlim = c(1.25,xmax-0.5), ylim = c(1,ymax+0.1));
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#99999950",border=NA);
  if(ncol(data_ranked_lonely)) {
    axis(2,at=ymax:1,labels=c(colnames(data_ranked)[column_order_to_plot],"(PROBLEMS LISTED ONCE)"),lwd=1, las = 1, cex.axis = 0.75,hadj =0.96)
  } else {
    axis(2,at=ymax:1,labels=colnames(data_ranked)[column_order_to_plot],lwd=1, las = 1, cex.axis = 0.75,hadj =0.96)
  }
  
  mtext(side = 1,"Count",line=2.5, cex = 1.2);
  abline(v = seq(2, xmax, by = 4), col = "#FFFFFF");
  abline(v = seq(4, xmax, by = 4), col = "#FFFFFF", lty = "44");
  axis(1,las = 1,at = seq(2,xmax,by = 4),cex.axis = 1.2);
  axis(3,las = 1,at = seq(2,xmax,by = 4),labels = NA, cex.axis = 1.2);
  y = ymax + 1;
  for(prob in column_order_to_plot) {
    y = y - 1;
    index = data_ranked %>% pull(prob);
    index[which(index>4)] = 5;
    cumsum_width = 0;
    for(k in unique(sort(index))) {
      width = sum(index == k, na.rm = T);
      points(x = cumsum_width + (1:width),
             y = rep(y,width),
             cex = plot_cex[k],
             pch = plot_pch[k],
             col = plot_colors[k]);
      cumsum_width = cumsum_width + max(width);
    }
  }
  if(ncol(data_ranked_lonely)) {
    y = y - 1;
    cumsum_width = 0;
    for(prob in lonely_column_order_to_plot) {
      index = data_ranked_lonely %>% pull(prob);
      index[which(index>4)] = 5;
      for(k in unique(sort(index))) {
        width = sum(index == k, na.rm = T);
        points(x = cumsum_width + (1:width),
               y = rep(y,width),
               cex = plot_cex[k],
               pch = plot_pch[k],
               col = plot_colors[k]);
        cumsum_width = cumsum_width + max(width);
      }
    }
  }
  legend(x = xmax - 0.5,y = 9,xjust = 1, legend=c("1", "2", "3", "4", ">4"),title=expression(underline(Rank)),ncol=5,pt.cex = max(plot_cex),pch = plot_pch,col=plot_colors,bg="#BBBBBB");
  title(main = title_text, outer = T, adj = 0);
  dev.off();
}
