# This script contains functions related to the calculation of the LDRBO
# metric proposed in Krauss, et al (2015). 


seq_size_intersection = function(x, y, max_d = NULL) {
  if(is.null(max_d)) {
    max_d = max(min(length(x), which(is.na(x))),
                min(length(y), which(is.na(y))));
  }
  seq_max_d = seq_len(max_d);
  agreement = numeric(max_d);
  agreement[1] = x[1] == y[1];
  for(d in seq_max_d[-1]) {
    agreement[d] =
      (!is.na(x[d])) * (x[d] %in% y[seq_len(d)]) + 
      (!is.na(y[d])) * (y[d] %in% x[seq_len(d-1)])
  }
  cumsum(agreement);
}

ldrbo = function(dat_new, 
                 psi = 1, 
                 dat_ref = dat_new, 
                 verbose_results = TRUE) {
  
  if(class(dat_new)!="matrix" || class(dat_ref)!="matrix") {
    stop("'dat_new' and 'dat_ref' both must be matrices");
  }
  while(ncol(dat_new) < ncol(dat_ref)) {dat_new = cbind(dat_new,NA);}
  while(ncol(dat_ref) < ncol(dat_new)) {dat_ref = cbind(dat_ref,NA);}
  max_depth = ncol(dat_new);
  n_new_samps = nrow(dat_new);
  n_ref_samps = nrow(dat_ref);
  
  # psi2d is an n*n*max_depth array, where each n*n submatrix has psi^d in 
  # every element with d = 1,...,max_depth
  psi2d = array(matrix(psi ^ (1:max_depth),
                       nrow = n_new_samps * n_ref_samps,
                       ncol = max_depth, byrow = T),
                c(n_new_samps, n_ref_samps, max_depth));  
  
  # agreement is an n*n*max_depth array, calculating the size of each pairwise 
  # problem list intersection for increasing depths this is a modification of 
  # A_d (equation 3) in Webber, et al. that takes into account the length
  # of the lists as well. 
  # agreement at d = [overlap(list1,list) + (d-max(length(list1),list2))_+]/d
  agreement = array(0,c(n_new_samps, n_ref_samps, max_depth));
  for(i in 1:n_new_samps) {
    for(j in 1:n_ref_samps) {
      max_d = max(sum(!is.na(dat_new[i,])),
                  sum(!is.na(dat_ref[j,])));
      
      agreement[i,j,1:max_d] = 
        seq_size_intersection(x = dat_new[i,], 
                              y = dat_ref[j, ], 
                              max_d = max_d) / (1:max_d);
      
      psi2d[i,j,-(1:max_d)] = 0;
    }
  }
  
  rbo = rowSums(agreement * psi2d, dims = 2) / rowSums(psi2d, dims = 2);
  
  if(!is.null(rownames(dat_new))) {
    rownames(rbo) = rownames(dat_new);
  }
  if(!is.null(rownames(dat_ref))) {
    colnames(rbo) = rownames(dat_ref);
  }
  
  if(verbose_results) {
    list(rbo = rbo, 
         psi = psi, 
         agreement = agreement, 
         psi2d = psi2d);
  } else {
    rbo;
  }
  
}

# function name and purpose: 'consensus_ldrbo' estimates the hypothetical o
# ordered list that maximizes the median or mean LDRBO across all input 
# ordered lists. 
#
# author: Phil Boonstra (philb@umich.edu)
#
# date: 10/17/2019, 10:00am EST
# 
#
# Function inputs:
#
### dat: (matrix, integer valued) each row of the matrix corresopnds to an 
# ordered list. Thus, no integers should appear more than once in a row. Also, 
# supposing there are v unique items that were ranked at least once by at 
# least one ranker, then the item labels should be consecutive integers 
# starting at 1 and ending at v. Missing data are allowed but should only 
# appear at the end of a list, i.e. if the s-th element of a row is NA, then 
# all elements beyond the s-th element should also be NA. 
# NOTE: play close attention to the definition of an ordered list. The s-th 
# entry of an ordered list is the integer label of the item that is ranked 
# s-th (i.e. items appearing early in the list are ranked higher).
# 
### psi: (numeric, positive) the weight parameter that controls the top-
# weightedness of the function. Smaller values of psi give more importance
# to agreement at higher ranks 
#
### max_size: (integer, positive) an upper bound on the consensus
# problem list to consider. Smaller values will make the algorithm run faster
# at the cost of potentially missing out a longer consensus list that improves
# overall agreement. If missing, this defaults to the number of unique items
# observed, i.e. the maximum possible. NOTE: this doesn't mean that the 
# returned consensus list will always have length 'max_size'; rather, this value
# excludes from consideration any list that is larger than 'max_size'
#
### min_size: (integer, positive but no larger than 'max_size') a lower bound
# on the consensus problem list to consider. A consensus list that has length
# less than this value will never be the returned consensus list, even if it
# optimizes the objective. Smaller values of 'min_size' will *not* make 
# the algorithm run faster. Generally, you can leave this at its default 
# value of 1; however, another logical choice would be to set it equal to your
# choice of 'max_size' if you want the best consensus list subject to having
# length exactly equal to 'min_size' (which equals 'max_size')
#
### look_beyond_init, look_beyond_final: (integer, positive) Tuning parameters
# that control the greediness of the algorithm. An equally spaced sequence with 
# length equal to the maximum possible consensus list size ('max_size') is 
# formed using these parameters as the start and and end points. At each iteration, 
# the value from this sequence is the maximum number of suboptimal candidate 
# lists that will be retained for further exploration (if there are fewer candidate 
# lists than this value, then they will all be retained) Typically 'look_beyond_init'
# is greater that 'look_beyond_final' based on the assumption that it is more 
# important to consider currently suboptimal lists at the beginning of a 
# consensus list than the end. Smaller values make the algorithm run faster 
# by being greedier, i.e. focusing on lists that have already proven to
# have good overall consensus at the expense of lists that have the potential to
# have good overall consensus. Larger values make the algorithm run slower by 
# considering a larger number of possible lists that are not currently optimal 
# but may turn out to be optimal. 
# 
### window_seq: (vector, positive integers) a vector of integers with length 
# equal to max_size. At the kth step, the algorithm will look at the next 
# window_seq[k] elements from 'initial_order' (see next entry) for possible 
# inclusion. Thus, larger numbers mean that the algorithm considers more 
# possibilities. Together, the choice of 'look_beyond' and 'window_seq' control 
# what might be considered the "greediness-slowness tradeoff" of the algorithm.
# If 'window_seq' is not provided, the algorithm constructs its own default value. 
#
### initial_order: (vector, positive integers) a permutation of the vector of
# unique item labels that have been ranked. Not all items can be considered 
# at all steps of the algorithm always due to time and memory constraints, thus 
# 'initial_order', in conjunction with 'window_seq', prioritizes the items in 
# terms of when they should be considered. 
#
### objective: (character equal to "median" or "mean"): should the median or 
# mean pairwise LDRBO be maximized? Anything besides 'median' is assumed to 
# be 'mean'
#
### Value
#
# A list with the following named objects:
# 
# consensus_list: An ordered list that the algorithm identified to have the 
# maximum median or mean pairwise LDRBO with each ordered list in dat
#
# consensus_total_rbo: the maximized median or mean pairwise LDRBO
# 
# consensus_tie: a matrix of alternative consensus lists that had the same
# maximized pairwise LDRBO. 
# 
# psi: the value of psi that was used by the algorithm
# 
# candidates: a matrix of candidates that were close-ish to the
# maximum. These include the best consensus list that was found at each iteration, 
# thus you will see at least one row for each length between 1 and whenever
# the function stopped. Rows are ordered with respect to their objective function.
# 
# candidates_total_rbo: a vector with length equal to number of rows of 'candidates'
# that gives the median or mean pairwise LDRBO betweeen that candidate and each
# ordered list in the data. Elements are sorted in decreasing order and correspond
# to the respective row of 'candidates'
#
# Details: 
# The kth step of the iteration, starting with k = 2, has three substeps. First, 
# grow all lists currently under consideration, each of which is exactly length 
# k-1, by exactly one item. For a given list under consideration, choose the 
# first 'window_seq[k]' items from the vector 'initial_order' that have not yet 
# been ranked in that list, and create'window_seq[k]' new lists from that list. 
# Thus, if the number of lists under consideration at the start of this step was
# m, then the number of lists under consideration at th end of this step is 
# m * window_seq[k], and each of these new lists is now length-(k). Second, 
# calculate the values x(obs) and x(max) for each of the m * window_seq[k] 
# lists x under consideration. This involves comparing each list to the observed 
# data, evaluating each pairwise LDRBO value, and calculating the mean or median 
# of all pairwise LDRBO values. Third, order all lists according to the criterion 
# above, dropping all lists for which the criterion equals 0. If all lists have
# a stricitly positive value of the criterion, then they are all candidates. 
# However, we can only keep at most 'look_beyond[k]' candidates for the next 
# iteration, and so we drop the lists with the smallest values of the criterion
# until we have only 'look_beyond[k]' candidates remaining. Then, 
# increment k, and go back to substep 1. Stop when k == max_size and return 
# list y, which achieves a consensus of y(obs) and is, by definition, the largest
# LDRBO we observed
# 
# This function is not optimized for speed or efficiency.  There are a few 
# calculations and loops that can get very expensive both in  terms of memory 
# and time. Specifically, a lower bound on the memory required for this 
# function (in bytes) is 40 * x * y * z, where
# x = max(c(0, look_beyond) * c(window_seq, 0))
# y = nrow(dat);
# z = max(max_size, ncol(dat));
# 
# Suggestions on improving the speed or efficiency of this algorithm are always
# welcomed via github pull requests. 

consensus_ldrbo = function(dat,
                           psi,
                           max_size,
                           min_size = 1,
                           look_beyond_init = 1e3,
                           look_beyond_final = 1e3,
                           window_seq = NULL,
                           initial_order = NULL,
                           objective = "median", 
                           verbose = FALSE)  {
  
  require(gtools);# we require gtools::permutations()
  num_obs = nrow(dat);
  num_uniq = sum(!is.na(unique(as.numeric(dat))));
  tiny_positive = sqrt(.Machine$double.eps);
  if(missing(max_size)) {max_size = num_uniq;}
  longest_list = max(max_size, ncol(dat));
  all_list_lengths = rowSums(!is.na(dat));
  
  seq_longest_list = seq_len(longest_list);
  
  #Check for invalid / problematic arguments
  if(psi <= 0) {stop("'psi' should be positive");}
  if(min_size > num_uniq || min_size < 1) {
    stop("'min_size' should be less than or equal to the number of
         unique items and greater than or equal to 1");
  }
  if(max_size > num_uniq || max_size < min_size) {
    stop("'max_size' should be less than or equal to the number of
         unique items and greater than or equal to 'min_size'");}
  if(num_uniq != max(dat,na.rm = TRUE)) {
    stop("Item labels in 'dat' must be consecutively labeled 
         integers starting with 1");
  }
  if(!isTRUE(all(apply(apply(is.na(cbind(dat,NA)),1,cumsum),2,diff) | 
                 (apply(is.na(dat),1,cumsum) == 0))) | 
     isTRUE(any(apply(dat,1,duplicated) & !apply(dat,1,is.na)))) {
    stop("'dat' must contain uninterrupted left-to-right sets of unique integers, i.e. no NAs between items");
  }
  if(max(colSums(is.na(dat))) == num_obs) {
    stop("One or more columns of 'dat' contain only NAs");
  }
  if(any(apply(dat,1,duplicated,incomparables=NA))) {stop("At least one list has multiple ranks for one item")} 
  
  
  if(min(look_beyond_init,look_beyond_final) < 1) {
    stop("'look_beyond_init' and 'look_beyond_final' must be 
         positive integers");
  }
  
  # The first value doesn't matter because the number of one-item candidates
  # to keep is determined by window_seq[1];
  look_beyond = c(0, round(seq(from = look_beyond_init, 
                               to = look_beyond_final, 
                               length = max_size - 1)));
  
  # Initial ranking is a simplistic point system: each item in each list gets
  # a certain number of points (between 0.04 and 25) for each occurrence in a
  # list, with higher rankings getting it more points
  
  points_seq = exp(seq(from = log(25), to = log(1/25), length = ncol(dat)));
  rank_points = numeric(num_uniq);
  for(i in seq_len(num_uniq)) {
    rank_points[i] = 
      sum(colSums(dat == i, na.rm = T) * points_seq)
  }
  if(is.null(initial_order)) {
    initial_order = order(rank_points, decreasing = TRUE);
  }
  if(!isTRUE(all.equal(sort(initial_order), seq_len(num_uniq)))) {
    stop("'initial_order' be a permutation of the integers 1 to 'num_uniq', 
         where 'num_uniq' is the number of unique objects ranked");
  }
  
  rm(i, points_seq);
  
  # If not provided, window_seq is an increasing-then-decreasing vector. 
  # the first value is the larger of 5 and the number of items that have 
  # at least 1/2 as many rank_points as the item with the most rank points. 
  # Then window_seq increases up to 'num_uniq' and finally decreases down 
  # to 1 at the end. Note that, there is a natural upper bound on what 
  # window_seq can take on: at the k-th iteration, i.e. after the lists
  # have been grown to length k, then there are only num_uniq - k possible 
  # items left to consider. 
  
  if(is.null(window_seq)) {
    
    window_init = 
      pmin(num_uniq, pmax(5, sum(rank_points > 0.5 * max(rank_points))));
    window_seq = 
      pmin(1 + num_uniq - (1:max_size),
           ceiling(seq(from = window_init, to = num_uniq, length = max_size)));
    rm(window_init);
  }
  if(length(window_seq) != max_size) {
    stop("'window_seq' should have length equal to 'max_size'");
  }
  
  # Initialize big matrices to try to save on expensive calculations
  
  max_rows = max(c(0, look_beyond) * c(window_seq, 0), 
                 prod(window_seq[1:2]))
  
  if((gB_required <- 5 * max_rows * num_obs * longest_list * 8 / (1024^3)) > 0.5 &&
     verbose) {
    cat(formatC(gB_required,format="f", digits = 3),
        "gigabytes of memory required by 'consensus_ldrbo'\n");
  }
  
  baseline_agreement = 
    array(0,c(max_rows,
              num_obs,
              longest_list));
  
  baseline_psi2d = 
    array(matrix(psi ^ seq_longest_list,
                 nrow = num_obs * max_rows,
                 ncol = longest_list,
                 byrow = T),
          c(max_rows,
            num_obs,
            longest_list));  
  
  old_search = 
    matrix(initial_order[seq_len(window_seq[1])], ncol = 1);
  
  if(objective == "median") {
    curr_search_total_rbo = apply(
      ldrbo(dat_new = old_search,
            psi = psi, 
            dat_ref = dat,
            verbose_results = FALSE),1,median);
  } else {
    curr_search_total_rbo = 
      rowMeans(ldrbo(dat_new = old_search,
                     psi = psi, 
                     dat_ref = dat,
                     verbose_results = FALSE));
  }
  top_candidates = matrix(which.max(curr_search_total_rbo), nrow = 1)
  
  k=1;
  
  while(k < max_size) {
    k = k + 1; 
    if(verbose) {cat("k = ", k, "\n");}
    curr_window = window_seq[k];
    curr_search = NULL;
    seq_extend = seq_len(min(num_uniq, k - 1 + curr_window));
    n_permutations = min(num_uniq - k + 1, curr_window);
    seq_n_permutations = seq_len(n_permutations);
    curr_search = 
      cbind(old_search[rep(seq_len(nrow(old_search)), each = n_permutations), , drop = FALSE], NA);
    for(j in seq_len(nrow(old_search))) {
      curr_index = seq_n_permutations + (j - 1) * n_permutations;
      
      curr_search[curr_index, k] = 
        gtools::permutations(
          n = n_permutations,
          r = 1,
          v = setdiff(initial_order[seq_extend],
                      old_search[j,]));
    }
    rm(j, curr_index);
    if(verbose) {cat("nrow(curr_search) = ", nrow(curr_search), "\n\n");}
    
    obs_agreement = 
      max_agreement = 
      baseline_agreement[1:nrow(curr_search),,,drop=FALSE];
    
    psi2d = 
      baseline_psi2d[1:nrow(curr_search),,,drop=FALSE];
    
    seq_nrow_curr_search = seq_len(nrow(curr_search));
    for(j in seq_len(num_obs)) {
      length_j = all_list_lengths[j];
      # max_d is the longer of the current consensus list and the current
      # ranker's list
      max_d = max(k, length_j);
      min_max_size_max_d = min(max_size, max_d);
      seq_max_d = seq_len(max_d);
      
      # Can only measure agreement up to the longer of the two lists
      psi2d[,j,-seq_max_d] = 0;
      for(i in seq_nrow_curr_search) {
        
        
        obs_agreement[i,j,1:max_d] =
          max_agreement[i,j,1:max_d] = 
          seq_size_intersection(x = curr_search[i,], 
                                y = dat[j, ], 
                                max_d = max_d) / (1:max_d);
        
        # Extrapolate to the maximum possible agreement
        if(k < min_max_size_max_d) {
          max_agreement[i,j,(k+1):min_max_size_max_d] = 
            k * (obs_agreement[i,j,k]-1) / ((k+1):min_max_size_max_d) + 1;
        }
        
        # But beyond maximum allowed size, agreement must go down
        if(max_d < max_size) {
          max_agreement[i,j,(max_d+1):max_size] = 
            max_d * (max_agreement[i,j,max_d]) / ((max_d+1):max_size)
        }
      }
    }
    rm(i, j);
    
    rowSums_psi2d = rowSums(psi2d, dims = 2);
    curr_search_rbo = rowSums(obs_agreement * psi2d, dims = 2) / rowSums_psi2d;
    
    if(objective == "median") {
      curr_search_total_rbo = apply(curr_search_rbo,1,median);
      max_possible_rbo = apply(rowSums(max_agreement * psi2d, dims = 2) / rowSums_psi2d,1,median);
    } else {
      curr_search_total_rbo = rowMeans(curr_search_rbo);
      max_possible_rbo = rowMeans(rowSums(max_agreement * psi2d, dims = 2) / rowSums_psi2d);
    }
    top_candidates = 
      rbind(cbind(top_candidates,NA),
            curr_search[which((max(curr_search_total_rbo) - curr_search_total_rbo) < tiny_positive),])
    
    # This is the step that looks beyond the current consensus list size
    # to determine the branches that can be pruned because it would be 
    # impossible for them to exceed the current achieved level of consensus
    for(extend_d in which(colSums(psi2d[1,,] == 0) != 0)) {
      psi2d[,,extend_d] = psi ^ extend_d;
      if(objective == "median") {
        max_possible_rbo = pmax(max_possible_rbo,apply(rowSums(max_agreement * psi2d, dims = 2) / rowSums(psi2d, dims = 2),1,median));
      } else {
        max_possible_rbo = pmax(max_possible_rbo,rowMeans(rowSums(max_agreement * psi2d, dims = 2) / rowSums(psi2d, dims = 2)));
      }
    }
    rm(extend_d);
    
    if(all(max_possible_rbo <= max(curr_search_total_rbo) + 1e-4) &&
       k >= min_size) {
      break;
    }
    
    
    # To order the candidate methods in the current search
    # (i) evalute the additional rbo that each row must gain to achieve its 
    # full potential relative to what the current best observed rbo would need 
    # to achieve to get that same maximum. 
    # (ii) exclude any rows that cannot exceed the current rbo value that has 
    # already been achieved by a list (this is a fairly weak exclusion)
    
    criterion = 
      ((max_possible_rbo - max(curr_search_total_rbo)) / 
         # The small number below is arbitrary: we just need a strictly positive
         # denominator for the cases that the numerator is negative, so that the
         # total value is a large but finite negative number instead of -Inf. 
         # Otherwise, we will end up with -Inf * 0 
         pmax(1e-8,(max_possible_rbo - curr_search_total_rbo))) * 
      ((max_possible_rbo > max(curr_search_total_rbo) + 1e-4) | 
         (near(curr_search_total_rbo, max(curr_search_total_rbo))));
    
    if(look_beyond[k] < sum(criterion > 0)) {
      
      order_criterion = 
        order(criterion, decreasing = TRUE)[1:look_beyond[k]];
      
    } else {
      
      order_criterion = 
        order(criterion, decreasing = TRUE)[1:sum(criterion > 0)];
      
    }
    
    old_search = 
      curr_search[order_criterion,,drop = FALSE];
    
    rm(curr_search, criterion, order_criterion, max_possible_rbo, curr_search_total_rbo);
  }
  
  curr_search = 
    top_candidates[which(rowSums(!is.na(top_candidates)) >= min_size),,drop = F];
  curr_search_rbo = ldrbo(dat_new = curr_search,
                          psi = psi,
                          dat_ref = dat);
  
  if(objective=="median") {
    curr_search_total_rbo = apply(curr_search_rbo$rbo,1,median);
  } else {
    curr_search_total_rbo = rowMeans(curr_search_rbo$rbo);
  }
  
  curr_search = 
    curr_search[order(curr_search_total_rbo, decreasing = TRUE),,drop = FALSE];
  curr_search_total_rbo = sort(curr_search_total_rbo, decreasing = TRUE);
  
  candidate_ranking = order(-curr_search_total_rbo);
  consensus_tie = which(max(curr_search_total_rbo)-curr_search_total_rbo < tiny_positive);
  
  #This is just to reflect that all lists at step 1 are always retained
  look_beyond[1] = Inf; 
  
  list(consensus_list = curr_search[1,1:sum(!is.na(curr_search[1,]))],
       consensus_total_rbo = curr_search_total_rbo[1], 
       consensus_tie = curr_search[which(curr_search_total_rbo[1] <= curr_search_total_rbo + 0.001),,drop = FALSE],
       control = list(dat = dat, 
                      psi = psi, 
                      initial_order = initial_order, 
                      min_size = min_size, 
                      max_size = max_size,
                      window_seq = window_seq, 
                      look_beyond = look_beyond, 
                      objective = objective),
       candidates = curr_search, 
       candidates_total_rbo = curr_search_total_rbo);
}




minimum_ldrbo <- function(size_library, 
                          length_list1 = length(library_size),
                          length_list2 = length(library_size), 
                          psi = 1) {
  
  if(any(pmax(length_list1, length_list2) > size_library) ||
     any(pmin(length_list1, length_list2) < 1)) {
    stop("'length_list1' and 'length_list2' must both be integers
         in [1, size_library]")
  }
  
  list1 = matrix(seq(1, length_list1, by = 1), nrow = 1);
  list2 = matrix(NA, nrow = length(length_list2), ncol = max(length_list2))
  for(i in seq_along(length_list2)) {
    list2[i,seq_len(length_list2[i])] = seq(size_library, size_library - length_list2[i] + 1, by = -1);
  }
  
  tibble(
    size_library = size_library,
    length_list1 = length_list1, 
    length_list2 = length_list2, 
    min_possible_ldrbo = drop(ldrbo(dat_new = list1,
                                    psi = psi,
                                    dat_ref = list2,
                                    verbose = FALSE)))
}


