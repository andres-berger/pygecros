#' @title Complete given initial values by sampling values from bounds 
#'
#' @param init_values A data.frame containing some initial
#' values for the parameters (or NULL)
#' 
#' @param nb_rep Number of repetition of the minimization process
#'
#' @param lb Lower bounds of the parameters
#' 
#' @param ub Upper bounds of the parameters
#' 
#' @param ranseed Set random seed so that each execution give the same results.
#' If you want randomization, set it to NULL
#' 
#' @return A data.frame containing initial values for all parameters and repetitions 
#' of the minimization process that includes the initial values given in `init_values` 
#' plus some randomly chosen ones in bounds `lb` and `ub`
#' 
complete_init_values <- function(init_values, nb_rep, lb, ub, ranseed) {
  
  if (!is.null(lb) && !is.null(ub)) {
    param_names <- names(lb)
    init_values <- CroptimizR:::get_params_init_values(list(lb=lb, ub=ub, init_values=init_values))
    sampled_values <- as.data.frame(CroptimizR:::sample_params(list(lb=lb, ub=ub),nb_rep,ranseed))
    for (param in param_names) {
      idx <- which(!is.na(init_values[,param]))
      sampled_values[idx,param] <- init_values[idx,param]
    }
  } else {
    if (nrow(init_values)<nb_rep || any(is.na(init_values))) {
      stop("Init_values must contain initial values for all parameters and repetitions if ub and lb are not provided.")
    }
    sampled_values <- init_values
  }
  
  return(sampled_values)
}
