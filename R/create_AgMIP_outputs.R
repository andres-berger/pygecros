#' @title Fill a data.frame with the optim_resultsults required by the AgMIP calibration PhaseIII protocol
#'
#' @param candidate_params Vector of names of current candidate parameters 
#'
#' @param obs_list List of observed values to use for parameter estimation
#' 
#' @param optim_results Results of the parameter estimation, as given by estim_param, 
#' for the current candidate parameters
#' 
#' @param model_function Crop Model wrapper function to use.
#' 
#' @param model_options List of options for the Crop Model wrapper (see help of
#' the Crop Model wrapper function used).
#' 
#' @param digits Number of digits to take into account for outputs printing format
#' 
#' @param info_crit_func Function that computes the information criterion (AICc or BIC)
#' 
#' @param info_crit_name Name of the information criterion to use for parameter selection ("AICc" or "BIC")
#' 
#' @return A data.frame containing for each set of candidate parameters, the names of the parameters 
#' and the initial and final values of the parameters, the OLS criterion and of the information criterion.
#' 

create_AgMIP_outputs <- function(candidate_params, obs_list, optim_results, model_function, 
                                 model_options, info_crit_func, info_crit_name, digits) {
  
  # Compute initial value of criterion and information criterion (for the "best" repetition)
  model_results <- model_function(model_options = model_options, 
                                  param_values = optim_results$init_values[optim_results$ind_min_crit,], 
                                  sit_names = names(obs_list),
                                  sit_var_dates_mask = obs_list)    
  
  obs_sim_list <- CroptimizR:::intersect_sim_obs(model_results$sim_list, obs_list)
  init_crit <- crit_ols(obs_sim_list$sim_list, obs_sim_list$obs_list)
  init_info_crit <- info_crit_func(obs_sim_list$obs_list, init_crit, param_nb=length(candidate_params))
  final_info_crit <- info_crit_func(obs_list, optim_results$min_crit_value, param_nb=length(candidate_params))
  
  # Gather optim_results
  df <- setNames(data.frame(paste(candidate_params,collapse = ", "), 
                            paste(format(optim_results$init_values[optim_results$ind_min_crit,], 
                                         scientific=FALSE, digits=digits, nsmall=2), collapse=", "), 
                            #ini_crit=format(init_crit, scientific=FALSE, digits=digits, nsmall=2), 
                            #init_info_crit=format(init_info_crit, scientific=FALSE, digits=digits, nsmall=2), 
                            paste(format(optim_results$final_values, scientific=FALSE, digits=digits, nsmall=2),collapse=", "), 
                            format(optim_results$min_crit_value, scientific=FALSE, digits=digits, nsmall=2), 
                            format(final_info_crit, scientific=FALSE, digits=digits, nsmall=2)),
                 c("Estimated parameters","Initial parameter values","Final values",
                   "Sum of squared errors",info_crit_name))
  return(df)
  
}
