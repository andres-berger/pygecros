#' @title Apply the AgMIP phase III protocol (Crop model calibration using phenology data)
#'
#' @param optim_options List of options of the parameter estimation method:
#'    o `path_results` The path where to store the results (optional, default=getwd())
#'    o `nb_rep` The number of repetition of the minimization
#'    o `ranseed`  Set random seed so that each execution give the same results
#       If you want randomization, don't set it.
#' 
#' @param oblig_param_list Vector of names of parameters that must be estimated 
#' (list of obligatory parameters)
#'
#' @param add_param_list Vector of names of additional parameters candidate to the 
#' estimation
#' 
#' @param param_info_tot Information on all the candidate parameters to estimate.
#' A list containing:
#'    - `ub` and `lb`: named vectors of upper and lower bounds, THAT WILL BE USED 
#'    FOR CONSTRAINING THE RANGE OF PARAMETERS DURING MINIMIZATION PROCESS. 
#'    Bounds must be defined FOR ALL parameters, set -Inf and Inf if you don't need 
#'    to constrain the bounds.
#'    - `init_values`, A data.frame containing initial values for the parameters.
#'    One column per parameter, one line per repetition of the minimization.
#'    Set it to NULL if you want all the initial values to be automatically sampled 
#'    (within bounds `lb_initV` and `ub_initV`).
#'    If you want to provide initial values for only a subpart of the parameters or repetitions, 
#'    set NA for parameters and/or repetitions for which you do not provide values.
#'
#' @param lb_initV Named vector of parameters' lower bounds to use for sampling initial values
#' 
#' @param ub_initV Named vector of parameters' upper bounds to use for sampling initial values
#' 
#' @param obs_list List of observed values to use for parameter estimation
#' A `named list` (names = situations names) of data.frame containing
#' one column named Date with the dates (Date or POSIXct format) of the different observations
#' and one column per observed variables with either the measured values or NA, if
#' the variable is not observed at the given date.
#' 
#' @param model_function Crop Model wrapper function to use.
#' 
#' @param model_options List of options for the Crop Model wrapper (see help of
#' the Crop Model wrapper function used).
#' 
#' @param digits Number of digits to take into account for outputs printing format
#' 
#' @param info_crit_name Name of the information criterion to use for parameter selection ("AICc" or "BIC")
#' 
#' @return A data.frame containing for each set of candidate parameters, the names of the parameters 
#' and the initial and final values of the parameters, the OLS criterion and of the information criterion used.
#' 
#' This data.frame is also saved in csv and Rdata formats (AgMIP_outputs.* files, 
#' saved in optim_options$path_results)

main_function_guidelines <- function(optim_options, oblig_param_list, add_param_list, 
                                     param_info_tot, lb_initV=NULL, ub_initV=NULL, 
                                     obs_list, model_function, model_options,
                                     info_crit_name, digits=3) {
  
  # Checks
  if (is.null(param_info_tot$lb) || is.null(param_info_tot$ub)) {
    stop("param_info_tot$lb and param_info_tot$ub must be defined!")
  } else if (any(sort(names(param_info_tot$lb))!=sort(names(param_info_tot$ub))) ||
             any(sort(names(param_info_tot$lb))!=sort(c(oblig_param_list, add_param_list)))) {
    stop("param_info_tot$lb and param_info_tot$ub must be defined for all candidate parameters (i.e. obligatory and additional ones)")
  }
  if (!is.null(param_info_tot$init_values)) {
    if (any(sort(names(param_info_tot$init_values))!=sort(c(oblig_param_list, add_param_list)))) {
      stop("If param_info_tot$init_values is provided, it must have a column for each candidate parameter 
      (i.e. obligatory and additional ones). Fill this column with NA if you don't want to define initial values for this parameter.")
    }
  }
  if ((is.null(lb_initV) || is.null(ub_initV)) && 
      (nrow(param_info_tot$init_values)<model_options$nb_rep || any(is.na(param_info_tot$init_values)))) {
    stop("Init_values must contain initial values for all parameters and repetitions if lb_initV and ub_initV are not provided.")
  }
  if (info_crit_name=="AICc") {
    info_crit_func <- AICc
  } else if (info_crit_name=="BIC") {
    info_crit_func <- BIC
  } else {
    stop("info_crit_name must be AICc or BIC.")
  }
  
  
  param_info_tot$init_values <- complete_init_values(param_info_tot$init_values, 
                                                     optim_options$nb_rep, 
                                                     lb_initV, ub_initV, 
                                                     optim_options$ranseed)
  
  candidate_params <- oblig_param_list
  best_final_values <- setNames(data.frame(matrix(data=NA, ncol=length(candidate_params),nrow=0)),
                                candidate_params)
  prev_info_crit <- Inf
  count<-1
  df_outputs<-NULL
  
  while(!is.null(candidate_params)) {
    
    print(paste("Estimated parameters:",paste(candidate_params,collapse=" ")))
    
    # initialize already estimated parameters with the values leading to the best criterion obtained so far
    param_info <- lapply(param_info_tot,function(x) x[candidate_params])
    if (nrow(best_final_values)>0 && length(best_final_values)>0) {
#      param_info$init_values <- setNames(data.frame(matrix(data=NA, ncol=length(candidate_params),
#                                                           nrow=optim_options$nb_rep)),
#                                         candidate_params)
      param_info$init_values[names(best_final_values)] <- best_final_values[rep(1,optim_options$nb_rep),]
      new_candidate <- candidate_params[length(candidate_params)]
      param_info$init_values[, new_candidate] <- param_info_tot$init_values[,new_candidate]
    }
    
    optim_results <-     estim_param(obs_list=obs_list,
                                     model_function=model_function,
                                     model_options=model_options,
                                     optim_options=optim_options,
                                     crit_function=crit_ols,
                                     param_info=param_info)
    
    
    optim_results$info_crit <- info_crit_func(obs_list, optim_results$min_crit_value, 
                             param_nb=length(candidate_params))
    if (optim_results$info_crit<min(prev_info_crit)) {
      best_final_values <- tibble::tibble(!!!optim_results$final_values)
    }
    
    # Save and move results
    save(optim_results, file = file.path(optim_options$path_results,
                                         paste0("optim_results_set",count,".Rdata")))
    file.copy(from="EstimatedVSinit.pdf",
              to=paste0("EstimatedVSinit_set",count,".pdf"), overwrite = TRUE)
    
    print(paste("Values for the estimated parameters:",paste(optim_results$final_values,
                collapse=" ")))
    print(paste0(info_crit_name,"=",optim_results$info_crit))
    
    df_outputs<-bind_rows(df_outputs,create_AgMIP_outputs(candidate_params, obs_list, 
                                                          optim_results, model_function, 
                                                          model_options, info_crit_func, 
                                                          info_crit_name, digits))
    
    candidate_params <- select_param_FwdReg(oblig_param_list, add_param_list, 
                                            candidate_params, optim_results$info_crit, 
                                            prev_info_crit)
    prev_info_crit <- c(prev_info_crit, optim_results$info_crit)
    
    count <- count+1
  }
  
  # Save the results
  save(df_outputs, file = file.path(optim_options$path_results,
                                    paste0("phase3_",info_crit_name,"_Table4.Rdata")))
  
  df <- as.data.frame(sapply(df_outputs, as.character, simplify = FALSE))
  write.table(df, sep=";", file=file.path(optim_options$path_results,
                                          paste0("phase3_",info_crit_name,"_Table4.csv")), row.names=FALSE)
 
  return(df) 
}
