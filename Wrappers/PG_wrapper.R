#' @title My model wrapper for CroptimizR
#'
#' @description This function runs my crop model on a set of situations (i.e. environments) 
#' using the values of the parameters defined in the param_values argument. It returns 
#' the values of the simulated outputs.
#'
#' @param param_values (optional) a named vector that contains the value(s) and name(s) 
#' of the parameters to force for each situation to simulate. If not provided (or if is 
#' NULL), the simulations will be performed using default values of the parameters 
#' (e.g. as read in the model input files).
#'
#' @param sit_names Vector of situations names for which results must be returned. 
#'
#' @param model_options List containing any information needed to run the model 
#' (e.g. path to model input files and executable, ...) 
#'
#' @return A list containing:
#'     o `sim_list`: a `named list` (names = situations names) of data.frames (or tibbles) of 
#` simulated output values (one element of the list per simulated situation) 
#'     o `error`: an error code indicating if at least one simulation ended with an error.

#setup model
source_python("cal.py")
s<-sim()
get_obs <- function() {
  obs_list <- sapply(s$obs(), function(x)  py_to_r(x) , simplify = FALSE)
  obs_list <- sapply(obs_list, function(x) data.frame(Date=as.POSIXct(x$Date),Zadok65=x$Zadok65,Zadok89=x$Zadok89), simplify = FALSE)
  return(obs_list)
}

PG_wrapper <- function( param_values=NULL, sit_names, model_options, ...) {
  #print(paste0("parameters ",param_values))
  d<-s$run_sites(sit_names,param_values)
  d <- sapply(d, function(x)  py_to_r(x) , simplify = FALSE)
  d <- sapply(d, function(x) data.frame(Date=as.POSIXct(x$Date),Zadok65=x$Zadok65,Zadok89=x$Zadok89), simplify = FALSE)
  #print(paste0("Model wrapper otput ",d))
  results <- list(sim_list = d, error=FALSE)
  return(results)
}


