#' @title Computes AIC for ordinary least squares
#'
#' @param obs_list List of observed values to use for parameter estimation
#' (same as the one given to `estim_param` function)
#' A `named list` (names = situations names) of data.frame containing
#' one column named Date with the dates (Date or POSIXct format) of the different observations
#' and one column per observed variables with either the measured values or NA, if
#' the variable is not observed at the given date.
#' @param crit_value Final value of the estimated criterion
#' @param param_nb Number of estimated parameters
#'
#' @return Value of the AIC criterion
#'

AIC <- function(obs_list, crit_value, param_nb) {

  # Total number of observations
  n <- sum(sapply(obs_list, function(x)  sum(!is.na(x %>% select(-Date)))))
  
  return(n*log(crit_value/n)+2*param_nb)

}
