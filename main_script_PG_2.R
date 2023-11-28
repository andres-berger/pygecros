# Install / Load needed libraries and functions
library("CroptimizR")
library(reticulate)
library(dplyr)
source("R/AICc.R")
source("R/BIC.R")
source("R/select_param_FwdReg.R")
source("R/main_function_guidelines.R")
source("R/create_AgMIP_outputs.R")
source("R/complete_init_values.R") 
source("Wrappers/PG_wrapper.R")

# Define your model and associated options
# i.e. model_options argument given to estim_param that depends on your model wrapper 

model_function <- PG_wrapper
model_options<-c()  
# Select the observations for the parameter estimation
# i.e. set obs_list here
sit_names <- c("Eradu-2011-TOS1", "Eradu-2011-TOS2", "Eradu-2011-TOS3", 
                "LakeBolac-2011-TOS1", "LakeBolac-2011-TOS2", "LakeBolac-2011-TOS3", 
                "Minnipa-2011-TOS1", "Minnipa-2011-TOS2", "Minnipa-2011-TOS3", 
                "SpringRidge-2011-TOS1", "SpringRidge-2011-TOS2", "SpringRidge-2011-TOS3")

# Select the observations for the parameter estimation
# i.e. set obs_list here
# obs_list must be a named list of data.frames or tibbles (similar to sim_list 
# returned by the model wrapper, see ? estim_param for more information).
# BE CAREFUL: data.frames/tibbles in obs_list MUST ONLY CONTAIN ONE COLUMN PER OBSERVED VARIABLE,
# and ONE COLUMN "Date".The presence of any other column will perturbe the computation of AICc and BIC.
obs_list <- get_obs()
# Give information on the parameters to estimate : 

## Names of obligatory parameters
oblig_param_list <- c("MTDV", "MTDR")

## Names of additional candidate parameters
add_param_list <- c("PSEN","VSEN","VDSA")

param_names <- c(oblig_param_list,add_param_list)

## Bound constraints
param_info_tot <- list(lb=c(MTDV=20 , MTDR=20, PSEN=0. , VSEN=0,   VDSA=0 ),
                       ub=c(MTDV=80 , MTDR=60, PSEN=0.3, VSEN=0.1, VDSA=60 ))

## Lower and upper bounds for initial values sampling 
lb_initV <- c(MTDV=20 , MTDR=20, PSEN=0,   VSEN=0, VDSA=0)
ub_initV <- c(MTDV=80 , MTDR=60, PSEN=0.3, VSEN=0.1, VDSA=60)


# Define parameter estimation algorithm options
optim_options=list(path_results = getwd(), # path where to store the results (graph and Rdata)
                   nb_rep = 4,             # Number of repetitions of the minimization
                   ranseed = 1234)         # set random seed so that each execution give the same results
                                           # If you want randomization, don't set it.
optim_options$maxeval=1000

##########################################

main_function_guidelines(optim_options, oblig_param_list, add_param_list, 
                         param_info_tot, lb_initV, ub_initV, obs_list, 
                         model_function, model_options, info_crit_name="AICc")

main_function_guidelines(optim_options, oblig_param_list, add_param_list, 
                         param_info_tot, lb_initV, ub_initV, obs_list, 
                         model_function, model_options, info_crit_name="BIC")
						 
#print(paste0("Results saved in ",optim_options$path_results))

#load(file.path(…, "optim_results_set***.Rdata"))
#model_results <- model_function(model_options = model_options, param_values = optim_results$final_values), var_names = c(”Zadok65”,”Zadok89")) 
