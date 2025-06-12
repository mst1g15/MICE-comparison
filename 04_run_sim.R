#-------------------------------------------------------------------------------
#This R script contains a runs the simulation on an HPC
#-------------------------------------------------------------------------------

source("00_init.R")
source("01_generate_complete_data.R")
source("01_simulate_missing_data.R")
source("02_impute_functions.R")
source("03_simulation_function.R")

nsim=5000

set.seed(1000)
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(task_id)

all_settings <- readRDS("settings/all_settings.RDS") %>% filter(interact==TRUE) %>% filter(method=="superlearner")
setting_i <- all_settings[task_id,]
res_sim <- run_sim(td=setting_i$true_diff, ss= setting_i$sample_size, miss=setting_i$miss_mech, nsim=nsim, method=setting_i$method, tuning=setting_i$tuning, response=setting_i$response, interact=setting_i$interact)


saveRDS(res_sim, paste0("output/res_seq", task_id, ".RDS"))
