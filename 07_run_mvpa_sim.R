#-----------------------------------------------------------------------------
#this R script runs the simulation based on the inopulse trial on an HPC
#----------------------------------------------------------------------------



source("00_init.R")
source("07_mvpa_functions.R")
source("07_mvpa_simulation_function.R")


nsim=5000

# Get number of cores from SLURM
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Register parallel backend
cl <- makeCluster(cores)
registerDoParallel(cl)

cat("Registered", getDoParWorkers(), "workers\n")





set.seed(1000)
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(task_id)


all_settings <- readRDS("settings/MVPA_settings.RDS")
setting_i <- all_settings[task_id,]
res_sim <- run_sim_mvpa(method=setting_i$method, miss=setting_i$miss, nsim=nsim)

stopCluster(cl)


saveRDS(res_sim[[1]], paste0("output_MVPA/ANOVA", task_id, ".RDS"))
saveRDS(res_sim[[2]], paste0("output_MVPA/ANOVA_means", task_id, ".RDS"))

saveRDS(res_sim[[3]], paste0("output_MVPA/marginal", task_id, ".RDS"))
saveRDS(res_sim[[4]], paste0("output_MVPA/collapsed", task_id, ".RDS"))
saveRDS(res_sim[[5]], paste0("output_MVPA/MMRM_means", task_id, ".RDS"))

