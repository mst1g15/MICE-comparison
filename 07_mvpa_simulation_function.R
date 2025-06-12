#-----------------------------------------------------------------------------
#this R script contains a function to run the simulation based on the 
#inopulse trial. This is set up to run with parallelisation on the HPC
#----------------------------------------------------------------------------



safe_extract <- function(x, name) {
  if (is.list(x) && !is.null(x[[name]])) {
    return(x[[name]])
  } else {
    return(NULL)
  }
}


run_sim_mvpa <- function(method, miss, nsim){
  
  
  results <-  foreach(i=1:nsim) %dopar% {
    
    tryCatch({
    source("00_init.R")
    source("07_mvpa_functions.R")
    source("07_mvpa_simulation_function.R")
    
    
  dat <- generate_daily_data(ss=145, weekno=4, rand_ratio=NULL)
  dat_long = dat[[1]]
  dat_wide = dat[[2]]
  
  dat_miss <- induce_missingness(dat_wide, miss, miss_prop_n=0.3)
  
  results <- impute_and_analyse(dat_miss, method)
  
  ANOVA_results <- results[[1]]
  ANOVA_means <- results[[2]]
  marginal_results <- results[[3]]
  collapsed_results <- results[[4]]
  MMRM_means <- results[[5]]
  
  ANOVA_results$miss=miss
  ANOVA_results$sim=i
  ANOVA_results$analysis="ANOVA"
  
  ANOVA_means$miss=miss
  ANOVA_means$sim=i
  ANOVA_means$analysis="ANOVA"
  
  
  marginal_results$miss=miss
  marginal_results$sim=i
  marginal_results$analysis="marginal"
  
  collapsed_results$miss=miss
  collapsed_results$sim=i
  collapsed_results$analysis="collapsed"
  
  MMRM_means$miss=miss
  MMRM_means$sim=i
  MMRM_means$analysis="MMRM"

  list(ANOVA_results=ANOVA_results, ANOVA_means=ANOVA_means, 
       marginal_results=marginal_results, collapsed_results=collapsed_results, MMRM_means=MMRM_means)
    }, error = function(e) {
      message(sprintf("Iteration %d failed: %s", i, e$message))
      return(NULL)
    })
  }
  
  
  results <- Filter(Negate(is.null), results)
  
  
  all_results_ANOVA <- do.call(rbind.data.frame, lapply(results, safe_extract, name = "ANOVA_results"))
  all_results_ANOVA_means <- do.call(rbind.data.frame, lapply(results, safe_extract, name = "ANOVA_means"))
  
  all_results_marginal <- do.call(rbind.data.frame, lapply(results, safe_extract,name= "marginal_results"))
  all_results_collapsed <- do.call(rbind.data.frame, lapply(results, safe_extract, name="collapsed_results"))
  all_results_MMRM_means <- do.call(rbind.data.frame, lapply(results, safe_extract, name="MMRM_means"))
  
  
  return(list(all_results_ANOVA, all_results_ANOVA_means, all_results_marginal, all_results_collapsed, all_results_MMRM_means))
  
}
