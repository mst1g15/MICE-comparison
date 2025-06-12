#-------------------------------------------------------------------------------
#This R script contains a function run the simulation. It is set up for 
#parallelisation within the HPC 
#-------------------------------------------------------------------------------

run_sim <- function(td=NULL, ss, miss, nsim, method, tuning, response, interact=NULL){
  
  
  #each repetition of simulation 
  results <-  foreach(i=1:nsim) %dopar% {
    
    source("00_init.R")
    source("01_generate_complete_data.R")
    source("01_simulate_missing_data.R")
    source("02_impute_functions.R")
    source("03_simulation_function.R")
    
    
    #generate full dataset ---------------------------------------------------------------------
    if(interact==TRUE){
      #single covariate with treatment-covariate interaction 
      df <- gen_cts_inter(ss,  true_mean=FALSE, po=FALSE)
    }else{
      # single covariate, no interaction 
      df <- gen_cts_normal(n=ss, sd=42, true_diff=td, true_mean=FALSE, po=FALSE)
    }
    
    
    #generate missingness----------------------------------------------------------------------- 
    df_miss <- gen_miss(df, response, alpha0=as.numeric(miss_mech_list[miss][[1]][1]), 
                        alpha1=as.numeric(miss_mech_list[miss][[1]][2]), 
                        alpha2=as.numeric(miss_mech_list[miss][[1]][3]))
    
    
    #split by arm for imputation ----------------------------------------------------------
    df_miss_0 <- df_miss %>% filter(treat==0)
    df_miss_1 <- df_miss %>% filter(treat==1)
    
    #use try in case of imputation failing (e.g. with rf)
    final_results <- try(impute_by_arm_method(df_miss_0, df_miss_1, response, method, tuning, interact))
    
    #if this did not run (for example with rf), outcomes are all NA 
    if (!is.data.frame(final_results)){
      final_results <- data.frame(term="treat", estimate=NA, std.error=NA, p.value=NA, "2.5 %"= NA, "97.5 %"= NA, 
                                  response=response, method=method, tuning=tuning)
      colnames(final_results)[5] = "2.5 %"
      colnames(final_results)[6] = "97.5 %"
      
    }
    
    final_results$sim=i
    final_results
    
  }
  
  
  
  all_results <- do.call(rbind.data.frame, results)
  
  #include settings in table of results 
  all_results$miss = miss 
  all_results$ss = ss
  all_results$true_diff = td
  all_results$response=response
  
  
  
  return(all_results)
  
}



