#-----------------------------------------------------------------------------
#this R script contains functions required to run the simulation based on the 
#inopulse trial 
#----------------------------------------------------------------------------

calc_lognorm_param <- function(mu, sigma){
  #Purpose: given a log-normal distribution X with mean mu and standard deviation sigma,
  # provide the location (mu_x) and scale (sigma_x) parameters of X.
  
  mu_x=log((mu^2)/sqrt(mu^2+sigma^2))
  sigma_x=sqrt(log(1+(sigma^2)/(mu^2)))
  
  return(list(mu_x, sigma_x))
}

generate_daily_data <- function(ss, weekno=4, rand_ratio=NULL){
  #Purpose: Generate daily time spent in MVPA for a setting inspired by the Bellerophon Phase II trial on INOPulse 
  #Inputs:  season_prop: fraction of participants who experience a seasonal effect 
  #         weekno: number of weeks for follow-up
  #         rand_ratio: randomization ratio, assumed by default as 2:1 for treatment:control and 1:1 otherwise. 
  #         season_effect: scalar for the effect of season 
  #         observer effect: scalar for the effect of observer 
  #         miss_type: set either to "MCAR" (Missing Completely at Random) or "MNAR" (Missing Not at Random)
  #         miss_prop_n: fraction of participants who have missing data
  #         miss_prop_days: fraction of days in the measurement period which are non-compliant 
  #Outputs: dataset for one run of simulation with daily MVPA 
  
  
  #generate ID: 44 patients. Assume 28 days (four weeks) of data per patient 
  dayno=7*weekno
  ID <- rep(1:ss, each=dayno)
  
  #tratment assignment vector 
  if(is.null(rand_ratio)){
    #assume by default a 2:1 ratio to treatment:control
    treat_allocation= sort((sample(c(1, 0), ss, replace=TRUE, prob=c(2/3, 1/3))), decreasing=TRUE)
  }else{
    #use a 1:1 ratio for treatment:control
    treat_allocation=sort((sample(c(1, 0), ss, replace=TRUE, prob=c(1/2, 1/2))), decreasing=TRUE)
  }
  
  # based on the distribution of MVPA at baseline in Table 1 of paper 
  baseline <- rep(rlnorm(ss, 4.155, 0.6129), each=dayno)
  
  treat <- rep(treat_allocation, each=dayno)
  week <- as.factor(rep(rep(1:weekno, each=7), ss))
  
  #treat_effect <- rep(12.5*treat, each=dayno)
  #treat_time <-5* treat_effect*(week%in% c(3, 4))
  
  mu <- baseline + 12*treat+#treat_time 
    rep(rnorm(ss,0, 2), each=dayno) 
  
  
  #vanilla version of the outcome: 
  
  params <- calc_lognorm_param(mu, 46)
  
  daily_outcome <- rlnorm(n=dayno*ss, params[[1]], params[[2]])
  
  #generate complete data
  dat<- tibble( ID=ID, daily_outcome=daily_outcome, 
                treat=treat, baseline=baseline,
                week=week, day=rep(1:dayno, ss)) %>% group_by(ID)
  
  
  dat_wide <- dat[, !names(dat) %in% "week"]  %>%
    pivot_wider(names_from = day, values_from = daily_outcome,  names_prefix = "Day_")
  return(list(dat, dat_wide))
}



induce_missingness <- function(dat, miss_type, miss_prop_n=0.3, miss_week_prop=c(0.25, 0.5, 0.75, 0.85)){
  #dat in wide format 
  ss <- nrow(dat)
  
  num_ind <- round(miss_prop_n*ss)
  
  
  if(miss_type=="MCAR"){
    
    
    select_id <- sample(1:ss, num_ind)
    
    for(p in select_id){
      
      miss_weeks <- rbinom(4, 1, miss_week_prop)
      
      select_dates <- c(miss_weeks[1] * (1:7), 
                        miss_weeks[2] * (8:14),
                        miss_weeks[3] * (15:21),  
                        miss_weeks[4] * (22:28))
      
      
      cols <- which(names(dat) %in% paste0("Day_", select_dates[select_dates!=0]))
      
      dat[p, cols] <- NA
      
    }
    
  }
  
  
  if(miss_type=="MAR_baseline"){
    
    mean_baseline <- mean(dat$baseline)
    
    
    base0 <- which(dat$baseline<mean_baseline)
    base1 <- which(dat$baseline>mean_baseline)
    
    
    select0 <- round((2/3)*num_ind)
    select1 <- num_ind - select0
    
    select_id <- c(sample(base0, select0), sample(base1, select1))
    
    
    for(p in select_id){
      
      miss_weeks <- rbinom(4, 1, miss_week_prop)
      
      select_dates <- c(miss_weeks[1] * (1:7), 
                        miss_weeks[2] * (8:14),
                        miss_weeks[3] * (15:21),  
                        miss_weeks[4] * (22:28))
      
      
      cols <- which(names(dat) %in% paste0("Day_", select_dates[select_dates!=0]))
      
      dat[p, cols] <- NA
      
    }
    
    
    
  }
  
  if(miss_type=="MAR_treat"){
    
    arm0 <- which(dat$treat==0)
    arm1 <- which(dat$treat==1)
    
    select0 <- round((2/3)*num_ind)
    select1 <- num_ind - select0
    
    select_id <- c(sample(arm0, select0), sample(arm1, select1))
    
    
    for(p in select_id){
      
      miss_weeks <- rbinom(4, 1, miss_week_prop)
      
      select_dates <- c(miss_weeks[1] * (1:7), 
                        miss_weeks[2] * (8:14),
                        miss_weeks[3] * (15:21),  
                        miss_weeks[4] * (22:28))
      
      
      cols <- which(names(dat) %in% paste0("Day_", select_dates[select_dates!=0]))
      
      dat[p, cols] <- NA
      
    }
    
    
  }
  
  
  return(dat)
  
}

safe_lme <- function(fixed, data, random, ..., correlation = NULL, weights = NULL, control = lmeControl()) {
  tryCatch({
    lme(
      fixed = fixed,
      data = data,
      random = random,
      correlation = correlation,
      weights = weights,
      control = control,
      ...
    )
  }, error = function(e) {
    message("Error during model fitting: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("Warning during model fitting: ", w$message)
    invokeRestart("muffleWarning")
  })
}
impute_and_analyse <- function(dat, method, dat_true){
  #Purpose: run repetitions under a simulation scenario and save results 
  #Inputs:  reps: number of repetitions 
  #         scenario_name: character for scenario name 
  #         additional inputs: include inputs to generate_daily_data
  #Outputs: mean and standard error of the treatment effect of the primary analysis,  
  #         and their Monte Carlo Errors 
  
  #IMPUTE STEP 
  #wide format data 
  
  
  ss <- nrow(dat)
  
  if(method=="cc"){
    
    #convert to long format 
    summary_dat <- dat %>% pivot_longer(cols=starts_with("Day_"), names_to="Day", values_to="Daily_outcome")
    
    #ANOVA analysis, complete case ---------------------------------------------------------------
    month_summary <- summary_dat %>% group_by(ID, baseline, treat) %>% summarize(month_average=mean(Daily_outcome, na.rm=TRUE))
    
    model <- lm(month_average~baseline+treat, data=month_summary)
    
    summary1 <- summary(model)$coef[,c("Estimate", "Std. Error", "Pr(>|t|)")]
    summary2 <- confint(model)
    summary_ANOVA <- data.frame(rownames(summary1), summary1, summary2)
    colnames(summary_ANOVA) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
    summary_ANOVA$method="cc"
    
    
    ANOVA_means_table <- data.frame(emmeans(model, ~ treat))
    colnames(ANOVA_means_table) <- c("term", "estimate", "std.error", "df", "2.5 %", "97.5 %")
    ANOVA_means_table$method <- "default"
    
    
    #MMRM, default: missing weeks do not contribute to likelihood-------------------------------------
    summary_dat$week <- rep(rep(1:4, each=7), ss)
    summary_dat$week <- as.factor(summary_dat$week)
    week_summary <- summary_dat %>% group_by(week, ID, treat, baseline) %>% summarise(week_average = mean(Daily_outcome))
    
    #discard single weeks (not entire subjects) with missing data: 
    week_summary<- week_summary[complete.cases(week_summary),]
    
    week_summary$treat <- as.factor(week_summary$treat)
    week_summary$treat <- relevel(week_summary$treat, ref="1")
    
    #run MMRM model 
    mmrm_model <- safe_lme(
      week_average ~ baseline + treat * week,
      data = week_summary,
      random = ~ 1 | ID,
      weights = varIdent(form = ~ 1 | week),
      correlation = corSymm(form = ~ 1 | ID),
      control = lmeControl(maxIter = 100, msMaxIter = 100),
      method = "REML"
    )
    
    #create dummy matrices if model does not run - to help with analyses
    if(is.null(mmrm_model)){
      
      marginal_means <-  data.frame(term=1:4, matrix(NA, nrow = 4, ncol = 5))
      colnames(marginal_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      marginal_means$method <- "default"
      
      collapsed_means <-  data.frame(term="treat", matrix(NA, nrow = 1, ncol = 5))
      colnames(collapsed_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      collapsed_means$method <- "default"
      
      MMRM_means_table <-  data.frame(treat=rep(c(0, 1), 4), week=rep(1:4, each=2), matrix(NA, nrow = 8, ncol = 5))
      colnames(MMRM_means_table) <- c("term", "week", "estimate", "std.error", "df", "2.5 %", "97.5 %")
      MMRM_means_table$method <- "default"
      
    }else{
      
      #estimated marginal means 
      emm <- contrast(emmeans(mmrm_model, ~ treat | week, data=week_summary),  method = "pairwise")
      marginal_means <- data.frame(emm)[,c("week", "estimate", "SE", "p.value")]
      marginal_means_ci <- data.frame(confint(emm))[,c("lower.CL", "upper.CL")]
      marginal_means <- cbind(marginal_means, marginal_means_ci)
      colnames(marginal_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      marginal_means$method <- "default"
      
      #marginal means collapsed over weeks
      emm_collapsed <- emmeans(mmrm_model, ~ treat, data=week_summary)
      collapsed_means <- contrast(emm_collapsed, method = "pairwise")
      collapsed_means_df <- data.frame(collapsed_means)[, c("contrast", "estimate", "SE", "p.value")]
      collapsed_means_ci <- data.frame(confint(collapsed_means))[,c("lower.CL", "upper.CL")]
      collapsed_means <- cbind(collapsed_means_df, collapsed_means_ci)
      colnames(collapsed_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      collapsed_means$term="treat"
      collapsed_means$method <- "default"
      
      #estimated means for two arms per week
      MMRM_means_table <- data.frame(emmeans(mmrm_model, ~ treat*week, data=week_summary))
      colnames(MMRM_means_table) <- c("term", "week", "estimate", "std.error", "df", "2.5 %", "97.5 %")
      MMRM_means_table$method <- "default"
      
    }
    
  }
  
  
  #Imputation methos-------------------------------------------------------------
  if(method %in% c("rf", "cart", "pmm", "norm")){
    
    
    #impute at week level 
    dat_wk <-dat %>%  mutate(week1 = rowMeans(pick(Day_1:Day_7)), 
                             week2 = rowMeans(pick(Day_8:Day_14)), 
                             week3 = rowMeans(pick(Day_15:Day_21)),
                             week4 = rowMeans(pick(Day_22:Day_28)))
    dat_wk <- dat_wk %>% dplyr::select(-starts_with("Day_"))

    
    #imputation step----------------------------------------------------------
    dat_miss_0 <- dat_wk %>% filter(treat==0)
    dat_miss_1 <- dat_wk %>% filter(treat==1)
    
    if(method %in% c("pmm", "norm")){
      
      imp_0 <- mice(dat_miss_0, m = 30, method=method, print=F)
      imp_1 <- mice(dat_miss_1, m = 30, method=method, print=F)
      
    }else if (method =="rf"){
      
      imp_0 <- mice(dat_miss_0, m = 30, method=method, ntree=5, print=F)
      imp_1 <- mice(dat_miss_1, m = 30, method=method, ntree=5, print=F)
      
    }else if (method=="cart"){
      
      imp_0 <- mice(dat_miss_0, m = 30, method=method, minbucket=5, print=F)
      imp_1 <- mice(dat_miss_1, m = 30, method=method, minbucket=5, print=F)
      
    }
    
    
    
    complete_0 <- complete(imp_0, include=TRUE, action="long")
    complete_1 <- complete(imp_1, include=TRUE, action="long")
    
    
    complete_0_imp <- list()
    
    complete_0_imp<-lapply(0:30, function(i){
      complete_0 %>%
        subset(.imp==i) %>%
        pivot_longer(cols=starts_with("week"), names_to="Week", values_to="Week_outcome") %>%
        mutate(.id = 1:nrow(.))
      
      
    })
    
    imputed_0_long <- as.mids(do.call(rbind, complete_0_imp))
    
    imputed_0_long$imp$Week_outcome
    
    complete_1_imp <- list()
    complete_1_imp<-lapply(0:30, function(i){
      complete_1 %>%
        subset(.imp==i) %>%
        pivot_longer(cols=starts_with("week"), names_to="Week", values_to="Week_outcome") %>%
        mutate(.id = 1:nrow(.))
      
      
    })
    
    imputed_1_long <- as.mids(do.call(rbind, complete_1_imp))
    
    imputed_all <- rbind(imputed_0_long, imputed_1_long)
    imputed_all$data$treat <- as.factor( imputed_all$data$treat)
    imputed_all$data$treat <- relevel( imputed_all$data$treat , ref="1")
    
    fit <- with(data= imputed_all, {
      tryCatch(lme(
                  Week_outcome ~ baseline + treat * Week,
                  random = ~ 1 | ID,
                  correlation = corSymm(form = ~ 1 | ID),
                  weights = varIdent(form = ~ 1 | Week),
                  method = "REML", 
                  control = lmeControl(maxIter = 100, msMaxIter = 100)), 
                 error = function(e) {
                   message("Error during model fitting: ", e$message)
                   return(NULL)
                 },
                 warning = function(w) {
                   message("Warning during model fitting: ", w$message)
                   invokeRestart("muffleWarning")  # suppress warning if you want
                 }
      )
    })
    
    fit_list <- fit$analyses[!sapply(fit$analyses, is.null)]
    
    if(is.null(fit)){
       
      marginal_means <-  data.frame(term=1:4, matrix(NA, nrow = 4, ncol = 5))
      colnames(marginal_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      marginal_means$method <- method
      
      collapsed_means <-  data.frame(term="treat", matrix(NA, nrow = 1, ncol = 5))
      colnames(collapsed_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      collapsed_means$method <- method
      collapsed_means$term="treat"
      
      MMRM_means_table <-  data.frame(treat=rep(c(0, 1), 4), week=rep(1:4, each=2), matrix(NA, nrow = 8, ncol = 5))
      colnames(MMRM_means_table) <- c("term", "week", "estimate", "std.error", "df", "2.5 %", "97.5 %")
      MMRM_means_table$method <- "default"
      
    }else{

      
      #get marginal means -----------------------------------------------
      contrast_list <- lapply(fit_list, function(model) {
        emm <- emmeans(model, ~ treat | Week)
        contrast(emm, method = "pairwise", adjust = "none")
      })
      
      # Step 2: Extract results with p-values and CIs from each contrast
      estimates_list <- lapply(contrast_list, function(res) as.data.frame(res)$estimate)
      se_list <- lapply(contrast_list, function(res) as.data.frame(res)$SE)
      se_list_squared <- lapply(se_list, function(x) x^2)
      df_list <- lapply(contrast_list, function(res) as.data.frame(res)$df)
      
      # Step 3: Pool using MIcombine
      pooled_results <- MIcombine(results = estimates_list, variances = se_list_squared)
      
      estimate=pooled_results$coefficients
      se=sqrt(diag(pooled_results$variance))
      df = pooled_results$df
      ci_lower = estimate - qt(0.975, df) * se
      ci_upper = estimate + qt(0.975, df) * se
      pval = 2 * pt(abs(estimate / se), df, lower.tail = FALSE)
      
      marginal_means <- data.frame(1:4, estimate, se, pval, ci_lower, ci_upper)
      colnames(marginal_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      marginal_means$method <- method
      
      #get marginal means collapsed over weeks---------------------------------------
      
      contrast_list <- lapply(fit_list, function(model) {
        emm <- emmeans(model, ~ treat)
        contrast(emm, method = "pairwise", adjust = "none", reverse=TRUE)
      })
      
      # Step 2: Extract results with p-values and CIs from each contrast
      estimates_list <- lapply(contrast_list, function(res) as.data.frame(res)$estimate)
      se_list <- lapply(contrast_list, function(res) as.data.frame(res)$SE)
      se_list_squared <- lapply(se_list, function(x) x^2)

      
      # Step 3: Pool using MIcombine
      pooled_results <- MIcombine(results = estimates_list, variances = se_list_squared)
      estimate = pooled_results$coefficients
      se=sqrt(diag(pooled_results$variance))
      df = pooled_results$df
      ci_lower = estimate - qt(0.975, df) * se
      ci_upper = estimate + qt(0.975, df) * se
      pval = 2 * pt(abs(estimate / se), df, lower.tail = FALSE)
      
      collapsed_means <- data.frame("mean", estimate, se, pval, ci_lower, ci_upper)
      colnames(collapsed_means) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
      collapsed_means$method <- method
      
      
      #get estimated means per treatment/week -----------------------------------------------
      means_list <- lapply(fit_list, function(model) {
        emm <- emmeans(model, ~ treat | Week)
      })
      
      # Step 2: Extract results with p-values and CIs from each contrast
      estimates_list <- lapply(means_list, function(res) data.frame(res)$emmean)
      se_list <- lapply(means_list, function(res) data.frame(res)$SE)
      se_list_squared <- lapply(se_list, function(x) x^2)
      df_list <- lapply(means_list, function(res) data.frame(res)$df)
      
      # Step 3: Pool using MIcombine
      pooled_results <- MIcombine(results = estimates_list, variances = se_list_squared)
      
      estimate=pooled_results$coefficients
      se=sqrt(diag(pooled_results$variance))
      df = pooled_results$df
      ci_lower = estimate - qt(0.975, df) * se
      ci_upper = estimate + qt(0.975, df) * se
      pval = 2 * pt(abs(estimate / se), df, lower.tail = FALSE)
      
      MMRM_means_table <- data.frame(treat=data.frame(means_list[[1]])[,1], 
                                     week=data.frame(means_list[[1]])[,2],
                                     estimate, se, df, ci_lower, ci_upper)
      colnames(MMRM_means_table) <- c("term", "week", "estimate", "std.error", "df", "2.5 %", "97.5 %")
      MMRM_means_table$method <- method
      
    
    }
    
    ###ANOVA
    
    
    complete_0_month <-complete_0 %>% mutate(month_outcome=rowMeans(pick(week1:week4), na.rm=FALSE)) %>% 
      dplyr::select(-c(week1, week2, week3, week4))
    complete_1_month <- complete_1 %>% mutate(month_outcome=rowMeans(pick(week1:week4), na.rm=FALSE)) %>% 
      dplyr::select(-c(week1, week2, week3, week4)) 
    
    imputed_0_month <- as.mids(complete_0_month)
    imputed_1_month <- as.mids(complete_1_month)
    
    imputed_all_month <- rbind(imputed_0_month, imputed_1_month)
    
    
    fit <- with(imputed_all_month, lm(month_outcome~baseline+treat))
    summary_ANOVA <- summary(pool(fit), conf.int=TRUE)[,c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")]
    summary_ANOVA$method <- method
  
    
    ANOVA_means_table <- data.frame(emmeans(fit, ~ treat))
    colnames(ANOVA_means_table) <- c("term", "estimate", "std.error", "df", "2.5 %", "97.5 %")
    ANOVA_means_table$method <- method
    }
  
  
  return(list(summary_ANOVA, ANOVA_means_table, marginal_means, collapsed_means, MMRM_means_table))
  
}
