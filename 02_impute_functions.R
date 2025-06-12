#-------------------------------------------------------------------------------
#This R script contains a function to which performs MI separately by arm 
# methods include: "pmm", "norm, "rf", "cart" and "superlearner"
#-------------------------------------------------------------------------------


impute_by_arm_method <- function(df_miss_0, df_miss_1, y_var, method, tuning=NULL, interact=NULL){
  #input: df_miss_0: dataframe where arm is 0 
  #       df_miss_1: dataframe where arm is 1 
  #       y_var: character for the particular yvar to use (between 1 and 4)
  #       method: character for imputation method: "pmm", "norm" "rf", "cart", "superlearner
  #       tuning: optional input for cart and rf - ML tuning parameters 
  
  
  #assume that the final analysis model will regress on 
  #covariate and treatment only (not the true relationship)---------------------
  
  #complete cases 
  if(method=="cc"){
    
    totaldat  <- rbind(df_miss_0, df_miss_1)
    totaldat <- totaldat[,c("treat", "x", y_var)]
    
    cc_analysis <- lm(as.formula(paste(y_var,  "~ x + treat")), data=totaldat)
    
    
    summary1 <- summary(cc_analysis)$coef[,c("Estimate", "Std. Error", "Pr(>|t|)")]
    summary2 <- confint(cc_analysis)
    summary_results <- data.frame(rownames(summary1), summary1, summary2)
    
    colnames(summary_results) <- c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")
    
    summary_results$response <- y_var
    summary_results$method <- "cc"
    summary_results$tuning <- NA
  }
  
  
  
  if(method %in% c("superlearner", "pmm_default", "norm_default", "rf", "cart")){
    #assume linear model in imputation model and outcome model 
    
    #imputation by superlearner 
    if(method=="superlearner"){
      
      SL.lib <- c("SL.mean", "SL.glm", "SL.gam", "SL.randomForest")
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("x",  y_var)], m=30, method="SuperLearner", SL.lib=SL.lib, kernel="gaussian", bw=c(0.25, 1, 5), print=FALSE)
      imp_1 <- mice(df_miss_1[, c("x",  y_var)], m=30, method="SuperLearner", SL.lib=SL.lib, kernel="gaussian", bw=c(0.25, 1, 5), print=FALSE)
      
    }else if(method=="pmm_default"){
      
      imp_0 <- mice(df_miss_0[, c("x",  y_var)], m = 30, method="pmm", print=F)
      imp_1 <- mice(df_miss_1[, c("x",  y_var)], m = 30, method="pmm", print=F)
      
    }else if(method=="norm_default"){
      
      imp_0 <- mice(df_miss_0[, c("x",  y_var)], m = 30, method="norm", print=F)
      imp_1 <- mice(df_miss_1[, c("x",  y_var)], m = 30, method="norm", print=F)
      
    }else if(method == "rf"){
      
      imp_0 <- mice(df_miss_0[, c("x",  y_var)], m = 30, method="rf", print=F, ntree=tuning)
      imp_1 <- mice(df_miss_1[, c("x",  y_var)], m = 30, method="rf", print=F, ntree=tuning)
      
    }else if(method=="cart"){
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("x",  y_var)], m = 30, method="cart", minbucket=tuning, print=F)
      imp_1 <- mice(df_miss_1[, c("x",  y_var)], m = 30, method="cart", minbucket=tuning, print=F)
      
    }
    
    #impute in teach arm, then add treatment variable and merge datasets from two arms 
    complete_0 <- complete(imp_0, include=TRUE, action="long")
    complete_0$treat = 0
    complete_0 <- as.mids(complete_0)
    
    complete_1 <- complete(imp_1, include=TRUE, action="long")
    complete_1$treat = 1
    complete_1 <- as.mids(complete_1)
    
    imputed_dat <- rbind(complete_0, complete_1)
    
    
    # Fit the linear model for the current y variable - adjusting for the covariate and treatment 
    fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ x + treat")))) 
    
    # Pool the results of the model
    pooled_results <- pool(fit)
    
    # Summarize the pooled results and add a column for the response variable
    summary_results <- summary(pooled_results, conf.int=TRUE)[,c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")]
    summary_results$response <- y_var
    summary_results$method <- method
    summary_results$tuning <- tuning
    
  }
  
  
  if(method=="pmm_nonlinear" | method=="norm_nonlinear"){
    #now assume the true response in the imputation and outcome model 
    
    impute_method = ifelse(method=="pmm_nonlinear", "pmm", "norm")
    
    if(interact=="FALSE" & y_var=="y2"){
      
      df_miss_0$xsq <-exp(-(df_miss_0$x))
      df_miss_1$xsq <-exp(-(df_miss_1$x))
      
      imp_0 <- mice(df_miss_0[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ xsq+treat")))) 
      
    }else if(interact==FALSE &y_var=="y3"){
      
      df_miss_0$xsq <-df_miss_0$x>0
      df_miss_1$xsq <-df_miss_1$x>0
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ xsq+treat")))) 
      
    }else if(interact==FALSE & y_var=="y4"){
      
      df_miss_0$xsq <-  cos(2*pi*3*df_miss_0$x*0.15 + 4)
      df_miss_1$xsq <-  cos(2*pi*3*df_miss_1$x*0.15 + 4)
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ xsq+treat")))) 
      
      
    }else if(interact == FALSE & y_var=="y5"){
      
      df_miss_0$xsq <-  df_miss_0$x^2
      df_miss_1$xsq <-  df_miss_1$x^2
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("x", "xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("x", "xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ x + xsq+treat")))) 
      
    }else if (interact == FALSE & y_var=="y6"){
      
      df_miss_0$xsq <-  df_miss_0$x^2
      df_miss_1$xsq <-  df_miss_1$x^2
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ xsq+treat")))) 
      
    }else if(interact==TRUE & y_var=="y3"){
      
      
      df_miss_0$xsq <-df_miss_0$x^2     #######not sure about this! 
      df_miss_1$xsq <-df_miss_1$x^2
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("x", "xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("x", "xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ x + treat"))))
      
    }else if(interact==TRUE & y_var=="y4"){
      
      
      df_miss_0$xsq <-df_miss_0$x^4     #######not sure about this! 
      df_miss_1$xsq <-df_miss_1$x^4
      
      # Perform the multiple imputation for the current y variable
      imp_0 <- mice(df_miss_0[, c("x", "xsq", y_var)], m = 30, method=impute_method, print=F)
      imp_1 <- mice(df_miss_1[, c("x", "xsq", y_var)], m = 30, method=impute_method, print=F)
      
      complete_0 <- complete(imp_0, include=TRUE, action="long")
      complete_0$treat = 0
      complete_0 <- as.mids(complete_0)
      
      complete_1 <- complete(imp_1, include=TRUE, action="long")
      complete_1$treat = 1
      complete_1 <- as.mids(complete_1)
      
      imputed_dat <- rbind(complete_0, complete_1)
      
      # Fit the linear model for the current y variable
      fit <- with(data= imputed_dat, expr= lm(as.formula(paste(y_var,  "~ x + treat")))) 
      
      
    }
    
    # Pool the results of the model
    pooled_results <- pool(fit)
    
    # Summarize the pooled results and add a column for the response variable
    summary_results <- summary(pooled_results, conf.int = TRUE)[,c("term", "estimate", "std.error", "p.value", "2.5 %", "97.5 %")]
    summary_results$response <- y_var
    summary_results$method <- method
    summary_results$tuning=NA
    
  }  
  
  
    return(summary_results)
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #
  # 
  # #k nearest neighbors-----------------------------------------------------------
  # #note that this is a method that uses single imputation only.
  # #k, the number of neighbours, included as an option
  # if(method=="knn"){
  #   
  #   imp_knn_0 <- kNN(df_miss_0, k=tuning)
  #   imp_knn_0$treat=0
  #   imp_knn_1 <- kNN(df_miss_1, k=tuning)
  #   imp_knn_1$treat <- 1
  #   imp_knn_dat <- rbind(imp_knn_0, imp_knn_1)
  #   
  #   knn_result <- lm(as.formula(paste(y_var, "~ x+treat")), data=imp_knn_dat)
  #   
  #   summary_results <- summary(knn_result)$coefficients[,-3]
  #   summary_results <- rownames_to_column(as.data.frame(summary_results), var="term")
  #   colnames(summary_results) <- c("term", "estimate", "std.error", "p.value")
  #   intervals <- confint(knn_result)
  #   summary_results <- cbind(summary_results, intervals)
  #   summary_results$response <- y_var
  #   summary_results$method <- "knn"
  #   summary_results$tuning <- tuning
  #   
  # }
  # 
  # 
  