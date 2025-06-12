#-------------------------------------------------------------------------------
#This R script contains a function to generate missingness in the dataset 
#-------------------------------------------------------------------------------

gen_miss <- function(df, y, alpha0, alpha1, alpha2){
  

  des_vec <- cbind(1, df$x, df$treat) #one row of design matrix 
  
  f <-   as.matrix(des_vec) %*% (matrix(c(alpha0, alpha1, alpha2))) #linear predictor
  probs <-  1/(1+exp(-f))
  indic <- rbinom(length(f), 1, 1/(1+exp(-f)))

  
  df[[y]] <-  ifelse(indic, NA, df[[y]])
  
  return(df) 
  
}
