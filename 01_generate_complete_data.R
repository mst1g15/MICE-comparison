#-------------------------------------------------------------------------------
#This R script contains functions for the data generating mechanisms 
#-------------------------------------------------------------------------------

#single covariate setting 
gen_cts_normal<-function(n, sd=42, true_diff=40, true_mean=FALSE, po=FALSE){
  
  # Randomised treatment 
  treat <- rep_len(0:1,n)
  
  # Covariate
  x  <- rnorm(n, 0, 1)
  
  # Outcomes
  
  # Generic error term 
  e1 <- rnorm(n, 0, sd)
  
  # Outcome 1: Linear
  y1_0_mean <- 395 + 110*x
  y1_1_mean <- 395 + 110*x + true_diff 
  
  y1_0 <- y1_0_mean + e1 
  y1_1 <- y1_1_mean + e1 
  
  y1 <- y1_0*(1-treat) + y1_1*treat
  
  
  # Outcome 2: Flattening-off  
  y2_0_mean <- 700 - exp(-(x-4))   
  y2_1_mean <- 700 - exp(-(x-4)) + true_diff
  
  y2_0 <- y2_0_mean + e1 
  y2_1 <- y2_1_mean + e1 
  
  y2 <- y2_0*(1-treat) + y2_1*treat
  
  # Outcome 3: Two-tier  
  y3_1_mean <- 180 + 470*(x>0) + true_diff
  y3_0_mean <- 180 + 470*(x>0)
  
  y3_0 <- y3_0_mean + e1 
  y3_1 <- y3_1_mean + e1 
  
  y3 <- y3_0*(1-treat) + y3_1*treat 
  
  
  # Outcome 4: Harmonic
  y4_1_mean <- 400 + 300*cos(2*pi*3*x*0.15 + 4) + true_diff
  y4_0_mean <- 400 + 300*cos(2*pi*3*x*0.15 + 4) 
  
  y4_0 <- y4_0_mean + e1 
  y4_1 <- y4_1_mean + e1 
  
  y4 <- y4_0*(1-treat) + y4_1*treat
  
  
  # Outcome 4: J-shaped 
  y5_0_mean <- 200 + 100*x + 60*x*x 
  y5_1_mean <- 200 + 100*x + 60*x*x + true_diff
  
  
  y5_0 <- y5_0_mean + e1 
  y5_1 <- y5_1_mean + e1 
  
  y5 <- y5_0*(1-treat) + y5_1*treat 
  
  # Outcome 5: Full (i.e. symmetric) quadratic 
  y6_1_mean <- 100 + 104*x*x    + true_diff
  y6_0_mean <- 100 + 104*x*x 
  
  y6_0 <- y6_0_mean + e1 
  y6_1 <- y6_1_mean + e1 
  
  y6 <- y6_0*(1-treat) + y6_1*treat 
  
  # Keep required variables
  d <- as.data.frame(cbind(treat, x, y1, y2, y3, y4, y5, y6))
  
  if (po==TRUE) {
    d<- as.data.frame(cbind(d, y1_0, y2_0, y3_0, y4_0, y5_0, y6_0,  
                            y1_1, y2_1, y3_1, y4_1, y5_1, y6_1))
  }
  if (true_mean==TRUE) {
    d<- as.data.frame(cbind(d, y1_0_mean, y2_0_mean, y3_0_mean, y4_0_mean, y5_0_mean, y6_0_mean, 
                            y1_1_mean, y2_1_mean, y3_1_mean, y4_1_mean, y5_1_mean, y6_1_mean))
  }
  
  # Return the dataset needed for analysis            
  return(d)
}




#interaction setting
gen_cts_inter <- function(n, true_mean=FALSE, po=FALSE){
  
  # Randomised treatment 
  treat <- rep_len(0:1,n)
  
  # covariate
  x  <- rnorm(n, 0, 1)
  c  <- x**2 
  
  # Outcome variables
  
  # Generic error term 
  e1 <- rnorm(n, 0, 60)
  
  # Outcome 1: Linear in both arms, quantitative interaction
  y1_0_mean <- 395 + 80*x  
  y1_1_mean <- 395 + 80*x + 20*x + 40
  
  y1_0 <- y1_0_mean + e1
  y1_1 <- y1_1_mean + e1
  
  y1 <- y1_0*(1-treat) + y1_1*treat
  
  # Outcome 2: Linear in both arms, qualitative interaction (i.e. switch direction)
  y2_0_mean <- 500 + 20*x 
  y2_1_mean <- 500 + 20*x + 90*x - 65
  
  y2_0 <- y2_0_mean + e1
  y2_1 <- y2_1_mean + e1
  
  y2 <- y2_0*(1-treat) + y2_1*treat
  
  # Outcome 3: Linear in one arm, exponential in other 
  y3_0_mean <- 450 + 100*x 
  y3_1_mean <- 450 + 100 - exp(-(x-3.5)) 
  
  y3_0 <- y3_0_mean + e1
  y3_1 <- y3_1_mean + e1
  
  y3 <- y3_0*(1-treat) + y3_1*treat
  
  
  # Outcome 4: Exponential in one arm, absent in other
  y4_0_mean <- 0
  y4_1_mean <- 13*(c**2) 
  
  y4_0 <- y4_0_mean + e1
  y4_1 <- y4_1_mean + e1
  
  y4 <- y4_0*(1-treat) + y4_1*treat
  
  
  ### Keep required variables
  d <- as.data.frame(cbind(treat, x, y1, y2, y3, y4))
  
  if (po==TRUE) {
    d<- as.data.frame(cbind(d, y1_0, y2_0, y3_0, y4_0, 
                            y1_1, y2_1, y3_1, y4_1))
  }
  if (true_mean==TRUE) {
    d<- as.data.frame(cbind(d, y1_0_mean, y2_0_mean, y3_0_mean, y4_0_mean, 
                            y1_1_mean, y2_1_mean, y3_1_mean, y4_1_mean))
  }
  
  
  # Return the dataset needed for analysis            
  return(d)
}


