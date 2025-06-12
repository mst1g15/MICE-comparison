#-------------------------------------------------------------------------------
#This R script loads libraries and parameter settings for the ML MICE simulation study
#-------------------------------------------------------------------------------

#load libraries ------------------------------
library(mice)
library(tidyverse)
library(VIM)
library(rlang)
library(purrr)
library(gridExtra)
library(doParallel)
library(foreach)
library(superMICE)
library(randomForest)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(nlme)
library(emmeans)
library(mitools)


#parameter settings ------------------------------
#values for missing data mechanism: intercept, coefficient of covariate
# and coefficient of response (in case of MNAR mechanism)
miss_mech_list <- list(
  mcar10 = c(-2.197, 0, 0),
  mcar20= c(-1.386 , 0, 0),
  mcar30 = c(-0.847 , 0, 0), 
  mar_x_10 = c(-2.5, -1, 0),
  mar_x_30 = c(-1, -1, 0),
  mar_t_10= c(-1.55, 0, -1.9),
  mar_t_30 = c(-0.4, 0, -1)
)


#true values for interaction terms---------------
true_inter <- list(y1=40.00227, 
                   y2=-64.98979, 
                   y3=45.39995, 
                   y4=39.03279)


