
miss_mech_list <- list(
  mcar10 = c(-2.197, 0, 0),
  mcar20= c(-1.386 , 0, 0),
  mcar30 = c(-0.847 , 0, 0), 
  mar10 = c(-2.5, -1, 0),
  mar30 = c(-1, -1, 0),
  mnar10 = c(-6.25, 0, 5),
  mnar30= c(-6.25, 0, 6.95)
)

df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
miss="mar_x_30"
df_miss <- gen_miss(df, "y1", alpha0=as.numeric(miss_mech_list[miss][[1]][1]), 
                    alpha1=as.numeric(miss_mech_list[miss][[1]][2]), 
                    alpha2=as.numeric(miss_mech_list[miss][[1]][3]))

sum(is.na(df_miss$y1) &(df_miss$x>0))/1000000
sum(is.na(df_miss$y1) &(df_miss$x<0))/1000000


sum(is.na(df_miss$y1) &(df_miss$t==0))/1000000
sum(is.na(df_miss$y1) &(df_miss$t==1))/1000000


df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", 
                    alpha0=-1, 
                    alpha1= -1, 
                    alpha2=0)
sum(is.na(df_miss$y1) &(df_miss$x>0))/1000000
sum(is.na(df_miss$y1) &(df_miss$x<0))/1000000

df_miss <- gen_miss(df, "y5", 
                    alpha0=-1, 
                    alpha1= -1, 
                    alpha2=0)
sum(is.na(df_miss$y5) &(df_miss$x>0))/1000000
sum(is.na(df_miss$y5) &(df_miss$x<0))/1000000



df_miss <- gen_miss(df, "y1", 
                    alpha0=-1.55, 
                    alpha1= 0, 
                    alpha2=-1.9)
sum(is.na(df_miss$y1))/1000000

sum(is.na(df_miss$y1) &(df_miss$treat==0))/1000000
sum(is.na(df_miss$y1) &(df_miss$treat==1))/1000000

df_miss <- gen_miss(df, "y1", 
                    alpha0=-0.4, 
                    alpha1= 0, 
                    alpha2=-1)
sum(is.na(df_miss$y1))/1000000

sum(is.na(df_miss$y1) &(df_miss$x>0))/1000000
sum(is.na(df_miss$y1) &(df_miss$x<0))/1000000

sum(is.na(df_miss$y1) &(df_miss$treat==0))/1000000
sum(is.na(df_miss$y1) &(df_miss$treat==1))/1000000

df_miss <- gen_miss(df, "y5", 
                    alpha0=-0.4, 
                    alpha1= 0, 
                    alpha2=-1)
sum(is.na(df_miss$y5) &(df_miss$x>0))/1000000
sum(is.na(df_miss$y5) &(df_miss$x<0))/1000000
sum(is.na(df_miss$y5) &(df_miss$treat==0))/1000000
sum(is.na(df_miss$y5) &(df_miss$treat==1))/1000000







df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", 
                    alpha0=-6.25, 
                    alpha1=  0, 
                    alpha2=5)
sum(is.na(df_miss$y1))/1000000


df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", 
                    alpha0=-1, 
                    alpha1=  -1, 
                    alpha2=0)
sum(is.na(df_miss$y1))/1000000
sum(is.na(df_miss$y1 & df_miss$x>0))/1000000
sum(is.na(df_miss$y1 & df_miss$x<0))/1000000





df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", 
                    alpha0=-1.55, 
                    alpha1=  0, 
                    alpha2=-1.9)
sum(is.na(df_miss$y1))/1000000
sum(is.na(df_miss$y1 & df_miss$treat==1))/1000000
sum(is.na(df_miss$y1 & df_miss$treat==0))/1000000




df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", 
                    alpha0=-0.4, 
                    alpha1=  0, 
                    alpha2=-1)
sum(is.na(df_miss$y1))/1000000
sum(is.na(df_miss$y1 & df_miss$treat==1))/1000000
sum(is.na(df_miss$y1 & df_miss$treat==0))/1000000



hist(df$y1)
hist(df$y1[is.na(df_miss$y1)], add=T, col="blue")




df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", 
                    alpha0=-6.25, 
                    alpha1=  0, 
                    alpha2=6.95)
sum(is.na(df_miss$y1))/1000000


hist(df$y1)
hist(df$y1[is.na(df_miss$y1)], add=T, col="blue")




sum(is.na(df_miss$y1)/1000000)

df$y1

#test proportion of missing under each mechanism. 
df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", alpha0=-2.5, 
                    alpha1= -1, 
                    alpha2=0)
sum(is.na(df_miss$y1))/1000000

df <- gen_cts_normal(1000000, sd=42, true_diff=40, true_mean=FALSE, po=FALSE)
df_miss <- gen_miss(df, "y1", alpha0=-0.45, 
                    alpha1= -0.3, 
                    alpha2= -0.001)
sum(is.na(df_miss$y1))/1000000

1/(1+exp(-(-2.2-0.3-300*0.001)))
1/(1+exp(-(-2.2-0.3)))

1/(1+exp(2.944+0.5))

1/(1+exp(1.386+1))
1/(1+exp(1.386))


pmm_settings <- expand.grid(method="pmm", y_var_all)



#create all_settings----------------------------------------------------------------------------------
#single covariate setting, no interaction
response <- c("y1", "y2", "y3", "y4", "y5", "y6")
miss_mech <- c("mcar10", "mcar30", "mar_x_10", "mar_x_30",  "mar_t_10", "mar_t_30")

sample_size <- c(50, 100, 200, 500)
true_diff <- c(0, 40)

settings_include_mice <- data.frame(method=c("norm_nonlinear", "pmm_nonlinear", "norm_default", "pmm_default", "cc"), tuning=NA)
settings_include_rf <- expand.grid(method=c("rf"), tuning = c(5, 10, 20))
settings_include_cart <- expand.grid(method=c("cart"), tuning = c(5, 10, 20))
settings_include_superlearner <- data.frame(method="superlearner", tuning=NA)

settings_temp <- rbind(settings_include_mice, settings_include_rf, settings_include_cart, settings_include_superlearner)

settings_temp1 <- expand.grid(sample_size=sample_size, miss_mech=miss_mech, true_diff=true_diff, response=response)

all_settings <- merge(settings_temp, settings_temp1)
all_settings <- all_settings %>% filter(!(response=="y1" & method=="norm_nonlinear")) %>%    #for y1, no need to look at non-linear 
                                 filter(!(response=="y1" & method=="pmm_nonlinear"))

all_settings$interact=FALSE
all_settings$response <- as.character(all_settings$response)

#interaction settings 

all_settings_interact <- all_settings
all_settings_interact$interact <- TRUE
all_settings_interact <- all_settings_interact %>% filter(response %in% c("y1", "y2", "y3", "y4")) %>%
  filter(true_diff==40)
all_settings_interact <- all_settings_interact %>%  filter(!(response=="y2" & method=="norm_nonlinear")) %>%    #y2 is nonlinear for the interaction
                                  filter(!(response=="y2" & method=="pmm_nonlinear"))

                                                          

all_settings_combined <- rbind(all_settings, all_settings_interact)                                                           


saveRDS(all_settings_combined, "settings/all_settings.RDS")
