# IMPORT LIBRARIES
library(epiDisplay); library(lubridate) ; library(dplyr); library(tidyr); library(plyr)
library(splines); library(gnm); library(Epi); library(tsModel); library(meta); library(sf) 
library(data.table); library(viridis); library(dlnm)
library(mgcv) ; library(corrplot) ; library(mixmeta)
library(ggplot2) ; library(patchwork)


# 100K CA_CAMS AND PRE-PROCESS 
#load regional file
data <- read.csv('CAMS_CA_100k_v1.csv')  
data <- as.data.table(data)


# Define the public holidays in Milan
milan_holidays <- c(
  "2016-01-01", "2016-01-06", "2016-03-28", "2016-04-25", "2016-06-02", "2016-08-15", "2016-11-01", "2016-12-08", "2016-12-26",
  "2017-01-06", "2017-04-17", "2017-04-25", "2017-05-01", "2017-06-02", "2017-08-15", "2017-11-01", "2017-12-08", "2017-12-25", "2017-12-26",
  "2018-01-01", "2018-04-02", "2018-04-25", "2018-05-01", "2018-08-15", "2018-11-01", "2018-12-25", "2018-12-26",
  "2019-01-01", "2019-04-22", "2019-04-25", "2019-05-01", "2019-08-15", "2019-11-01", "2019-12-25", "2019-12-26"
)


# chose data only till 2019 
data$Date <- as.Date(data$Date)
data <- data[data$Date <= as.Date("2019-12-31"), ]
# Edit date format and add additional time related parameters
data$Date <- as.Date(data$Date)
data$year <- as.factor(substr(data$Date,1,4))
data$month <- as.factor(months(data$Date,abbr=TRUE))
data$day <- format(data$Date, format = "%d")
data$dow <- wday(data$Date) 
data$doy <- yday(data$Date)
# Update the "adult" column
data[, adult := adult + young]
# Generate holiday dates for Italy from 2016 to 2019
data$holidays <- as.integer(data$Date %in% as.Date(milan_holidays))
# define stratum
data[, stratum:=factor(paste(DISTRICT, year, month, dow, sep=":"))]


# 2-STAGE MODEL WITH DISTRICT STRATIFICATION AND META-ANALYSIS

# load parameters 
districts <- unique(data$DISTRICT)
lags <- 0:7
argvar <- list(fun = "lin")
arglag <- list(knots=logknots(7, 2)) 

pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO") 
inputs <- c("all_CA", "senior", "adult", "male", "female")  

# Initialize a nested list to store results
results_list <- list()
all_results_lag <- data.frame()
# Process through multiple loops
for (input in inputs) {
  print(paste("Processing input:", input))
  results_list[[input]] <- list()
  for (pollutant in pollutants) {
    print(paste("Processing pollutant:", pollutant))
    coef_vcov_list <- list()
    results_lag <- data.frame()
    for (dist in districts)  {
      data_sub <- data[data$DISTRICT == dist, ]
      data_sub[, keep := sum(get(input)) > 0, by = stratum]
      # accounting for long-term trend
      spldoy <- onebasis(data_sub$doy, "ns", df=3) 
      # adjusting for confounders 
      cb_temp <- crossbasis(data_sub$Temp, lag=14, argvar=list(fun="ns", knots=quantile(data_sub$Temp, c(0.10,0.75,0.90), na.rm=T), df=3), arglag=list(fun='strata', breaks=1) , group = data_sub$year) 
      cb_rh <- onebasis(data_sub$RH, "ns", df=4)
      # onebasis for exposure - response model 
      cb_pol_var <- onebasis(runMean(data_sub[[pollutant]],0:3), "ns", df=3)
      # crossbasis for lag - response model 
      cb_pol_lag <- crossbasis(data_sub[[pollutant]], lag = 7, argvar = argvar, arglag = arglag, group = data_sub$year)
      
      # count holidays with keep is TRUE and keep factor(holidays) only when it is >1
      add_holidays <- ifelse(data_sub[keep == TRUE & holidays == 1, .N] > 1, "+ factor(holidays)", "")
      
      # run exposure and lag - response models
      model_var <- gnm(as.formula(paste(input, "~ cb_pol_var + spldoy:factor(year) + cb_temp + cb_rh + factor(dow)", add_holidays)), 
      eliminate = stratum, data = data_sub, family = quasipoisson, subset = keep)
      
      model_lag <- gnm(as.formula(paste(input, "~ cb_pol_lag + spldoy:factor(year) + cb_temp + cb_rh + factor(dow)", add_holidays)), 
      eliminate = stratum, data = data_sub, family = quasipoisson, subset = keep)
     
      # record overall estimates
      coef_vcov_list[[dist]] <- list()
      df_pol <- 3
      ind <- seq(df_pol)
      coefall <- coef(model_var)[ind]  
      vcovall <- vcov(model_var)[ind,ind] 
      coef_vcov_list[[dist]]$coefall <- coefall
      coef_vcov_list[[dist]]$vcovall <- vcovall
      
      # record lag-wise estimates 
      # crossreduce model by lag 
      mod_pred_red <- crossreduce(cb_pol_lag, model_lag, type="var", value = 10, cen=0) 
      
      # extract estimates
      coeflag <- coef(mod_pred_red)
      vcovlag <- vcov(mod_pred_red)
      coef_vcov_list[[dist]]$coeflag <- coeflag
      coef_vcov_list[[dist]]$vcovlag <- vcovlag
      
    }
    # attach district numbers to the dataframe
    names(coef_vcov_list) <- districts
    # filter coefs and vcovs for overall meta-analysis
    coefall <- do.call(rbind, lapply(coef_vcov_list, function(x) x$coefall))
    vcovall <- lapply(coef_vcov_list, function(x) x$vcovall)
    
    # Run overall meta-analysis 
    meta_tot <- mixmeta(coefall, vcovall , method="ml", random=~1|names(coef_vcov_list), control=list(igls.inititer=10))
    
    # Predict coef and vcov from overall meta-analysis
    max_pol <- ceiling(max(data[[pollutant]], na.rm = TRUE) / 5) * 5
    meta_pred_tot <- crosspred(do.call(onebasis, c(list(x=0:max_pol), list(fun="ns", df=3))), coef=coef(meta_tot), vcov=vcov(meta_tot), cen=0, model.link="log")
    
    # Store overall meta-analysis results in a list
    results_list[[input]][[pollutant]]$meta_tot <- meta_tot
    results_list[[input]][[pollutant]]$meta_pred_tot <- meta_pred_tot
    
    # Run lag-wise meta-analysis 
    coeflag <- do.call(rbind, lapply(coef_vcov_list, function(x) x$coeflag))
    vcovlag <- lapply(coef_vcov_list, function(x) x$vcovlag)
    meta_lag <- mixmeta(coeflag, vcovlag, method="ml" , random=~1|names(coef_vcov_list), control=list(igls.inititer=10))
    
    # do prediction of meta lag
    meta_pred_lag <- crosspred(do.call(onebasis, c(list(x=seq(0,7)), attr(cb_pol_lag,"arglag"))), coef=coef(meta_lag), vcov=vcov(meta_lag), model.link="log", at=0:7)
    
    # record and store effects
    for (lag in 0:7) {
      RR <- meta_pred_lag$matRRfit[lag + 1]
      ci_low <- meta_pred_lag$matRRlow[lag + 1]
      ci_high <- meta_pred_lag$matRRhigh[lag + 1]
      SE <- (log(ci_high) - log(ci_low)) / (2 * 1.96)
      Z <- log(RR) / SE
      pvalue <- 2 * (1 - pnorm(abs(Z)))
      effect <- data.frame(
        input = input,
        pollutant = pollutant,
        lag = lag,
        RR = RR,
        ci_low = ci_low,
        ci_high = ci_high,
        pvalue = pvalue
      )
      # Save individual results
      results_lag <- rbind(results_lag, effect)
    }
    results_list[[input]][[pollutant]]$results_lag <- results_lag
    results_list[[input]][[pollutant]]$meta_lag <- meta_lag
    results_list[[input]][[pollutant]]$meta_pred_lag <- meta_pred_lag
    all_results_lag <- rbind(all_results_lag, results_lag)
  }
}
# Save results_list to an RDS file
saveRDS(results_list, file = "2_step_results_100k.rds")

# Save all_results_lag to a CSV file
write.csv(all_results_lag, file = "2_step_results_100k.csv", row.names = FALSE)


# CHECK HETEROGENEITY STATS 
summary(results_list$all_CA$PM25$meta_tot)
summary(results_list$all_CA$PM10$meta_tot)
summary(results_list$all_CA$NO2$meta_tot)
summary(results_list$all_CA$O3$meta_tot)
summary(results_list$all_CA$SO2$meta_tot)
summary(results_list$all_CA$CO$meta_tot)