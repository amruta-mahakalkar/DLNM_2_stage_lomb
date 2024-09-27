# PARAMETRES SELECTION
# Import libraries 
library(dlnm)
library(epiDisplay); library(lubridate) ; library(dplyr); library(tidyr); library(plyr); library(scales)
library(splines); library(gnm); library(Epi); library(tsModel);  library(sf) 
library(data.table); library(viridis)
library(mgcv) ; library(corrplot); library(mixmeta)
library(ggplot2) ; library(patchwork)

#load 100k dist file
data_dist <- read.csv('CAMS_CA_100k_v1.csv')
data_dist <- as.data.table(data_dist)

# DATA PRE-PROCESSING
# chose data only till 2019 
data_dist$Date <- as.Date(data_dist$Date)
data_dist <- data_dist[data_dist$Date <= as.Date("2019-12-31"), ]
# Edit date format and add additional time related parameters
data_dist$Date <- as.Date(data_dist$Date)
data_dist$year <- as.factor(substr(data_dist$Date,1,4))
data_dist$month <- as.factor(months(data_dist$Date,abbr=TRUE))
data_dist$day <- format(data_dist$Date, format = "%d")
data_dist$dow <- wday(data_dist$Date) 
data_dist$doy <- yday(data_dist$Date)
# Update the "adult" column
data_dist[, adult := adult + young]
# Generate holiday dates for Italy from 2016 to 2019
data_dist$holidays <- as.integer(data_dist$Date %in% as.Date(milan_holidays))
# define stratum
data_dist[, stratum:=factor(paste(DISTRICT, year, month, dow, sep=":"))]

# MAKE QAIC FUNCTION
fqaic_single <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

## PARAMETER SELECTION FOR TIME CURVE 
# Define the functions and degrees of freedom
fun <- c("ns") 
df <- 3:10
inputs <- c("all_CA", "senior", "adult", "male", "female")
# make storages
spydoy_qaic <- list()
spydoy_rmse <- list()
# Loop through each combination of fun and df
for (input in inputs) {
  for (f in fun) {
    for (d in df) {
      # make the name
      name <- paste0("spldoy_", input, "_", d)
      
      data_dist[, keep := sum(get(input)) > 0, by = stratum]
      
      # make the basis function
      spldoy <- onebasis(data_dist$doy, fun = f, df = d)
      
      model <- gnm(as.formula(paste(input, "~ spldoy")), eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
      
      # get their Q-AIC
      fqaic <- fqaic_single(model)
      rmse <- sqrt(mean(model$residuals^2))
      
      # Store the qaic values in a list 
      spydoy_qaic[[name]] <- fqaic
      spydoy_rmse[[name]] <- rmse
    }
  }
}

write.csv(spydoy_qaic, "spydoy_qaic.csv", row.names =TRUE)


## QAIC for argvar (fun, df) for temperature and RH
# Define the parameters
knots_argvar <- list(c(10, 90)/100, c(25, 75)/100, c(10, 75, 90)/100, c(25, 75, 90)/100, c(33.33, 66.67)/100)
lags_temp <- c(0,14)
df_arg <- 3:6
vars <- c("Temp", "RH") 

# make storage
argvar_dfs_qaic <- list()
for (var in vars) {
  # Loop through each combination
  argvar_dfs_qaic[[var]] <- list()
  
  # try linear argvar
  cb_lin <- onebasis(data_dist[[var]], fun = "lin")
  
  # Fit the linear model
  model_lin <- gnm(all_CA ~ cb_lin, eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
  
  # Calculate the Q-AIC for the linear model
  fqaic_lin <- fqaic_single(model_lin)
  argvar_dfs_qaic[[var]][[paste0("qaic_", "lin")]] <- fqaic_lin 
  
  # try non-linear argvar with dfs 
  for (df_var in df_arg) {
    cb_df <- onebasis(data_dist[[var]], fun = "ns", df = df_var)
    
    # Fit the model for the current df
    model_df <- gnm(all_CA ~ cb_df, eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
    
    fqaic_df <- fqaic_single(model_df)
    
    # Store the Q-AIC value for the current df model
    argvar_dfs_qaic[[var]][[paste0("qaic_", "ns", df_var)]] <- fqaic_df
  }
  # try non-linear argvar with dfs + specific knots 
  for (df_var in df_arg) {      
    for (knot in knots_argvar) {
      cb_knot <- onebasis(data_dist[[var]], fun = "ns", knots=quantile(data_dist[[var]], knot, na.rm=T), df = df_var)
      
      model_knot <- gnm(all_CA ~ cb_knot, eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)  
      # Calculate the Q-AIC
      fqaic_knots <- fqaic_single(model_knot)
      # Store the Q-AIC values in a list
      argvar_dfs_qaic[[var]][[paste0("qaic_", "ns_", df_var, "_knots_", paste(knot, collapse = "_"))]] <- fqaic_knots
    }
  }
}

## QAIC for arglag (fun, df) per pollutant
# Define the parameters
fun_arglag <- c("ns", "bs", "strata")
df_arglag <- 3:6
vars <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO")
# make storage
arglag_dfs_qaic <- list()
for (var in vars) {
  # Loop through each combination
  arglag_dfs_qaic[[var]] <- list()
  # try linear varlag
  cb_lin <- crossbasis(data_dist[[var]], lag = 7, argvar = list(fun = "ns", df = 3), arglag = list(fun = "lin"))
  
  # Fit the linear model
  model_lin <- gnm(all_CA ~ cb_lin, eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
  
  # Calculate the Q-AIC for the linear model
  fqaic_lin <- fqaic_single(model_lin)
  arglag_dfs_qaic[[var]][[paste0("qaic_", "lin")]] <- fqaic_lin 
  
  # try strata varlag
  cb_stra <- crossbasis(data_dist[[var]], lag = 7, argvar = list(fun = "ns", df = 3), arglag = list(fun = "strata", breaks=1))
  
  # Fit the strata model
  model_stra <- gnm(all_CA ~ cb_stra, eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
  
  # Calculate the Q-AIC for the strata
  fqaic_stra <- fqaic_single(model_stra)
  arglag_dfs_qaic[[var]][[paste0("qaic_", "stra")]] <- fqaic_stra
  
  # try non-linear varlag with functions + dfs
  for (fun_lag in fun_arglag) {
    for (df_lag in df_arglag) {
      cb_df <- crossbasis(data_dist[[var]], lag = 7, argvar = list(fun = "ns", df = 3),arglag = list(fun = fun_lag, df = df_lag))
      
      # Fit the model for the current df
      model_df <- gnm(all_CA ~ cb_df, eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
      
      fqaic_df <- fqaic_single(model_df)
      
      # store qaic
      arglag_dfs_qaic[[var]][[paste0("qaic_", fun_lag, df_lag)]] <- fqaic_df
    }
  }
}
unlisted_arglag <- unlist(var_lag_qaic)
write.csv(unlisted_arglag, "arglag_pol.csv", row.names = TRUE)


## QAIC for arglag for temp
vars <- c("Temp")
lags <- 1:14
knot_var <- list("Temp" = c(0.10, 0.75, 0.90)) 
fun = c("ns", "bs", "strata") 
knot_lag <- 1:3
var_lag_qaic <- list()
data_dist[, keep := sum(get(input)) > 0, by = stratum]
spldoy <- onebasis(data_dist$doy, "ns", df=3)

for (var in vars) {
  var_lag_qaic[[var]] <- list()
  # use best argvar from earlier analysis
  argvar = list(fun = "ns", knots = quantile(data_dist[[var]], knot_var[[var]], na.rm=TRUE), df = 3)
  # try varlag with lags, function and knots 
  for (lag in lags) {
    for (f in fun){
      for (k in knot_lag){
        # try strata with breaks = 1
        cb_lag_break <- crossbasis(data_dist[[var]], lag = lag, argvar = argvar, arglag = list(fun = "strata", breaks = 1))
        # trying other combinations
        cb_lag_knots <- crossbasis(data_dist[[var]], lag = lag, argvar = argvar, arglag = list(knots=logknots(lag, fun = f, k)))
        
        # Fit the linear model for cb_lag_break
        model_lag_break <- gnm(all_CA ~ cb_lag_break + spldoy:factor(year) + factor(dow), eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
        fqaic_lag_break <- fqaic_single(model_lag_break)
        var_lag_qaic[[var]][[paste0("qaic_break_", lag)]] <- fqaic_lag_break
        
        # Fit the linear model for cb_lag_knots
        model_lag_knots <- gnm(all_CA ~ cb_lag_knots + spldoy:factor(year) + factor(dow), eliminate = stratum, data = data_dist, family = quasipoisson, subset = keep)
        fqaic_lag_knots <- fqaic_single(model_lag_knots)
        var_lag_qaic[[var]][[paste0("qaic_knots_", lag, "_", f, "_", k)]] <- fqaic_lag_knots
      }
    }
  }
}
unlisted_arglag <- unlist(var_lag_qaic)
write.csv(unlisted_arglag, "arglag_temp.csv", row.names = TRUE)