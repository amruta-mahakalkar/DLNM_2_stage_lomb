# Import libraries 
library(dplyr); library(tidyr); library(splines); library(gnm); 
library(data.table); library(dlnm); library(lubridate)

#load regional file
data_dist <- read.csv('file path')
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
data_dist[, stratum:=factor(paste(year, month, dow, sep=":"))]  
data_dist[, RH_03 := frollmean(RH, n = 4, align = "right")] 

# Construct qaic formula
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

# Built param grid vectors
# doy spline paramters
spldoy_df <- 2:3 
spldoy_fun <- list("lin", "ns") 
# temp variables and lag parameters
temp_var_knots <- list(c(0.10, 0.75, 0.90), c(0.10, 0.90)) 
temp_lag_fun <- c("strata", "logknots") 
# relative humidity variable params
rh_df <- 3:4
rh_var <- c("RH", "RH_03") 
# define pollution params
pol_fun <- c("lin", "ns")
pol_lag <- c(3, 7, 10, 14)
pol_df <- 2:4 
pol_lag_fun <- c("ns", "logknots") 

param_grid <- expand.grid(spldoy_df = spldoy_df, 
                          spldoy_fun = spldoy_fun,
                          temp_var_knots = temp_var_knots, 
                          temp_lag_fun = temp_lag_fun,
                          rh_df = rh_df, 
                          rh_var = rh_var,
                          pol_fun = pol_fun, 
                          pol_df = pol_df,
                          pol_lag = pol_lag, 
                          pol_lag_fun = pol_lag_fun,
                          stringsAsFactors = FALSE)

param_grid$pol_df[param_grid$pol_fun == "lin"] <- NA
param_grid$spldoy_df[param_grid$spldoy_fun == "lin"] <- NA
param_grid <- unique(param_grid)

############################ 
# PARAMETER SELECTION FOR TIME CURVE, TEMP AND RH
pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO")
results_list <- list()
best_param_df <- list()
results_all_df <- list()
for (pol in pollutants) {
  print(paste("Processing pollutant:", pol))
  
  data_sub <- data_dist[, .(all_CA, get(pol), doy, year, Temp, RH, RH_03, dow, holidays, stratum)] 
  setnames(data_sub, "V2", pol)
  
  # container for results
  all_results_dt <- list()
  
  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]
    res <- tryCatch({
      data_sub[, keep := sum(all_CA) > 0, by = stratum]
      
      spldoy <- if (params$spldoy_fun == "ns") {
        onebasis(data_sub$doy, "ns", df = params$spldoy_df)
      } else {
        data_sub$doy  
      }
      
      # Build temperature
      temp_var_knots_unlisted <- unlist(params$temp_var_knots)
      argvar_temp <- list(fun = "ns", knots = quantile(data_sub$Temp, temp_var_knots_unlisted, na.rm = TRUE))
      arglag_temp <- switch(params$temp_lag_fun,
                            "strata" = list(fun = "strata"),
                            "logknots" = list(knots=logknots(14, 2)))
      
      cb_temp <- crossbasis(data_sub$Temp, lag = 14, argvar = argvar_temp, arglag = arglag_temp) 
      
      # Build RH
      cb_rh <- onebasis(data_sub[[params$rh_var]], "ns", df = params$rh_df) 
      
      # Build Pollutant
      argvar_pol <- switch(params$pol_fun,
                           "ns" = list(fun = "ns", df = params$pol_df),
                           "lin" = list(fun = "lin"))
      
      arglag_pol<- switch(params$pol_lag_fun,
                          "ns" = list(fun = "ns", df = 3),
                          "logknots" = list(knots=logknots(params$pol_lag, 2)))
      
      cb_pol <- crossbasis(data_sub[[pol]], lag=params$pol_lag, argvar=argvar_pol, arglag=arglag_pol)
      
      # Build Model
      model_formula <- as.formula(paste("all_CA ~", paste(c(
        "cb_pol", "spldoy", "cb_temp", "cb_rh", "factor(dow)", "factor(holidays)"
      ), collapse = " + ")))
      
      model_pol <- gnm(formula = model_formula, eliminate = stratum, data = data_sub, family = quasipoisson, subset = keep)
      
      # Predictions
      pol_value <- if (pol == "SO2") 1 else 10
      
      pred_pol <- crosspred(cb_pol, model_pol, cumul = TRUE, cen = 0, at = pol_value)
      
      # Build one row with param info + QAIC + lag-wise columns
      # Lag RR columns
      lag_fit  <- setNames(as.list(as.numeric(pred_pol$matRRfit)),  
                           paste0("lagRR_", 0:(ncol(pred_pol$matRRfit)-1), "_fit"))
      lag_low  <- setNames(as.list(as.numeric(pred_pol$matRRlow)),  
                           paste0("lagRR_", 0:(ncol(pred_pol$matRRlow)-1), "_low"))
      lag_high <- setNames(as.list(as.numeric(pred_pol$matRRhigh)), 
                           paste0("lagRR_", 0:(ncol(pred_pol$matRRhigh)-1), "_high"))
      
      # Cumulative RR columns
      cum_fit  <- setNames(as.list(as.numeric(pred_pol$cumRRfit)),  
                           paste0("cumRR_", 0:(ncol(pred_pol$cumRRfit)-1), "_fit"))
      cum_low  <- setNames(as.list(as.numeric(pred_pol$cumRRlow)),  
                           paste0("cumRR_", 0:(ncol(pred_pol$cumRRlow)-1), "_low"))
      cum_high <- setNames(as.list(as.numeric(pred_pol$cumRRhigh)), 
                           paste0("cumRR_", 0:(ncol(pred_pol$cumRRhigh)-1), "_high"))
      
      # Combine all into one row
      data.table(pol = pol, QAIC = fqaic(model_pol))[, (names(params)) := params
      ][, (names(lag_fit))  := lag_fit
      ][, (names(lag_low))  := lag_low
      ][, (names(lag_high)) := lag_high
      ][, (names(cum_fit))  := cum_fit
      ][, (names(cum_low))  := cum_low
      ][, (names(cum_high)) := cum_high]
      
    }, error = function(e) {
      # Empty row if error
      data.table(pol = pol, QAIC = Inf)[
        , (names(params)) := params
      ]
    })
    
    all_results_dt[[i]] <- res
  }
  
  # Combine to one table per pollutant
  all_results_dt <- rbindlist(all_results_dt, fill = TRUE)
  results_list[[pol]] <- all_results_dt
  
  # Pick best QAIC row
  best_row <- all_results_dt[which.min(QAIC)]
  best_param_df[[pol]] <- best_row
}

# Final best param table
best_param_df <- rbindlist(best_param_df, fill = TRUE)
results_all_df <- rbindlist(results_list, fill = TRUE)

fwrite(best_param_df, file = "best_params.csv", row.names = FALSE)
fwrite(results_all_df, sfile = "params_qaic_rr.csv", row.names = FALSE)

# TEST DOW AND SPLDOY
all_results <- data.frame()   
results_list <- list() 

for (pol in pollutants) {
  print(paste("Processing pollutant:", pol))
  param_pol <- best_params_df[best_params_df$pol == pol, ]
  
  data_sub <- data_dist
  data_sub[, keep := sum(all_CA) > 0, by = stratum]
  
  spldoy <- if (param_pol$spldoy_fun == "ns") {
    onebasis(data_sub$doy, "ns", df = param_pol$spldoy_df)
  } else {
    data_sub$doy  
  }
  
  # Built temperature, rh and pol from best params
  temp_knots <- as.numeric(strsplit(param_pol$temp_var_knots, "\\|")[[1]])
  argvar_temp <- list(fun = "ns", knots = quantile(data_sub$Temp, temp_knots, na.rm = TRUE))
  arglag_temp <-  list(knots=logknots(14, 2))
  
  cb_temp <- crossbasis(data_sub$Temp, lag=14, argvar=argvar_temp, arglag=arglag_temp, group = data_sub$year) 
  
  cb_rh <- onebasis(data_sub[[param_pol$rh_var]], "ns", df =  param_pol$rh_df)
  
  argvar_pol <- if (param_pol$pol_fun == "lin") {list(fun = "lin")} else {list(fun = "ns", df = param_pol$pol_df)}
  arglag_pol <- if (param_pol$pol_lag_fun == "ns") {list(fun = "ns", df = 3)} else {list(knots=logknots(7, 2))}
  
  cb_pol <- crossbasis(data_sub[[pol]], lag = 7, argvar = argvar_pol, arglag = arglag_pol, group = data_sub$year)
  
  # Model formula
  model_formula <- as.formula(paste("all_CA ~", paste(c("cb_pol", "spldoy", "cb_temp",
                                                        "cb_rh", "factor(dow)", "factor(holidays)"),
                                                      collapse = " + ")))
  
  model_1 <- gnm(formula = model_formula, eliminate = stratum, data = data_sub, family = quasipoisson,
                 subset = keep)
  model_2  <- update(model_1, . ~ . - spldoy)  
  model_3  <- update(model_1, . ~ . - factor(dow)) 
  model_4  <- update(model_1, . ~ . - factor(holidays))
  model_5  <- update(model_1, . ~ . - factor(holidays) - factor(dow))
  model_6  <- update(model_1, . ~ . - spldoy - factor(dow))
  model_7  <- update(model_1, . ~ . - spldoy - factor(dow) - factor(holidays)) 
  
  # Make predictions
  pol_value <- if (pol == "CO") {100} else if (pol == "SO2") {1} else {10}
  
  preds <- list(
    model_1 = crosspred(cb_pol, model_1, cumul = TRUE, cen = 0, at = pol_value),
    model_2 = crosspred(cb_pol, model_2, cumul = TRUE, cen = 0, at = pol_value),
    model_3 = crosspred(cb_pol, model_3, cumul = TRUE, cen = 0, at = pol_value),
    model_4 = crosspred(cb_pol, model_4, cumul = TRUE, cen = 0, at = pol_value),
    model_5 = crosspred(cb_pol, model_5, cumul = TRUE, cen = 0, at = pol_value),
    model_6 = crosspred(cb_pol, model_6, cumul = TRUE, cen = 0, at = pol_value),
    model_7 = crosspred(cb_pol, model_7, cumul = TRUE, cen = 0, at = pol_value),
  )
  
  # Store results for each model
  for (i in seq_along(preds)) {
    prediction <- preds[[i]]
    model_name <- names(preds)[i]
    model_obj <- get(model_name)
    
    results_lag <- data.frame()
    
    for (lag in 0:7) {
      RR <- prediction$matRRfit[lag + 1]
      ci_low <- prediction$matRRlow[lag + 1]
      ci_high <- prediction$matRRhigh[lag + 1]
      SE <- (log(ci_high) - log(ci_low)) / (2 * 1.96)
      Z <- log(RR) / SE
      pvalue <- 2 * (1 - pnorm(abs(Z)))
      qaic <- fqaic(model_obj)
      
      results_lag <- rbind(results_lag, data.frame(
        input = "all_CA",
        pollutant = pol,
        model = model_name,
        lag = lag,
        RR = RR,
        ci_low = ci_low,
        ci_high = ci_high,
        pvalue = pvalue,
        qaic = qaic,
        stringsAsFactors = FALSE
      ))
    }
    # Store results
    results_list[[pol]] <- results_lag
    all_results <- rbind(results_lag, all_results)
  }
}
fwrite(all_results, file = "time_params_qaic_rr.csv", row.names = FALSE)
