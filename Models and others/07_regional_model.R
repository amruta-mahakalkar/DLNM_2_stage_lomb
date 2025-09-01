# LOMB CA_CAMS AND PRE-PROCESS 

data_lomb <- read.csv('path to regional data') 
data_lomb <- as.data.table(data_lomb)
data_lomb$Date <- as.Date(data_lomb$Date)
data_lomb <- data_lomb[data_lomb$Date <= as.Date("2019-12-31"), ]

# Edit date format and add additional time related parameters
data_lomb$Date <- as.Date(data_lomb$Date)
data_lomb$year <- as.factor(substr(data_lomb$Date,1,4))
data_lomb$month <- as.factor(months(data_lomb$Date,abbr=TRUE))
data_lomb$day <- format(data_lomb$Date, format = "%d")
data_lomb$dow <- wday(data_lomb$Date) 
data_lomb$doy <- yday(data_lomb$Date)
# Update the "adult" column
data_lomb[, adult := adult + young]
# Generate holiday dates for Italy from 2016 to 2019
data_lomb$holidays <- as.integer(data_lomb$Date %in% as.Date(milan_holidays))
# define stratum
data_lomb[, stratum:=factor(paste(year, month, dow, sep=":"))] 
data_lomb[, RH_03 := frollmean(RH, n = 4, align = "right")]


pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO") 
inputs <- c("all_CA", "adult", "senior", "male", "female", "urban", "rural") 

results_list <- list()
l_results <- data.frame()
# Process through loops
for (input in inputs) {
  print(paste("Processing input:", input))
  for (pol in pollutants) {
    print(paste("Processing pollutant:", pol))
    coef_vcov_list <- list()
    results_lag <- data.frame()
    param_pol <- selected_params[selected_params$pol == pol, ]
    data_sub <- as.data.table(data_lomb)
    data_sub[, keep := sum(get(input)) > 0, by = stratum]
    
    # adjusting for confounders 
    temp_knots <- as.numeric(strsplit(param_pol$temp_var_knots, "\\|")[[1]])
    argvar_temp <- list(fun = "ns", knots = quantile(data_sub$Temp, temp_knots, na.rm = TRUE))
    arglag_temp <- if (param_pol$temp_lag_fun == "logknots") {list(knots=logknots(14, 2))
    } else {
      list(fun = "strata")}
    
    cb_temp <- crossbasis(data_sub$Temp, lag=14, argvar=argvar_temp, arglag=arglag_temp) 
    cb_rh <- onebasis(data_sub[[param_pol$rh_var]], "ns", df =  param_pol$rh_df)
    
    spldoy <- if (param_pol$spldoy_fun == "ns") {
      onebasis(data_sub$doy, "ns", df = param_pol$spldoy_df)
    } else {
      data_sub$doy}
    
    argvar_pol <- if (param_pol$pol_fun == "lin") {list(fun = "lin")
    } else {
      list(fun = "ns", df = param_pol$pol_df)}
    arglag_pol <- if (param_pol$pol_lag_fun == "ns") {list(fun = "ns", df = 3)
    } else {
      list(knots=logknots(7, 2))}
    cb_pol <- crossbasis(data_sub[[pol]], lag = 7, argvar = argvar_pol, arglag = arglag_pol) 
    
    # run exposure and lag - response models 
    model <- gnm(as.formula(paste(input, "~ cb_pol + spldoy + cb_temp + cb_rh + factor(holidays)")), eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep)
    
    # make prediction
    pol_value <- if (pol == "CO") {100} else if (pol == "SO2") {1} else {10}
    prediction <- crosspred(cb_pol, model, cen = 0, at = pol_value) 
    
    # Extract RR, ci_low and ci_high
    results_lag <- data.frame()
    for (lag in 0:7) {
      RR <- prediction$matRRfit[lag + 1]
      ci_low <- prediction$matRRlow[lag + 1]
      ci_high <- prediction$matRRhigh[lag + 1]
      SE <- (log(ci_high) - log(ci_low)) / (2 * 1.96)
      Z <- log(RR) / SE
      pvalue <- 2 * (1 - pnorm(abs(Z)))
      # Append the results to df
      results_lag <- rbind(results_lag, data.frame(
        input = paste0("lomb_", input),
        pollutant = pol,
        lag = lag,
        RR = RR,
        ci_low = ci_low,
        ci_high = ci_high,
        pvalue = pvalue,
        stringsAsFactors = FALSE
      ))
      
    }
    l_results <- rbind(l_results, cbind(results_lag))
    results_list[[paste0("lomb_", input)]][[pol]]$results_lag <- results_lag
  }
}
saveRDS(results_list, file = "results_list_regional.rds")
write.csv(l_results, "RR_results_regional.csv", row.names = FALSE)