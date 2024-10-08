# LOMB CA_CAMS AND PRE-PROCESS 
#load regional file
data_lomb <- read.csv('CAMS_CA_lomb_v1.csv') 
data_lomb <- as.data.table(data_lomb)
# chose data only till 2019 
data_lomb$Date <- as.Date(data_lomb$Date)
data_lomb <- data_lomb[data_lomb$Date <= as.Date("2019-12-31"), ]
# Edit date format and add additional time related parameters
data_lomb$Date <- as.Date(data_lomb$Date)
data_lomb$year <- as.factor(substr(data_lomb$Date,1,4))
data_lomb$month <- as.factor(months(data_lomb$Date,abbr=TRUE))
data_lomb$day <- format(data_lomb$Date, format = "%d")
data_lomb$dow <- wday(data_lomb$Date) 
data_lomb$doy <- yday(data_lomb$Date) 
data_lomb$holy <- as.integer(data_lomb$Date %in% as.Date(milan_holidays))
# Update the "adult" column
data_lomb[, adult := adult + young]
# define stratum
data_lomb[, stratum:=factor(paste(year, month, dow, sep=":"))]

# CB - CONFOUNDERS 
# make spline of time to account for long-term exposure 
spldoy_l <- onebasis(data_lomb$doy, "ns", df=7)
# spline for temp and RH 
cb_temp_l <- crossbasis(data_lomb$Temp, lag=14, argvar=list(fun="ns", knots=quantile(data_lomb$Temp, c(0.10,0.75,0.90), na.rm=T), df=3), arglag=list(fun='strata', breaks=1), group = data_lomb$year)
cb_rh_l <- onebasis(data_lomb$RH, "ns", df=4)

# CB - POLLUTANTS AND MODEL 
pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO")
inputs <- c("all_CA", "adult", "senior", "male", "female", "urban", "rural")

argvar <- list(fun = "lin")
arglag <- list(fun = "ns", df=3)

# Initialize a data frame to store the results
l_results <- data.frame()
# Initialize lists for models and predictions
for (input in inputs) {
  assign(paste0("l_models_", input), list())
  assign(paste0("l_predictions_", input), list())
  
  data_lomb[, keep := sum(get(input)) > 0, by = stratum]
  
  for (pollutant in pollutants) {
    # make pol cb
    cb_pol <- crossbasis(data_lomb[[pollutant]], lag = 7, argvar = argvar, arglag = arglag, group = data_lomb$year)
    
    # make model
    model <- gnm(as.formula(paste(input, "~ cb_pol + spldoy_l:factor(year) + cb_temp_l + cb_rh_l + factor(dow) + factor(holy)")), eliminate = stratum, data = data_lomb, family = quasipoisson, subset = keep)
    # make prediction
    prediction <- crosspred(cb_pol, model, cen = 0, at = 10) 
    # Extract RR, ci_low and ci_high
    for (lag in 0:7) {
      RR <- prediction$matRRfit[lag + 1]
      ci_low <- prediction$matRRlow[lag + 1]
      ci_high <- prediction$matRRhigh[lag + 1]
      SE <- (log(ci_high) - log(ci_low)) / (2 * 1.96)
      Z <- log(RR) / SE
      pvalue <- 2 * (1 - pnorm(abs(Z)))
      # Append the results to the data frame
      l_results <- rbind(l_results, data.frame(
        input = input,
        pollutant = pollutant,
        lag = lag,
        RR = RR,
        ci_low = ci_low,
        ci_high = ci_high,
        pvalue = pvalue,
        stringsAsFactors = FALSE
      ))
    }
    models_list <- get(paste0("l_models_", input))
    models_list[[pollutant]] <- model
    assign(paste0("l_models_", input), models_list)
    predictions_list <- get(paste0("l_predictions_", input))
    predictions_list[[pollutant]] <- prediction
    assign(paste0("l_predictions_", input), predictions_list)
  }
}
write.csv(l_results, "RR_results_lombardia.csv", row.names = FALSE)
