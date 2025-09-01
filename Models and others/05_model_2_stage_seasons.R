# total OHCA stratified by seasons

pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO") 
districts <- as.character(unique(data$DISTRICT))

knot_temp <- list(
  Winter = c(0.1, 0.25),
  Autumn = c(0.1, 0.25),
  Summer = c(0.75, 0.90),
  Spring = c(0.75, 0.90))

lag_temp <- list(Winter = 14, Autumn = 14, Summer=14, Spring=14) 
seasons <- unique(data$Season)

results_list <- list()
all_results_lag <- data.frame()
# Process through loops
for (season in seasons) {
  data_input <- data[data$Season == season, ]
  print(paste("Processing season:", season))
  results_list[[season]] <- list()
  for (pol in pollutants) {
    print(paste("Processing pollutant:", pol))
    param_pol <- selected_params[selected_params$pol == pol, ]
    coef_vcov_list <- list()
    results_lag <- data.frame()
    for (dist in districts)  {
      data_sub <- data_input[data_input$DISTRICT == dist]
      data_sub[, keep := sum(all_CA)>0, by = stratum]
      
      # adjusting for confounders 
      temp_knots <- knot_temp[[season]]
      argvar_temp <- list(fun = "ns", knots = quantile(data_sub$Temp, temp_knots, na.rm = TRUE))
      arglag_temp <- if (param_pol$temp_lag_fun == "logknots") {list(knots=logknots(14, 2))
      } else {
        list(fun = "strata")}
      cb_temp <- crossbasis(data_sub$Temp, lag=lag_temp[[season]], argvar=argvar_temp, arglag=arglag_temp) 
      cb_rh <- onebasis(data_sub[[param_pol$rh_var]], "ns", df =  param_pol$rh_df)
      
      # crossbasis for exposure and lag - response model 
      argvar_pol <- if (param_pol$pol_fun == "lin") {list(fun = "lin")
      } else {
        list(fun = "ns", df = param_pol$pol_df)}
      arglag_pol <- if (param_pol$pol_lag_fun == "ns") {list(fun = "ns", df = 3)
      } else {
        list(knots=logknots(7, 2))}
      
      cb_pol_lag <- crossbasis(data_sub[[pol]], lag = 7, argvar = argvar_pol, arglag = arglag_pol) 
      
      spldoy <- if (param_pol$spldoy_fun == "ns") {
        onebasis(data_sub$doy, "ns", df = param_pol$spldoy_df)
      } else {
        data_sub$doy}
      
      tryCatch({
        model_lag <- gnm(all_CA ~ cb_pol_lag + cb_temp + cb_rh + factor(holidays),
                         eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep) #  
      }, error = function(e) {
        print(paste("Error processing pollutant:", pol, "in district:", dist, ":", e$message))
      })
      
      # record lag-wise estimates 
      # crossreduce model by lag 
      pol_value <- if (pol == "CO") {100} else if (pol == "SO2") {1} else {10}
      mod_pred_red <- crossreduce(cb_pol_lag, model_lag, type="var", value = pol_value, cen=0) 
      
      # extract estimates
      coeflag <- coef(mod_pred_red)
      vcovlag <- vcov(mod_pred_red)
      coef_vcov_list[[dist]]$coeflag <- coeflag
      coef_vcov_list[[dist]]$vcovlag <- vcovlag
      
    }
    # attach district numbers to the dataframe
    names(coef_vcov_list) <- districts
    
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
        input = season,
        pollutant = pol,
        lag = lag,
        RR = RR,
        ci_low = ci_low,
        ci_high = ci_high,
        pvalue = pvalue
      )
      # Save individual results
      results_lag <- rbind(results_lag, effect)
    }
    results_list[[season]][[pol]]$results_lag <- results_lag
    results_list[[season]][[pol]]$meta_lag <- meta_lag
    results_list[[season]][[pol]]$meta_pred_lag <- meta_pred_lag
    all_results_lag <- rbind(all_results_lag, results_lag)
  }
}
# Save results
saveRDS(results_list, file = "results_list_seasons.rds")
write.csv(all_results_lag, file = "RR_results_seasons.csv", row.names = FALSE)
