# URBAN - RURAL 2-STAGE MODEL
# Define parameters 
districts <- unique(data$DISTRICT)
lags <- 0:7
argvar <- list(fun = "lin")
arglag <- list(knots=logknots(7, 2)) 

pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO") 
inputs <- c("urban", "rural")  

# Initialize a nested list to store results
results_list <- list()
all_results_lag <- data.frame()
# Process through multiple loops
for (input in inputs) {
  data_input <- data[data$urbn_rurl == input, ]
  districts <- as.character(unique(data_input$DISTRICT))
  print(paste("Processing input:", input))
  results_list[[input]] <- list()
  for (pollutant in pollutants) {
    print(paste("Processing pollutant:", pollutant))
    coef_vcov_list <- list()
    results_lag <- data.frame()
    for (dist in districts)  {
      data_sub <- data_input[data_input$DISTRICT == dist, ]
      data_sub[, keep := sum(get(input)) > 0, by = stratum]
      # accounting for long-term trend
      spldoy <- onebasis(data_sub$doy, "ns", df=3) 
      # adjusting for confounders 
      cb_temp <- crossbasis(data_sub$Temp, lag=14, argvar=list(fun="ns", knots=quantile(data_sub$Temp, c(0.10,0.75,0.90), na.rm=T), df=3), arglag=list(fun='strata', breaks=1) , group = data_sub$year) 
      cb_rh <- onebasis(data_sub$RH, "ns", df=4)
      # crossbasis for exposure and lag - response model 
      cb_pol_var <- crossbasis(data_sub[[pollutant]], lag = 4, argvar = argvar, arglag = arglag, group = data_sub$year) 
      
      cb_pol_lag <- crossbasis(data_sub[[pollutant]], lag = 7, argvar = argvar, arglag = arglag, group = data_sub$year)
      
      # run exposure and lag - response models 
      model_var <- gnm(as.formula(paste(input, "~ cb_pol_var + spldoy:factor(year) + cb_temp + cb_rh + factor(dow) + factor(holidays)")), eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep)
      
      model_lag <- gnm(as.formula(paste(input, "~ cb_pol_lag + spldoy:factor(year) + cb_temp + cb_rh + factor(dow) + factor(holidays)")), eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep)
      
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
    meta_tot <- mixmeta(coefall, vcovall, method="ml", random=~1|names(coef_vcov_list), control=list(igls.inititer=10)) 
    
    # Predict coef and vcov from overall meta-analysis
    max_pol <- ceiling(max(data[[pollutant]], na.rm = TRUE) / 5) * 5
    meta_pred_tot <- crosspred(do.call(onebasis, c(list(x=0:max_pol), list(fun="ns", df=3))), coef=coef(meta_tot), vcov=vcov(meta_tot), cen=0, model.link="log")
    
    # Store overall meta-analysis results in a list
    results_list[[input]][[pollutant]]$meta_tot <- meta_tot
    results_list[[input]][[pollutant]]$meta_pred_tot <- meta_pred_tot
    
    # Run lag-wise meta-analysis 
    coeflag <- do.call(rbind, lapply(coef_vcov_list, function(x) x$coeflag))
    vcovlag <- lapply(coef_vcov_list, function(x) x$vcovlag)
    meta_lag <- mixmeta(coeflag, vcovlag, method="ml", random=~1|names(coef_vcov_list), control=list(igls.inititer=10)) #method="ml" added
    
    # do pred of meta lag
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
saveRDS(results_list, file = "2_step_urb_rur_100k.rds")

# Save all_results_lag to a CSV file
write.csv(all_results_lag, file = "2_step_urb_rur_100k.csv", row.names = FALSE)
