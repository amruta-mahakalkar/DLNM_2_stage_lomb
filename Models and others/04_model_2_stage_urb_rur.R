# total OHCA, age and sex stratified by urban and rural setting
inputs <- c("all_CA", "adult", "senior", "male", "female") 

results_list <- list()
all_results_lag <- data.frame()
# Process through loops
for (class in unique(data$urbn_rurl)) {
  data_class <- data[data$urbn_rurl == class, ] 
  districts <- as.character(unique(data_class$DISTRICT))
  print(paste("Processing class:", class))
  for (pol in pollutants) {
    print(paste("Processing pollutant:", pol))
    
    # get best params
    coef_vcov_list <- list()
    param_pol <- selected_params[selected_params$pol == pol, ]
    
    results_list[[pol]][[class]] <- list()
    
    for (input in inputs) {
      print(paste("Processing input:", input))
      results_lag <- data.frame()
      for (dist in districts)  {
        data_sub <- data_class[data_class$DISTRICT == dist, ]
        data_sub[, keep := sum(get(input)) > 0, by = stratum]
        
        # adjusting for confounders 
        temp_knots <- as.numeric(strsplit(param_pol$temp_var_knots, "\\|")[[1]])
        argvar_temp <- list(fun = "ns", knots = quantile(data_sub$Temp, temp_knots, na.rm = TRUE))
        arglag_temp <- list(fun = "strata")
        cb_temp <- crossbasis(data_sub$Temp, lag=14, argvar=argvar_temp, arglag=arglag_temp) 
        cb_rh <- onebasis(data_sub[[param_pol$rh_var]], "ns", df =  param_pol$rh_df)
        
        # crossbasis for exposure and lag - response model 
        argvar_pol <- if (param_pol$pol_fun == "lin") {list(fun = "lin")
        } else {
          list(fun = "ns", df = param_pol$pol_df)}
        arglag_pol <- if (param_pol$pol_lag_fun == "ns") {list(fun = "ns", df = 3)
        } else {
          list(knots=logknots(7, 2))}
        
        cb_pol_lag <- crossbasis(data_sub[[pol]], lag = 7, argvar = argvar_pol, arglag = arglag_pol) 
        
        data_sub[, pol_lag_highest := shift(get(pol), n = pol_lag_exp_high[[pol]], type = "lag")] 
        cb_pol_var_high <- do.call(onebasis, c(list(x = data_sub$pol_lag_highest), argvar_pol))
        
        spldoy <- if (param_pol$spldoy_fun == "ns") {
          onebasis(data_sub$doy, "ns", df = param_pol$spldoy_df)
        } else {
          data_sub$doy}
        
        # run exposure and lag - response models 
        model_var <- gnm(as.formula(paste(input, "~ cb_pol_var_high + cb_temp + cb_rh + factor(holidays)")), eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep) # 
        model_lag <- gnm(as.formula(paste(input, "~ cb_pol_lag + cb_temp + cb_rh + factor(holidays)")), eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep) 
        
        # record overall estimates
        coef_vcov_list[[dist]] <- list()
        df_pol <- param_pol$pol_df
        ind <- if (df_pol > 0) seq_len(df_pol) else 1
        coefall <- coef(model_var)[ind]  
        vcovall <- vcov(model_var)[ind,ind] 
        coef_vcov_list[[dist]]$coefall <- coefall
        coef_vcov_list[[dist]]$vcovall <- vcovall
        
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
      # filter coefs and vcovs for overall meta-analysis
      coefall <- do.call(rbind, lapply(coef_vcov_list, function(x) x$coefall))
      vcovall <- lapply(coef_vcov_list, function(x) x$vcovall)
      
      # Run overall meta-analysis 
      meta_tot <- mixmeta(coefall, vcovall , method="ml", random=~1|names(coef_vcov_list), control=list(igls.inititer=10))
      
      # Predict coef and vcov from overall meta-analysis
      max_pol <- ceiling(max(data[[pol]], na.rm = TRUE) / 5) * 5
      meta_pred_tot <- crosspred(do.call(onebasis, c(list(x=0:max_pol), argvar_pol)), 
                                 coef=coef(meta_tot), vcov=vcov(meta_tot), cen=0, model.link="log") 
      
      # Store overall meta-analysis results in a list
      results_list[[pol]][[class]][[input]]$meta_tot <- meta_tot
      results_list[[pol]][[class]][[input]]$meta_pred_tot <- meta_pred_tot
      
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
          class = class,
          input = paste0(class, "_", input), 
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
      results_list[[pol]][[class]][[input]]$results_lag <- results_lag
      results_list[[pol]][[class]][[input]]$meta_lag <- meta_lag
      results_list[[pol]][[class]][[input]]$meta_pred_lag <- meta_pred_lag
      all_results_lag <- rbind(all_results_lag, results_lag)
    }
  }
}
# Save results
saveRDS(results_list, file = "result_list_urb_rur.rds")
write.csv(all_results_lag, file = "RR_results_urb_rur.csv", row.names = FALSE)
