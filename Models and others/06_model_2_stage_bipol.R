# Bi-pol models of total OHCA by interaction and confounding

pollutants <- c("PM25", "O3","CO", "NO2", "SO2") 
pol_cat_cols <- paste0(pollutants, "_cat") 
pol_cats <- c("Low", "High")

districts <- unique(data$DISTRICT)

results_list <- list()
coef_vcov_list <- list() 
all_results_lag <- data.frame()
for (cat_col in pol_cat_cols) {
  pol_1 <- sub("_cat$", "", cat_col)
  print(paste("processing interaction pol:", pol_1))
  # choose pols excluding current one
  run_pollutants <- setdiff(pollutants, pol_1)
  
  for (pol_2 in run_pollutants) {
    print(paste("processing pol:", pol_2))
    param_pol1 <- selected_params[selected_params$pol == pol_1, ]
    for (dist in districts)  {
      coef_vcov_list[[dist]] <- list() 
      data_sub <- data[data$DISTRICT == dist, ]
      data_sub[, keep := sum(all_CA) > 0, by = stratum]
      
      # adjusting for confounders 
      temp_knots <- as.numeric(strsplit(param_pol1$temp_var_knots, "\\|")[[1]])
      argvar_temp <- list(fun = "ns", knots = quantile(data_sub$Temp, temp_knots, na.rm = TRUE))
      arglag_temp <- if (param_pol1$temp_lag_fun == "logknots") {list(knots=logknots(14, 2))
      } else {
        list(fun = "strata")}
      
      cb_temp <- crossbasis(data_sub$Temp, lag=14, argvar=argvar_temp, arglag=arglag_temp) 
      cb_rh <- onebasis(data_sub[[param_pol1$rh_var]], "ns", df =  param_pol1$rh_df)
      
      argvar_pol1 <- if (param_pol1$pol_fun == "lin") {list(fun = "lin")
      } else {
        list(fun = "ns", df = param_pol1$pol_df)}
      arglag_pol1 <- if (param_pol1$pol_lag_fun == "ns") {list(fun = "ns", df = 3)
      } else {
        list(knots=logknots(7, 2))}
      
      cb_pol1 <- crossbasis(data_sub[[pol_1]], lag = 7, argvar = argvar_pol1, arglag = arglag_pol1) 
      
      # get interaction terms
      cb_pol2 <- crossbasis(data_sub[[pol_2]], lag = 7, argvar = argvar_pol1, arglag = arglag_pol1, group = data_sub$year)
      
      cb_pol_int  <- cb_pol1*cb_pol2
      
      model_main <- gnm(all_CA ~ cb_pol1 + cb_temp + cb_rh + factor(holidays), eliminate = stratum, data = data_sub, family = quasipoisson, subset=keep) 
      model_int <- update(model_main, . ~ . + cb_pol_int) 
      model_plus <- update(model_main, . ~ . + cb_pol2) 
      
      pol_value <- if (pol_1 == "CO") {100} else if (pol_1 == "SO2") {1} else {10}
      mod_pred_int  <- crossreduce(cb_pol1, model_int, value=pol_value, type="var", cen=0)
      mod_pred_plus  <- crossreduce(cb_pol1, model_plus, value=pol_value, type="var", cen=0)
      
      coeflag1 <- coef(mod_pred_int)
      vcovlag1 <- vcov(mod_pred_int)
      coeflag2 <- coef(mod_pred_plus)
      vcovlag2 <- vcov(mod_pred_plus)
      coef_vcov_list[[dist]]$coeflag1 <- coeflag1
      coef_vcov_list[[dist]]$vcovlag1 <- vcovlag1
      coef_vcov_list[[dist]]$coeflag2 <- coeflag2
      coef_vcov_list[[dist]]$vcovlag2 <- vcovlag2
    }
    # attach district numbers to the dataframe
    names(coef_vcov_list) <- districts
    # filter coefs and vcovs for overall meta-analysis
    coeflag1 <- do.call(rbind, lapply(coef_vcov_list, function(x) x$coeflag1))
    vcovlag1 <- lapply(coef_vcov_list, function(x) x$vcovlag1)
    coeflag2 <- do.call(rbind, lapply(coef_vcov_list, function(x) x$coeflag2))
    vcovlag2 <- lapply(coef_vcov_list, function(x) x$vcovlag2)
    
    # Do meta-analysis 
    meta_lag_int <- mixmeta(coeflag1, vcovlag1, method="ml", random = ~1 | names(coef_vcov_list), control = list(igls.inititer = 10))
    meta_lag_plus <- mixmeta(coeflag2, vcovlag2, method="ml", random = ~1 | names(coef_vcov_list), control = list(igls.inititer = 10))
    
    # Predict lag-response from meta-analysis
    lag_basis <- do.call(onebasis, c(list(x = 0:7), attr(cb_pol1,"arglag")))
    meta_pred_lag <- crosspred(lag_basis, coef = coef(meta_lag_int), vcov = vcov(meta_lag_int), model.link = "log", at = 0:7)
    meta_pred_plus <- crosspred(lag_basis, coef = coef(meta_lag_plus), vcov = vcov(meta_lag_plus), 
                                model.link = "log", at = 0:7)
    meta_preds <- list(main = meta_pred_lag, plus = meta_pred_plus)
    for (pred_name in names(meta_preds)) {
      results_lag <- data.frame()
      pred <- meta_preds[[pred_name]]
      for (lag in 0:7) {
        RR      <- pred$matRRfit[lag + 1]
        ci_low  <- pred$matRRlow[lag + 1]
        ci_high <- pred$matRRhigh[lag + 1]
        SE      <- (log(ci_high) - log(ci_low)) / (2 * 1.96)
        Z       <- log(RR) / SE
        pvalue  <- 2 * (1 - pnorm(abs(Z)))
        
        results_lag <- rbind(results_lag,
                             data.frame(
                               input     = paste0(pol_1, ifelse(pred_name == "main", "x", "+"), pol_2),
                               pollutant = pol_1,
                               lag       = lag,
                               RR        = RR,
                               ci_low    = ci_low,
                               ci_high   = ci_high,
                               pvalue    = pvalue
                             )) 
      }
      # Store in results_list
      results_list[[pol_1]][[paste0(pol_1, ifelse(pred_name == "main", "x", "+"), pol_2)]]$results_lag <- results_lag
      all_results_lag <- rbind(all_results_lag, results_lag)
    }
  }
}

# Save results
saveRDS(results_list, file = "results_list_bipol.rds")
write.csv(all_results_lag, file = "RR_results_bipol.csv", row.names = FALSE)

