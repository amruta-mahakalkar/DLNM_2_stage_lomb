# Test RR differences
z_test <- function(data, group_members, ref_group) {
  data_subset <- data %>%
    filter(input %in% group_members) %>%
    mutate(
      log_rr = log(RR),
      se     = (log(ci_high) - log(ci_low)) / (2 * 1.96)
    )
  
  # Reference group
  ref_data <- data_subset %>%
    filter(input == ref_group) %>%
    select(pollutant, lag, ref_log_rr = log_rr, ref_se = se)
  
  # Join and calculate Z-test
  joined <- data_subset %>%
    left_join(ref_data, by = c("pollutant", "lag")) %>%
    filter(!is.na(ref_log_rr)) %>%
    mutate(
      z         = (log_rr - ref_log_rr) / sqrt(se^2 + ref_se^2),
      pvalue_z  = 2 * (1 - pnorm(abs(z))),
      ref_group = ref_group
    )
  
  return(joined)
}

# Run z-test
ages_all <- c("adult", "senior") 
gender_all <- c("female", "male") 
urbanity <- c("urban_all_CA", "rural_all_CA") 
seasons <- c("Summer", "Autumn", "Winter", "Spring") 
lomb_all <- c("all_CA", "lomb_all_CA")
lomb_age <- c("lomb_senior", "lomb_adult")
lomb_sex <- c("lomb_female", "lomb_male")

z_pol<- data.frame()
z_pol <- bind_rows(
  z_test(results_lag_df, ages_all, "adult"), 
  z_test(results_lag_df, gender_all, "male"),
  z_test(results_lag_df, lomb_all, "all_CA"),
  z_test(results_lag_df, lomb_age, "lomb_adult"), 
  z_test(results_lag_df, lomb_sex, "lomb_male"), 
  z_test(results_lag_df, urbanity, "urban_all_CA"),
  z_test(results_lag_df, seasons, "Winter"))

z_pol <- z_pol %>%
  dplyr::distinct(input, pollutant, lag, RR, ci_low, ci_high, .keep_all = TRUE)
write.csv(z_pol, file = "z_tests.csv", row.names = FALSE)