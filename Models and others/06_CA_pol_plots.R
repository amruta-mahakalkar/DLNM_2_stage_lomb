# CORRELATION ANALYSIS
# Calculate the correlation table
correlation_table <- cor(data[, .(PM25, PM10, NO2, O3, SO2, CO, Temp, RH, all_CA)], method = "spearman")

# Print the correlation table
corrplot(correlation_table, type = 'lower', method = 'square', tl.col = "black", tl.srt = 90, diag = FALSE, tl.cex=0.7, cl.cex = 0.7, col = rev(corrplot::COL2('RdBu')))

# EXPOSURE-LAG-RESPONSE GRAPHS
# MERGE ALL RESULTS
#load results_list
results_base <- readRDS("2_step_results_100k.rds")
results_list_ur <- readRDS("2_step_urb_rur_100k.rds")
results_biPM25 <- readRDS("2_step_bipol_PM25_100k.rds")
results_biO3 <- readRDS("2_step_bipol_O3_100k.rds")
results_biNO2 <- readRDS("2_step_bipol_NO2_100k.rds")
results_biSO2 <- readRDS("2_step_bipol_SO2_100k.rds")
results_biCO <- readRDS("2_step_bipol_CO_100k.rds")
results_season <- readRDS("2_step_season_100k.rds")

# Combine the all_CA sub-groups
combined_all_CA <- c(results_base$all_CA, results_biPM25$all_CA, results_biNO2$all_CA,results_biO3$all_CA, results_biSO2$all_CA,results_biCO$all_CA, results_season$all_CA)

# Combine the lists
results_list <- c(results_base, results_list_ur, results_biPM25, results_biNO2, results_biO3, results_biSO2, results_biCO, results_season)
# Remove all existing all_CA entries
results_list <- results_list[!names(results_list) %in% "all_CA"]
results_list$all_CA <- combined_all_CA

# MELT ALL RESULTS
extract_results_lag <- function(lst) {
  if (is.list(lst)) {
    if ("results_lag" %in% names(lst)) {
      return(lst$results_lag)
    } else {
      return(rbindlist(lapply(lst, extract_results_lag), fill = TRUE))
    }
  } 
}
results_lag_df <- extract_results_lag(results_list)

# Save all_results_lag to a CSV file
write.csv(results_lag_df, file = "2_step_all_results_100k.csv", row.names = FALSE)


# PLOT EXP-LAG GRAPHS 
for (pollutant in pollutants) {
  file_name <- paste0("graphs/", pollutant, "_age_sex.jpeg")
  plot_exposure_lag(
    results_list, 
    c("all_CA", "male", "female", "adult", "senior"), 
    pollutant, 
    c("all_CA" = "black", "male" = "deepskyblue", "female" = "purple3", "adult" = "darkolivegreen3", "senior" = "orange"), 
    "age-sex", 
    file_name
  )
}

for (pollutant in pollutants) {
  file_name <- paste0("graphs/", pollutant, "_urban_rural.jpeg")
  plot_exposure_lag(
    results_list, 
    c("all_CA", "urban", "rural"), 
    pollutant, 
    c("all_CA" = "black", "urban" = "salmon2", "rural" = "darkolivegreen3"), 
    "location", 
    file_name)
}

for (pollutant in pollutants) {
  file_name <- paste0("graphs/", pollutant, "_season.jpeg")
  plot_exposure_lag(
    results_list, 
    c("all_CA", "Summer", "Autumn", "Winter", "Spring"), 
    pollutant, 
    c("all_CA" = "black", "Summer" = "salmon2", "Autumn" = "brown", "Winter" = "skyblue2", "Spring" = "violet"), 
    "season", 
    file_name)
}

# PLOT EXP-RES GRAPHS
# plot urban-rural exp-res curve
pollutants <- c("PM25", "PM10", "NO2", "O3", "SO2", "CO")
colors <- c("black","tomato3", "green4")
groups <- c("all_CA", "urban", "rural")
ylim <- list("PM25"= c(0.995, 1.006), "PM10"= c(0.995, 1.006), "O3"= c(0.995, 1.005), "NO2"= c(0.993, 1.009), "SO2"= c(0.92, 1.07), "CO"= c(0.9994, 1.0006))  
for (pollutant in pollutants) {
  # Construct the file name dynamically
  file_name <- paste0("graphs/exposure_response/", pollutant, "_ur.jpeg")
  
  # Call the plot function
  plot_exposure_response(results_list, pollutant, groups, colors, ylim[[pollutant]], file_name)
}

# plot age exp-res curve
colors <- c("black","darkred", "darkturquoise")
groups <- c("all_CA", "adult", "senior")
ylim <- list("PM25"= c(0.996, 1.008), "PM10"= c(0.996, 1.008), "O3"= c(0.995, 1.005), "NO2"= c(0.992, 1.009), "SO2"= c(0.94, 1.09), "CO"= c(0.9993, 1.0009))  
for (pollutant in pollutants) {
  # Construct the file name dynamically
  file_name <- paste0("graphs/exposure_response/", pollutant, "_age.jpeg")
  
  # Call the plot function
  plot_exposure_response(results_list, pollutant, groups, colors, ylim[[pollutant]], file_name)
}

# plot sex exp-res curve
colors <- c("black","deepskyblue","salmon2")
groups <- c("all_CA", "male", "female")
ylim <- list("PM25"= c(0.996, 1.008), "PM10"= c(0.996, 1.005), "O3"= c(0.995, 1.005), "NO2"= c(0.992, 1.006), "SO2"= c(0.92, 1.06), "CO"= c(0.9993, 1.0009))  
for (pollutant in pollutants) {
  # Construct the file name dynamically
  file_name <- paste0("graphs/exposure_response/", pollutant, "_sex.jpeg")
  
  # Call the plot function
  plot_exposure_response(results_list, pollutant, groups, colors, ylim[[pollutant]], file_name)
}

# PLOT BI-POLLUTANTS LAG-RESPONSE GRAPHS 
# Define the pollutants
pollutants <- c("PM25", "NO2", "O3", "SO2", "CO")

for (pollutant in pollutants) {
  # Extract the results for the current pollutant and its combinations
  results <- results_list$all_CA[[pollutant]]$results_lag
  results_PM25 <- results_list$all_CA[[paste0(pollutant, "_PM25")]]$results_lag
  results_O3 <- results_list$all_CA[[paste0(pollutant, "_O3")]]$results_lag
  results_NO2 <- results_list$all_CA[[paste0(pollutant, "_NO2")]]$results_lag
  results_SO2 <- results_list$all_CA[[paste0(pollutant, "_SO2")]]$results_lag
  results_CO <- results_list$all_CA[[paste0(pollutant, "_CO")]]$results_lag
  
  # Add a column for each dataset indicating the pollutant
  results$pollutant <- pollutant
  results_PM25$pollutant <- paste0(pollutant, " + PM25")
  results_O3$pollutant <- paste0(pollutant, " + O3")
  results_NO2$pollutant <- paste0(pollutant, " + NO2")
  results_SO2$pollutant <- paste0(pollutant, " + SO2")
  results_CO$pollutant <- paste0(pollutant, " + CO")
  
  # Combine the results
  result_bipol <- bind_rows(results, results_PM25, results_O3, results_NO2, results_SO2, results_CO)
  result_bipol$label <- ifelse(result_bipol$pvalue < 0.05, "*", "")
  # Construct the file name dynamically
  file_name <- paste0("graphs/bipol_", pollutant, ".jpeg")
  # Create the plot
  
  p <- ggplot(result_bipol, aes(x = lag, y = RR, color = pollutant)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.85) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
    labs(title = paste0(pollutant, " bipol"), x = "Lag", y = "RR") +
    scale_x_continuous(breaks = unique(result_bipol$lag)) +
    geom_text(aes(label = label, y = ci_high), vjust=0, position = position_dodge(0.5), show.legend = FALSE) +
    theme_minimal() + 
    theme(panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey90"),  
          panel.grid.minor.y = element_line(color = "grey90")) +
    guides(color = guide_legend(title = NULL))
  ggsave(file_name, plot = p, width = 6.25, height = 3.75, dpi = 400, quality = 100)
}