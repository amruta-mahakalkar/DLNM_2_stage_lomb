# BI-POLLUTANT META-ANALYSIS (PM2.5)
# Define the pollutants
pollutants <- c("PM25", "NO2", "O3", "SO2", "CO")

for (pollutant in pollutants) {
  # Extract the results for the current pollutant and its combinations
  results <- results_list$all_CA[[pollutant]]$results_lag
  
  # Initialize an empty list to store results
  combined_results <- list()
  
  # Make paired pollutant names while excluding the pollutant in loop
  for (pol in pollutants) {
    if (pollutant != pol) {
      comb_results <- results_list$all_CA[[paste0(pollutant, "_", pol)]]$results_lag
      comb_results$pollutant <- paste0(pollutant, " + ", pol)
      combined_results <- append(combined_results, list(comb_results))
    }
  }
  
  # Add a column for each dataset indicating the pollutant
  results$pollutant <- pollutant
  combined_results <- append(combined_results, list(results))
  
  # Combine the results
  result_bipol <- bind_rows(combined_results)
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