plot_exposure_lag <- 
  function(results_list, groups, pollutant, colors, group_name, filename) {
  # Extract the results for each group
  group_results <- lapply(groups, function(group) results_list[[group]][[pollutant]]$results_lag)
  
  # Merge the results
  df_results <- do.call(rbind, group_results)
  
  # Add a new column to distinguish the groups
  df_results$group <- factor(rep(groups, each = nrow(group_results[[1]])), levels = groups)
  
  # Plot the data
  # Open a JPEG device
  jpeg(filename, width = 2500, height = 1500, res = 400, quality = 100)
  plot.new()
  p <- ggplot(df_results, aes(x = lag, y = RR, color = group)) +
    geom_point(position = position_dodge(width = 0.4), size = 1) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
    labs(title = paste("Exposure-lag to", pollutant, "by", group_name), x = "Lag", y = "RR", color = "Group") +
    scale_x_continuous(breaks = unique(df_results$lag)) +
    scale_color_manual(values = colors) +
    theme_minimal() + 
    theme(panel.grid.major.x = element_line(color = "grey90", size = 0.2),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey90", size = 0.2),  
          panel.grid.minor.y = element_line(color = "grey90", size = 0.2))
  
  print(p)
  # Close the JPEG device
  dev.off()
}