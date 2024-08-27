plot_exposure_response <- function(results_list, pollutant, groups, colors, ylim) {
  # Plot the overall predictions
  plot(results_list[[groups[1]]][[pollutant]]$meta_pred_tot, 'overall', lwd = 2, 
       xlab = pollutant, ylab = 'RR', ylim = ylim, col = colors[1], 
       ci.arg = list(density = 20, col = grey(0.2)))
  
  # Add the group predictions
  if (length(groups) > 1) {
    for (i in 2:length(groups)) {
      lines(results_list[[groups[i]]][[pollutant]]$meta_pred_tot, "overall", ci = "area", 
            col = colors[i], lwd = 2, ci.arg = list(col = alpha(colors[i], 0.2)))
    }
  }
  
  # Add a legend
  legend("topright", legend = groups, col = colors, lwd = 2, bty = "n", cex = 0.7)
  
  # Add grid
  grid()
}