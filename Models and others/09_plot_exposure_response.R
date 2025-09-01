plot_exposure_response <- function(results_list, pollutant, groups, colors, ylim) { 
  pollutant_label <- switch(pollutant,
                            "PM25" = bquote(PM[2.5] ~ " (in μg/m³)"),
                            "PM10" = bquote(PM[10] ~ " (in μg/m³)"),
                            "NO2"  = bquote(NO[2] ~ " (in μg/m³)"),
                            "O3"   = bquote(O[3] ~ " (in μg/m³)"),
                            "SO2"  = bquote(SO[2] ~ " (in μg/m³)"),
                            "CO"   = "CO (in μg/m³)"
  )
  
  par(mar = c(3, 3, 2, 1), mgp = c(2, 0.5, 0))
  
  # Plot overall predictions
  plot(results_list[[pollutant]][[groups[1]]]$meta_pred, 'overall', lwd = 1.5, 
       xlab = pollutant_label, ylab = 'Relative Risk (RR)', ylim = ylim, col = colors[1], cex.lab = 1,
       cex.axis = 1, ci.arg = list(density = 20, col = grey(0.3)))
  
  # Add grid
  abline(h = 1, col = "darkgrey", lty = 1, lwd = 1)
  abline(h = pretty(ylim)[pretty(ylim) != 1], col = "lightgrey", lty = 1, lwd = 0.5)
  abline(v = pretty(par("usr")[1:2]), col = "lightgrey", lty = 1, lwd = 0.5)
  
  # Add the group predictions
  if (length(groups) > 1) {
    for (i in 2:length(groups)) {
      lines(results_list[[pollutant]][[groups[i]]]$meta_pred_tot, "overall", ci = "area", 
            col = colors[i], lwd = 1.5, ci.arg = list(col = alpha(colors[i], 0.2)))
    }
  }
}