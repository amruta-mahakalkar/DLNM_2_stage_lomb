plot_exposure_lag <- 
  function(results_list, groups, pollutant, colors, group_labels) { 

  group_results <- lapply(groups, function(group) results_list[[pollutant]][[group]]$results_lag)

  # Merge the results
  df_results <- dplyr::bind_rows(group_results)
  
  # Add a new column to support plot into
  df_results$pvalue <- as.numeric(df_results$pvalue)
  df_results$sign_label1 <- ifelse(df_results$pvalue < 0.05, "*", "")
  
  pollutant_label <- switch(pollutant,
                            "PM25" = bquote(PM[2.5]),
                            "PM10" = bquote(PM[10]),
                            "NO2"  = bquote(NO[2]),
                            "O3"   = bquote(O[3]),
                            "SO2"  = bquote(SO[2]),
                            "CO"   = bquote(CO)
  )
  y_label <- if (pollutant == "SO2") {
    bquote("RR per 1" ~ µg/m^3 ~ "increase in" ~ .(pollutant_label))
  } else if (pollutant == "CO") {
    bquote("RR per 100" ~ µg/m^3 ~ "increase in" ~ .(pollutant_label))
  } else {
    bquote("RR per 10" ~ µg/m^3 ~ "increase in" ~ .(pollutant_label))
  }
  
  # reorder plots
  df_results <- df_results %>%
    mutate(group = factor(input, levels = groups)) 
  df_results <- df_results %>%
    mutate(pollutant = factor(pollutant, levels = pollutants))
  df_results$input <- factor(df_results$input, levels = groups)

  
  # Plot the data
  p <- ggplot(df_results, aes(x = lag, y = RR, color = input)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.75) + 
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 1, color = "black", linetype = 1, size=0.3) +
    labs(x = "Lag (in days)", y = y_label, color = "Group", axis.title=element_text(size=4)) +
    annotate("label", x = 6.5, y = Inf, 
            label = deparse(bquote(bold(.(pollutant_label)))),
             hjust = 0, vjust = 1, 
             label.size = 0.1, 
             label.r = unit(0, "lines"), # for sharp box edges
             fill = alpha("white", 0.9),  
             size = 4, # font size                  
             parse = TRUE) +
    scale_x_continuous(breaks = unique(df_results$lag)) +
    scale_color_manual(values = colors, labels = group_labels[groups]) +
    theme_minimal() + 
    theme(panel.grid.major.x = element_line(color = "grey90", size = 0.2),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey90", size = 0.2),  
          panel.grid.minor.y = element_line(color = "grey90", size = 0.2),
          legend.position = "bottom",
          legend.text=element_text(size=12),
          legend.justification = "center",
          axis.title.y = element_text(margin = margin(r = 5))
          ) + 
    guides(color = guide_legend(title = NULL)) +
    geom_text(aes(label = sign_label1, y = ci_high), vjust=0.1, position = position_dodge(0.5), size = 4, show.legend = FALSE)
   
  return(p)
}