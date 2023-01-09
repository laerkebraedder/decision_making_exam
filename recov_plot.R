# More informative plotting - code-courtesy of Lasse - fanx!
recov_plot <- function(true, infer, plot_lab, plot_col) {
  
  # library(ggplot2)
  
  df <- data.frame(true, infer)
  
  pl <- ggplot(df, aes(x = true,
                       y = infer,
                       color = plot_col)) + #Setting aesthetics for plot
    geom_point() + #Giving points a color each
    geom_smooth(method = "lm", se = T, formula = "y ~ x") +
    theme_minimal() + #Setting theme
    xlab(plot_lab[1]) + #Setting x label
    ylab(plot_lab[2]) + #Setting y label
    labs(color = "") + #Setting legend title
    ggtitle(paste0("'", plot_lab[2], "' compared to '", plot_lab[1], "'"))
  
  return(pl)
  
}