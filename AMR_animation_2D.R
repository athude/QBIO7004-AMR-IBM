# Import libraries
library(gganimate)
library(ggplot2)
library(dplyr)

# Import simulation
spatial_simulation_result <- read.csv("simulation_output/spatial_AMR_r4_2.csv")
p$antibiotic_list <- c(0, 4, 8, 32)
con_factor <- 60/20


plot_bacteria <- function(spatial_simulation_result, p) {
    simulation_result$x <- as.integer(spatial_simulation_result$x)
    simulation_result$y <- as.integer(spatial_simulation_result$y)
    ggplot() +
      geom_point(data = spatial_simulation_result, aes(x = x, y = y, size = n, color = avg_MIC, group = timestep), position = position_nudge(x = -0.5, y = -0.5)) +
      scale_x_continuous(limits = c(0, 20), labels = NULL, breaks = c(seq(-0.5, 20.5, 1))) +
      scale_y_continuous(limits = c(0, 10), labels = NULL, breaks = c(seq(0.5, 9.5, 1))) +
      scale_color_distiller(palette = "Reds", direction = 1) +
      geom_vline(xintercept = 0, linetype="dashed", color = "darkgrey") +
      geom_vline(xintercept = 5, linetype="dashed", color = "darkgrey") +
      geom_vline(xintercept = 10, linetype="dashed", color = "darkgrey") +
      geom_vline(xintercept = 15, linetype="dashed", color = "darkgrey") +
      geom_vline(xintercept = 20, linetype="dashed", color = "darkgrey") +
      scale_size(range = c(4,32), guide = 'none') +
      labs(title = "Evolution of antimicrobial resistance over time",
           caption = 'Time passed: {round(frame_time / 3, 0)} hrs',
           x = NULL,
           y = NULL,
           color = "Average MIC") +
      coord_cartesian(clip = "off") +
      geom_text(x=2.5, y=-0.75, size = 7, aes(label=paste(p$antibiotic_list[1], "µg/mL"))) +
      geom_text(x=7.5, y=-0.75, size = 7, aes(label=paste(p$antibiotic_list[2], "µg/mL"))) +
      geom_text(x=12.5, y=-0.75, size = 7, aes(label=paste(p$antibiotic_list[3], "µg/mL"))) +
      geom_text(x=17.5, y=-0.75, size = 7, aes(label=paste(p$antibiotic_list[4], "µg/mL"))) +
      theme(plot.title = element_text(hjust = 0.5, size = 35),
            plot.caption = element_text(vjust = -3, size = 20),
            legend.title = element_text(size = 20),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_line(linewidth = 0.1, linetype = "dashed"),
            panel.grid.minor.y = element_line(linewidth = 0.1, linetype = "dashed"),
            legend.key.size = unit(2, 'cm')) +
      transition_time(timestep)
}


animation <- animate(plot_bacteria(spatial_simulation_result, p), height = 600, width = 1200)
anim_save(paste0("spatial_AMR_evol_GIF_r", run, "_", batch, ".gif"), animation = animation)
