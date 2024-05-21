library(gganimate)
library(ggplot2)
library(dplyr)

# Import simulation
results_import <- read.csv("simulation_output/temp_AMR_r4_2.csv")

results_import <- results_import %>% 
  filter(!is.na(ID))

animate_bacteria <- function(results_import) {
  ggplot(data = results_import) +
    geom_point(data = results_import, 
               aes(x = x, y = y, color = MIC, group = ID), 
               position = position_jitter(width = 0.2, height = 0.2)) +
    scale_color_distiller(palette = "Reds", direction = 1) +
    labs(title = "Evolution of antimicrobial resistance over time",
         caption = 'Time passed: {round(frame_time / 3, 0)} hrs \n Antibiotics concentration: {results_import$antibiotic}',
         x = NULL,
         y = NULL,
         color = "Bacterial MIC") +
    theme(line = element_blank(),
          panel.background = element_rect(fill = '#f1fafb'),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          plot.caption = element_text(vjust = -0.1, size = 15),
          legend.title = element_text(size = 15)) +
    transition_time(timestep)
}

animation <- animate(animate_bacteria(results_import), height = 600, width = 600)
anim_save(paste0("temp_AMR_evol_GIF_r", run, "_", batch, ".gif"), animation = animation)
