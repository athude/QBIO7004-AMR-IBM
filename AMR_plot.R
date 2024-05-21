# Import libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Import data
spatial_simulation_data <- read.csv("simulation_output/temp_AMR_r4_15.csv")
con_factor = 60/20
p <- list(antibiotic_MIC = 4)
t_max <- max(simple_history$timestep)
run <- 4
i <- 15

# Define ranges of MIC to store results in
MIC_breaks = c(-1, 0, p$antibiotic_MIC, p$antibiotic_MIC*2, p$antibiotic_MIC*4, p$antibiotic_MIC*8, p$antibiotic_MIC*16, Inf)

# Filter into groups and get MIC averages
simplified_data <- spatial_simulation_data %>% 
  mutate(MIC_count = cut(spatial_simulation_data$MIC, breaks = MIC_breaks, right = FALSE)) %>% 
  group_by(timestep, MIC_count, antibiotic) %>% 
  count()

# Create a tibble with all MIC ranges
unique_timesteps <- length(unique(simplified_data$timestep))
all_ranges <- tibble(MIC_count = rep(cut(MIC_breaks, breaks = MIC_breaks, right = FALSE), unique_timesteps),
                     timestep = rep(0:(unique_timesteps-1), each = length(MIC_breaks)))

simple_history <- left_join(all_ranges, simplified_data, by = c("MIC_count", "timestep"), relationship = "many-to-many")
simple_history <- subset(simple_history, !is.na(MIC_count))
simple_history <- simple_history %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>% 
  fill(antibiotic)

# Create labels based on the MIC count brackets
labels <- unique(simple_history["MIC_count"])[3:7,]

# Plot number of bacteria against time
plot <- ggplot(data = simple_history, aes(x = timestep, y = n, color = MIC_count)) +
  geom_col(aes(y = antibiotic*20), color = "gray93", fill = "gray93") +
  geom_line(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, max(simple_history$timestep))) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(0, t_max * con_factor, 2*con_factor),
                     labels = seq(0, t_max, 2)) +
  scale_y_continuous(sec.axis = sec_axis(~.*0.015, name="Antibiotics concentration (µg/mL)")) +
  labs(title = "Evolution of antimicrobial resistance over time",
       x = "Hours",
       y = "Number of bacteria",
       color = "MIC (µg/mL)") +
  scale_color_manual(values = c("black", "chartreuse3", brewer.pal(6, name = "Reds")[2:6]),
                     labels = c("Dead", "Sensitive", labels)) +
  theme_light() +
  theme(axis.ticks.y.right = element_line(color = "darkgrey"),
        axis.text.y.right = element_text(color = "darkgrey"), 
        axis.title.y.right = element_text(color = "darkgrey")
  )

# Save the plot as .jpg
ggsave(file.path("simulation_output/good_results", paste0("temp_AMR_r", run, "_", i, ".jpg")), plot = plot, width = 15, height = 10, units = "cm")
