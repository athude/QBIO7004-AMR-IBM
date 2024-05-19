### QBIO7004 project 2024: Simulation of AMR evolution in bacteria
### Model 1: Temporal changes in antibiotics concentration
### By Alberte Thude (s4859822)

################################################################################
### Initial document setup ----

# Load libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gganimate)

# Retrieve arguments from batch input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Command input missing :(")
} else {
  batch <- as.integer(args[1])
  run <- as.integer(args[2])
}

################################################################################
### Initial values and parameters ----

# Conversion factor from time step in min to hours (given by 60/doubling time)
con_factor <- 60/20

# Setting number of time steps to simulate
t_max <- con_factor*1  # number of time steps to simulate

# Loading initial values and parameters from helper script
source(paste0("settings_r", run, ".R"))

# Initializing population
MIC_list <- numeric(ini$n_bac)

for (bac in 1:ini$n_bac) {
  MIC_list[[bac]] <-   if (runif(1) < ini$resistance_perc) {
    round(runif(1, 4, 7.9), 2)
  } else {
    round(runif(1, 0, 3.9), 2)} 
}

bac_population <- tibble(
  ID = c(1:ini$n_bac, rep(NA, ini$max_bac - ini$n_bac)),
  x = c(sample(1:ini$n_grid, ini$n_bac, replace = TRUE), rep(NA, ini$max_bac - ini$n_bac)),
  y = c(sample(1:ini$n_grid, ini$n_bac, replace = TRUE), rep(NA, ini$max_bac - ini$n_bac)),
  MIC = c(MIC_list,
          rep(NA, ini$max_bac - ini$n_bac)),
  resistance_status = ifelse(is.na(MIC), NA, ifelse(MIC < p$antibiotic_MIC, "sensitive", "resistant")),
  resistant_nbours = rep(NA, ini$max_bac),
  total_nbours = rep(NA, ini$max_bac),
  group_survival = rep(NA, ini$max_bac),
  MIC_count = rep(NA, ini$max_bac),
  mut_outcome = rep(NA, ini$max_bac)
) 


################################################################################
### Part functions ----

check_survival <- function(bac_population_alive, antibiotic_conc) {
  bac_population_alive %>% 
    group_by(x, y) %>% 
    mutate(total_nbours = as.integer(sum(!resistance_status == "dead"))) %>% 
    mutate(resistant_nbours = as.integer(sum(resistance_status == "resistant"))) %>%
    mutate(group_survival = ifelse(any(MIC > p$antibiotic_MIC & sum(MIC > antibiotic_conc)) > 0 &
                                     antibiotic_conc < 2*p$antibiotic_MIC,
                                   ifelse((runif(1) < p$protective_factor * resistant_nbours / total_nbours), 
                                          TRUE, FALSE),
                                   FALSE)) %>%
    ungroup() %>%
    mutate(resistance_status = ifelse(MIC < antibiotic_conc & group_survival == FALSE, "dead", resistance_status),
           MIC = ifelse(resistance_status == "dead", -1, MIC)) # categorize dead cells as -1 for MIC count in history

  return(bac_population_alive)
}

simulate_HGT <- function(cell) {
  # Generate probability of meeting a resistant neighbour bacterium
  cell_meeting <- as.numeric(ifelse(cell["resistant_nbours"] > 0, 
                                    max(runif(as.integer(cell["total_nbours"])))*(cell["resistant_nbours"]/cell["total_nbours"]), 
                                    0))
  
  # HGT happens if the probability of meeting a resistant bacterium times the conjugation probability is >0.5. 
  # HGT happens at a lower rate for more dense cells in the grid.
  # HGT is assumed to confer a resistance to the minimum level determined as resistant.
  HGT_outcome <- ifelse(cell_meeting * p$conjugation_prob > 0.5 &
                          cell["MIC"] < p$antibiotic_MIC & 
                          runif(1) > cell["total_nbours"] / ini$max_bac_cell, 
                        TRUE, FALSE)
  
  cell["MIC"] <- ifelse(HGT_outcome, p$antibiotic_MIC, cell["MIC"])
  cell["resistance_status"] <- ifelse(HGT_outcome, "resistant", cell["resistance_status"])
  
  return(cell)
}

simulate_movement <- function(cell, bac_population) {
  # A bacterium can attempt a random "jump" in +- 1 cell in either direction, including diagonals.
  x_move <- sample((max(c(0, cell$x - 1))):(min(c(cell$x + 1, p$ini$n_grid))), 1)
  y_move <- sample((max(c(0, cell$y - 1))):(min(c(cell$y + 1, p$ini$n_grid))), 1)
  
  # The number of cells already present in the attempted cell is calculated.
  n_cells_move <- bac_population %>% filter(bac_population$resistance_status != "dead" & bac_population$x == x_move & bac_population$y == y_move) %>% count()
  
  n_cells_adjoining <- bac_population %>% 
    filter(bac_population$resistance_status != "dead" & !is.na(bac_population$ID)) %>% 
    group_by(x, y) %>% 
    filter(x == x_move + 1 & y == y_move | 
             x == x_move + 1 & y == y_move + 1| 
             x == x_move + 1 & y == y_move - 1| 
             x == x_move - 1 & y == y_move |
             x == x_move - 1 & y == y_move + 1|
             x == x_move - 1 & y == y_move - 1|
             x == x_move & y == y_move + 1 |
             x == x_move & y == y_move - 1 |
             x == x_move & y == y_move) %>% 
    count()
  
  
  # If there is space in the new cell, the bacterium moves to it.
  movement_outcome <- ifelse(n_cells_move < ini$max_bac_cell &
                               sum(n_cells_adjoining[3]) > 1, TRUE, FALSE)

  cell["x"] <- ifelse(movement_outcome[1], as_tibble(x_move), cell["x"])
  cell["y"] <- ifelse(movement_outcome[1], as_tibble(y_move), cell["y"])
  
  return(cell)
}
  
simulate_reproduction <- function(cell, bac_population, empty_slot, antibiotic_conc) {
  # Check how many neighbours the bacterium currently has
  group <- which(bac_population$resistance_status != "dead" & bac_population$x == cell$x & bac_population$y == cell$y)
  group_nbours <- length(group)
  
  # If there is space within the cell, the bacterium may reproduce in a density-dependent manner
  if (empty_slot < ini$max_bac & runif(1) > group_nbours/ini$max_bac_cell) {
    # Give daughter cell an ID and same position as mother cell
    bac_population$ID[empty_slot] <- empty_slot
    bac_population$x[empty_slot] <- cell$x
    bac_population$y[empty_slot] <- cell$y
    
    # Determine whether daughter cell mutates and generate a random number to determine mutation type
    bac_population$mut_outcome[empty_slot] <- ifelse(runif(1) * (1 - ifelse(cell$MIC > 0, antibiotic_conc / 2 * cell$MIC, 0)) < p$mut_rate, runif(1), 0)
    
    # Generate MIC for daughter cell based on mutation outcome
    if (bac_population$mut_outcome[empty_slot] > 0) {
      if (bac_population$mut_outcome[empty_slot] < p$mut_ben_prob) {
        # Increase MIC for beneficial mutations
        bac_population$MIC[empty_slot] <- cell$MIC * (1 + abs(rnorm(1, 0, 0.05)))
        
      } else if (bac_population$mut_outcome[empty_slot] < p$mut_lethal_prob) {
        # Lethal mutations reduce MIC to -1, which is a proxy for dead bacteria
        bac_population$MIC[empty_slot] <- -1
        
      } else {
        # Decrease MIC for deleterious mutations
        bac_population$MIC[empty_slot] <- max(0, cell$MIC * (1 - abs(rnorm(1, 0, 0.2))))
      }
    } else {
      # If no mutations occur, the MIC remains the same as for the mother cell
      bac_population$MIC[empty_slot] <- cell$MIC
    }
    
    # Add resistance status for the daughter cell based on the MIC
    bac_population$resistance_status[empty_slot] <- ifelse(bac_population$MIC[empty_slot] == -1, "dead",
                                                           ifelse(bac_population$MIC[empty_slot] < 4, "sensitive", "resistant"))
  }
  
  return(bac_population[empty_slot,])
}

################################################################################
### Temporal simulation function ----

simulate_temp_AMR <- function(bac_population, t_max, ini, p) {
  # Create vectors for storage of the evolution history and antibiotics concentration
  simple_history <- vector("list", length = t_max)
  spatial_history <- vector("list", length = t_max)
  antibiotics_concentration <- vector("list", length = t_max)
  
  # Define ranges of MIC to store results in
  MIC_breaks = c(-1, 0, p$antibiotic_MIC, p$antibiotic_MIC*2, p$antibiotic_MIC*4, p$antibiotic_MIC*8, p$antibiotic_MIC*16, Inf)
  
  # Store the counts within MIC ranges for the population
  bac_population$MIC_count <- cut(bac_population$MIC, breaks = MIC_breaks, right = FALSE)
  MIC_history <- bac_population %>% count(MIC_count)
  
  # Create a tibble representing all ranges
  all_ranges <- tibble(MIC_count = cut(MIC_breaks, breaks = MIC_breaks, right = FALSE),
                       timestep = 0)
  
  # Load initial antibiotic concentration
  antibiotic_conc <- ini$antibiotic_conc
  
  # Store the initial MIC counts and antibiotics concentration in the simple history
  simple_history[[1]] <- left_join(all_ranges, MIC_history, by = "MIC_count") %>% 
    mutate(n = ifelse(is.na(n), 0, n))
  simple_history[[1]]$antibiotic <- antibiotic_conc
  
  # Store the initial bacterial population and antibiotics concentration in the detailed history
  spatial_history[[1]] <- bac_population[2:4] %>% 
    mutate(timestep = 0)
  spatial_history[[1]]$antibiotic <- antibiotic_conc
  
  # Create a vector with the time steps of antibiotics administration
  antibiotic_increase <- p$dose_first + seq(0, p$dose_number-1, 1) * p$dose_interval
  
  for (t in 1:t_max) {
    # Adjust antibiotic concentration
    antibiotic_conc <- ifelse(t %in% antibiotic_increase, 
                             antibiotic_conc * exp(-log(2)/p$half_life) + p$dose_conc,
                             antibiotic_conc * exp(-log(2)/p$half_life))
    
    # Find indices of all living bacteria
    alive <- which(bac_population$resistance_status != "dead")
    
    # Check survival
    bac_population[alive,] <- check_survival(bac_population[alive,], antibiotic_conc)
      
    # Update indices of surviving bacteria
    alive <- which(bac_population$resistance_status != "dead")
    
    # Update parameters for surviving bacteria
    bac_population[alive,] <- bac_population[alive,] %>% 
      group_by(x, y) %>% 
      mutate(total_nbours = as.integer(sum(!resistance_status == "dead"))) %>% 
      mutate(resistant_nbours = as.integer(sum(resistance_status == "resistant")))
    
    # Simulate HGT, movement, and reproduction for each alive cell
    for (cell in alive) {
      bac_population[cell,] <- simulate_HGT(bac_population[cell,])
      bac_population[cell,] <- simulate_movement(bac_population[cell,], bac_population)
      
      # Find next empty slot in the population tibble for reproduction
      empty_slot <- which(is.na(bac_population$ID))[1]
      bac_population[empty_slot,] <- simulate_reproduction(bac_population[cell,], bac_population, empty_slot, antibiotic_conc)
    }
    
    # Store status after every time step
    bac_population$MIC_count <- cut(bac_population$MIC, breaks = MIC_breaks, right = FALSE)
    MIC_history <- bac_population %>% 
      filter(!is.na(ID)) %>% 
      count(MIC_count)
    
    # Create a tibble representing all ranges at every time step
    all_ranges <- tibble(MIC_count = cut(MIC_breaks, breaks = MIC_breaks, right = FALSE),
                         timestep = t)
    
    # Store the MIC counts and antibiotics concentration in the simple history for every time step
    simple_history[[t+1]] <- left_join(all_ranges, MIC_history, by = "MIC_count") %>% 
      mutate(n = ifelse(is.na(n), 0, n))
    simple_history[[t+1]]$antibiotic <- antibiotic_conc
    
    # Store the bacterial population and antibiotics concentration in the detailed history for every time step
    spatial_history[[t+1]] <- bac_population[2:4] %>% 
      mutate(timestep = t)
    spatial_history[[t+1]]$antibiotic <- antibiotic_conc
    
  }
  
  simple_history <- bind_rows(simple_history)
  spatial_history <- bind_rows(spatial_history)
  return(c(simple_history, spatial_history))
}

################################################################################
### Run temporal simulation  ----

# Run simulation
simulation_result <- simulate_temp_AMR(bac_population, t_max, ini, p)

# Divide results into MIC count and spatial data
simple_simulation_data <- subset(as_tibble(simulation_result[1:4]), !is.na(MIC_count))
spatial_simulation_data <- subset(as_tibble(simulation_result[5:9]), !is.na(x) & MIC != -1)

# Save as .csv
write.csv(simulation_result, paste0("temp_AMR_r", run, "_", batch, ".csv"))

################################################################################
### Plot number of bacteria against time  ----

# Create labels based on the MIC count brackets
labels <- unique(simple_simulation_data["MIC_count"])[3:7,]

# Plot number of bacteria against time
plot <- ggplot(data = simple_simulation_data, aes(x = timestep, y = n, color = MIC_count)) +
  geom_col(aes(y = antibiotic*20), color = "gray93", fill = "gray93") +
  geom_line(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, max(simple_simulation_data$timestep))) +
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
ggsave(paste0("temp_AMR_r", run, "_", batch, ".jpg"), plot = plot, width = 15, height = 10, units = "cm")


################################################################################
### Plot and animate bacteria in the 2D grid ----

animate_bacteria <- function(spatial_simulation_data) {
  ggplot(data = spatial_simulation_data, aes(x = x, y = y, color = MIC)) +
    geom_point() +
    scale_color_distiller(palette = "Reds", direction = 1) +
    labs(title = "Evolution of antimicrobial resistance over time",
         caption = 'Time passed: {round(frame_time / 3, 0)} hrs \n Antibiotics concentration: {spatial_simulation_data$antibiotic}',
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

animation <- animate(animate_bacteria(spatial_simulation_data), height = 600, width = 600)
anim_save(paste0("temp_AMR_evol_GIF_r", run, "_", batch, ".gif"), animation = animation)
