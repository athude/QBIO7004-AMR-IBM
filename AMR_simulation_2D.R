### QBIO7004 project 2024: Simulation of AMR evolution in bacteria
### Model 2: Spatial changes in antibiotics concentration
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
source(paste0("settings_2D_r", run, ".R"))

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
  x = c(sample(1:p$cells_segment, ini$n_bac, replace = TRUE), rep(NA, ini$max_bac - ini$n_bac)),
  y = c(sample(1:p$cells_segment, ini$n_bac, replace = TRUE), rep(NA, ini$max_bac - ini$n_bac)),
  MIC = c(MIC_list,
          rep(NA, ini$max_bac - ini$n_bac)),
  resistance_status = ifelse(is.na(MIC), NA, ifelse(MIC < p$antibiotic_MIC, "sensitive", "resistant")),
  antibiotic_conc = rep(NA, ini$max_bac),
  resistant_nbours = rep(NA, ini$max_bac),
  total_nbours = rep(NA, ini$max_bac),
  group_survival = rep(NA, ini$max_bac),
  MIC_count = rep(NA, ini$max_bac),
  mut_outcome = rep(NA, ini$max_bac)
) 

################################################################################
### Part functions ----

check_survival <- function(bac_population_alive) {
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

simulate_reproduction <- function(cell, bac_population, empty_slot) {
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
    bac_population$mut_outcome[empty_slot] <- ifelse(runif(1) * (1 - ifelse(cell$MIC > 0, cell$antibiotic_conc / 2 * cell$MIC, 0)) < p$mut_rate, runif(1), 0)
    
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

simulate_spatial_AMR <- function(bac_population, t_max, ini, p) {
  # Create vector for storage of the evolution history
  history <- vector("list", length = t_max)
  
  # Store the initial bacterial population and antibiotics concentration in the detailed history
  history[[1]] <- bac_population %>%
    filter(!is.na(x) & MIC != -1) %>% 
    group_by(x, y) %>% 
    summarise(n = n(),
              avg_MIC = mean(MIC)) %>% 
    mutate(timestep = 0)
  
  for (t in 1:t_max) {
    # Find indices of all living bacteria
    alive <- which(bac_population$resistance_status != "dead")
    
    # Check antibiotic concentration for each bacterium based on location
    for (cell in alive) {
      if (bac_population[[cell, "x"]] <= p$cells_segment) {
        bac_population[[cell, "antibiotic_conc"]] <- as.integer(p$antibiotic_list[1])
      } else if (bac_population[[cell, "x"]] <= 2 * p$cells_segment) {
        bac_population[[cell, "antibiotic_conc"]] <- as.integer(p$antibiotic_list[2])
      } else if (bac_population[[cell, "x"]] <= 3 * p$cells_segment){
        bac_population[[cell, "antibiotic_conc"]] <- as.integer(p$antibiotic_list[3])
      } else {
        bac_population[[cell, "antibiotic_conc"]] <- as.integer(p$antibiotic_list[4])
      }
    }
    
    # Check survival
    bac_population[alive,] <- check_survival(bac_population[alive,])
    
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
      bac_population[empty_slot,] <- simulate_reproduction(bac_population[cell,], bac_population, empty_slot)
    }
    
    # Store the history for every time step
    history[[t+1]] <- bac_population %>%
      filter(!is.na(x) & MIC != -1) %>% 
      group_by(x, y) %>% 
      summarise(n = n(),
                avg_MIC = mean(MIC)) %>% 
      mutate(timestep = t)
  }
  
  history <- bind_rows(history)
  return(history)
}

################################################################################
### Run spatial simulation  ----

# Run simulation
spatial_simulation_result <- simulate_spatial_AMR(bac_population, t_max, ini, p)

# Save as .csv
write.csv(spatial_simulation_result, paste0("spat_AMR_r", run, "_", batch, ".csv"))

################################################################################
### Plot and animate bacteria in the 2D grid ----

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
