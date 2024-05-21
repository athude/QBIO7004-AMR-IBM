## QBIO7004: Antimicrobial resistance evolution with an individual-based model
*by Alberte Thude (s48598228) for QBIO7004, 2024.*


This repository contains the following files and folders:
- `AMR_simulation_main.R`: the main simulation function, including plot, in which antibiotics concentration is changed over time
- `AMR_simulation_2D.R`: the secondary simulation function in which antibiotics concentration is defined across space
- `AMR_plot.R`: a script for plotting data from the main simulation
- `AMR_animation.R`: a script for animating data from the main simulation
- `AMR_animation_2D.R`: a script for animating data from the secondary simulation
- `/batch_scripts`: a folder with the batch scripts used to run simulations on the HPC for different parameter settings
- `/parameter_scripts`: a folder with R scripts containing the initial values and parameter settings used for different simulation runs
- `/simulation_output`: a folder with .csv, .jpg, and .gif output from simulations

In order to run the simulations, this folder structure should be conserved. Run numbers indicate that changes have been made to initial conditions and/or parameters.
