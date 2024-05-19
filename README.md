## QBIO7004: Antimicrobial resistance evolution with an individual-based model
*by Alberte Thude (s48598228) for QBIO7004, 2024.*


This repository contains the following files and folders:
- `AMR_simulation_main.R`: the main simulation function, including plot and animation, in which antibiotics concentration is changed over time
- `AMR_simulation_2D.R`: the secondary simulation function, including animation, in which antibiotics concentration is defined across space
- `/batch_scripts`: a folder with the batch scripts used to run simulations on the HPC for different parameter settings
- `/parameter_scripts`: a folder with R scripts containing the initial values and parameter settings used for different simulation runs
- `/simulation_output`: a folder with .csv, .jpg, and .gif output from simulations
- `simulation_run_overview.html`: a HTML file containing settings for all the runs

In order to run the simulations, this folder structure should be conserved. Run numbers indicate that changes have been made to initial conditions and/or parameters.
