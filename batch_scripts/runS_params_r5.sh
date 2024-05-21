#!/bin/bash
#SBATCH --job-name=AT_mainP5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=qbiol
#SBATCH --mem=2GB              # memory (MB)
#SBATCH --time=0-6:00          # time (D-HH:MM)
#SBATCH -o mainP5.%N.%j.out     # STDOUT
#SBATCH -e mainP5.%N.%j.err     # STDERR
#
#SBATCH --array=1


echo "Start time: "; date

run="5"

module load applications/R/4.1.3

Rscript --vanilla "AMR_simulation_main_params2.R" $SLURM_ARRAY_TASK_ID $run

echo "End time: "; date
