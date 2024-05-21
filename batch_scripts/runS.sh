#!/bin/bash
#SBATCH --job-name=AT_main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=qbiol
#SBATCH --mem=2GB              # memory (MB)
#SBATCH --time=0-2:00          # time (D-HH:MM)
#SBATCH -o main.%N.%j.out     # STDOUT
#SBATCH -e main.%N.%j.err     # STDERR
#
#SBATCH --array=1


echo "Start time: "; date

module load applications/R/4.1.3

for r in 1 2 3
do 
   run="$r"

Rscript --vanilla "AMR_simulation_main.R" $SLURM_ARRAY_TASK_ID $run

done

echo "End time: "; date
