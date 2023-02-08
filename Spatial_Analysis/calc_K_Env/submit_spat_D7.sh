#!/bin/bash
#SBATCH --account=def-gmitsis 
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26        # number of processes per node
#SBATCH --mem=4G         # memory; default unit is megabytes
#SBATCH --time=2-00:00           # time (DD-HH:MM)
#SBATCH --mail-user=nikolaos.dimitriou@mail.mcgill.ca
#SBATCH --mail-type=ALL         

module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore calcEnv_D7.R
