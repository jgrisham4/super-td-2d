#!/bin/bash
#
#SBATCH -N 1                  # Number of nodes requested
#SBATCH -n 1                  # Number of MPI tasks
#SBATCH -p normal             # Queue name
#SBATCH -J triple-deck        # Job name
#SBATCH -o output.%j          # Outlput name
#SBATCH -t 20:00:00           # Runtime (hh:mm:ss)
#SBATCH -A Aero               # Project/allocation name
#SBATCH --mail-user=james.grisham@mavs.uta.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#---- end of options to batch -------
#
# submit with sbatch sbatch.sh
./triple-deck
