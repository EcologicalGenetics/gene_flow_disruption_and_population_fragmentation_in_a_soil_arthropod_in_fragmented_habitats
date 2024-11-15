#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 12G
#SBATCH --cpus-per-task 8
#SBATCH --time 24:00:00

/faststorage/project/EcoGenetics/people/Sarah/fsc28/fsc28_linux64/fsc28 -i *.par -k 1000000 -m -n 10 -s 0 -c 8
