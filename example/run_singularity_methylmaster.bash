#!/usr/bin/env bash

#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=methylmaster_singularity
#SBATCH --output=methylmaster_singularity_out_%j.txt

## Example Slurm script for use with MethylMasteR Singularity container

/usr/bin/singularity run \
-B /path/to/work/dir:/home/docker/work \
path/to/methylmaster_base_slim_script.sif
