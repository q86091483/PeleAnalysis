#!/bin/bash
#SBATCH --job-name=postMMX
#SBATCH --account=rrg-bsavard
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=40
#SBATCH --time=10:00:00
#SBATCH --output=mpi_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

module load gcc openmpi

./MicroMix.sh > MicroMix.log 2>&1

echo ">>ending:  " `date`
