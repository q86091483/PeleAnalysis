#!/bin/bash
#SBATCH --job-name=Touch
#SBATCH --account=rrg-bsavard
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=5:00:00
#SBATCH --output=mpi_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

module load gcc openmpi hdf5

./Touch.sh > Touch.log 2>&1

echo ">>ending:  " `date`
