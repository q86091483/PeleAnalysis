#!/bin/bash
#SBATCH --job-name=nonreact
#SBATCH --account=rrg-bsavard
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=8:00:00
#SBATCH --output=mpi_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

module load gcc openmpi

./FlameStructure_basic.sh > FlameStructure_basic.log 2>&1

echo ">>ending:  " `date`
