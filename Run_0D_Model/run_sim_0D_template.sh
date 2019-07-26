#!/bin/bash

# Name of your job 
#SBATCH --job-name=run0DModel
#SBATCH --partition=compute

# Specify the name of the output file. The %j specifies the job ID
#SBATCH --output=run0DModel.o

# Specify a name of the error file. The %j specifies the job ID
#SBATCH --error=run0DModel.e%j

# The walltime you require for your simulation
#SBATCH --time=0:30:00

# Number of nodes you are requesting for your job.
#SBATCH --nodes=1

# Number of processors per node
#SBATCH --ntasks-per-node=1


# Start in submit directory of main job

# Change to local scratch
export curr_dir_local=$(pwd)
export dir=$(basename "$PWD")
export SCRATCH_DIR_LOCAL=/scratch/${USER}/${SLURM_JOB_ID}/${dir}
mkdir -p ${SCRATCH_DIR_LOCAL} 
cp -pr * ${SCRATCH_DIR_LOCAL}
cd ${SCRATCH_DIR_LOCAL}

python DAKOTA_driver.py -d0D MODELNAME
wait

# Put output file where the main dakota job looks for it
cp -pr ${SCRATCH_DIR_LOCAL}/Outputs/* ${curr_dir_local}/Outputs/

# Remove the directory we made
cd ${curr_dir_local}
rm -rf ${SCRATCH_DIR_LOCAL}
