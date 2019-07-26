#!/bin/bash

# Name of your job 
#SBATCH --job-name=3Dcoarse
#SBATCH --partition=compute

# Specify the name of the output file. The %j specifies the job ID
#SBATCH --output=3Dcoarse.o%j

# Specify a name of the error file. The %j specifies the job ID
#SBATCH --error=3Dcoarse.e%j

# The walltime you require for your simulation
#SBATCH --time=6:00:00

# Number of nodes you are requesting for your job.
#SBATCH --nodes=1

# Number of processors per node
#SBATCH --ntasks-per-node=24

# Send an email to this address when you job starts and finishes
# CHANGE THIS EMAIL ADDRESS
#SBATCH --mail-user=name@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


# Start in submit directory of main job

# Change to local scratch
export curr_dir_local=$(pwd)
export up_dir_local=$(dirname "$PWD")
export dir=$(basename "${up_dir_local}")
export SCRATCH_DIR_LOCAL=/scratch/${USER}/${SLURM_JOB_ID}/${dir}
mkdir -p ${SCRATCH_DIR_LOCAL} 

cp -pr * ${SCRATCH_DIR_LOCAL}
cd ${SCRATCH_DIR_LOCAL}
mkdir Outputs/


./svpre MODELNAME.svpre
ibrun ./svsolver
python Postsolve_master.py
wait

# Put output file where the main dakota job looks for it
cp -pr ${SCRATCH_DIR_LOCAL}/Outputs/* ${up_dir_local}/Outputs/

# Remove the directory we made
cd ${up_dir_local}
rm -rf ${SCRATCH_DIR_LOCAL}
