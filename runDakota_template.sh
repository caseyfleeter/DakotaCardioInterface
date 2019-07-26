#!/bin/bash
#
# Name of your job 
#SBATCH --job-name=dakotaMLMF
#SBATCH --partition=compute
#
# Specify the name of the output file. The %j specifies the job ID
#SBATCH --output=dakota.o%j
#
# Specify a name of the error file. The %j specifies the job ID
#SBATCH --error=dakota.e%j
#
# The walltime you require for your simulation
#SBATCH --time=48:00:00
#
# Job priority.
#SBATCH --qos=normal
#
# Number of nodes you are requesting for your job.
#SBATCH --nodes=1
#
# Number of processors per node
#SBATCH --ntasks-per-node=1
#
# Send an email to this address when you job starts and finishes
#SBATCH --mail-user=name@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


# GET ALL SUBMISSION DIRECTORIES
export SUBMIT_DIR=${SLURM_SUBMIT_DIR}
export SCRATCH_DIR=/scratch/${USER}/${SLURM_JOB_ID}
mkdir -p ${SCRATCH_DIR}
echo "Created scratch directory successfully!"

cp dakota_MLMF_template.in dakota_MLMF.in
export WDLR=${SLURM_SUBMIT_DIR}/LF
sed -i -e "s@WORKDIR_LF_REPLACE@'${WDLR}'@g" dakota_MLMF.in
export CLR=${SLURM_SUBMIT_DIR}/Run_0D_Model/*
# export CLR=${SLURM_SUBMIT_DIR}/Run_1D_Model/*
sed -i -e "s@COPY_LF_REPLACE@'${CLR}'@g" dakota_MLMF.in
export WDHR=${SLURM_SUBMIT_DIR}/HF
sed -i -e "s@WORKDIR_HF_REPLACE@'${WDHR}'@g" dakota_MLMF.in
export CHR=${SLURM_SUBMIT_DIR}/Run_3D_Model/*
sed -i -e "s@COPY_HF_REPLACE@'${CHR}'@g" dakota_MLMF.in
echo "Replacement for dakota.in file complete!"


# RUN DAKOTA
./dakota -i dakota_MLMF.in -o dakota_MLMF.out -e dakota_MLMF.err -w dakota_MLMF.rst
# ./dakota -i dakota_MLMF.in -r dakota_MLMF.rst -o dakota_MLMF.out -e dakota_MLMF.err -w dakota_MLMF.rst
wait

# CLEAN UP FILES
rm -rf ${SCRATCH_DIR}

exit 0
