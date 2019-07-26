#!/bin/bash

# EXECUTED FROM INSIDE THE WORK DIRECTORY FOR THE LF MODEL
export curr_dir=$(pwd)
       
# Name of the executable you want to run
if [ ! -f Outputs/0DOutputVector.dat ]; then
        
    # pre-processing for the input format
    head -n 5 Inputs/MODELNAME_StochasticVector_DAKOTA.dat > Inputs/ao1.dat
    tail -n 4 Inputs/ao1.dat > Inputs/ao2.dat
    cut -c 23-100 Inputs/ao2.dat >  Inputs/MODELNAME_StochasticVector_Intermediary.dat
    rm Inputs/ao1.dat Inputs/ao2.dat

    # running the actual computation
    shopt -s extglob
    head -n 10 Inputs/MODELNAME_StochasticVector_DAKOTA.dat > temp.dat
    mesh="$(tail -n 1 temp.dat)"
    mesh="${mesh##*( )}"
    rm temp.dat
    shopt -u extglob

    if [ "${mesh}" == "-fine -MESHLEVEL" ]; then
        python run.py
        wait
        sbatch run_sim_1D.sh
    elif [ "${mesh}" == "-medium -MESHLEVEL" ]; then
        python run.py
        wait
        sbatch run_sim_0D.sh
    else 
        export dir=$(basename "$PWD")
        export SCRATCH_DIR_LF=/scratch/${USER}/${SLURM_JOB_ID}/${dir}
        mkdir -p ${SCRATCH_DIR_LF} 
        cp -pr * ${SCRATCH_DIR_LF}
        cd ${SCRATCH_DIR_LF}
        python DAKOTA_driver.py -d0D MODELNAME > 'run0DModel.o' 2> 'run0DModel.e'
        wait

        # Put output file where the main dakota job looks for it
        cp -pr ${SCRATCH_DIR_LF}/Outputs/* ${curr_dir}/Outputs/
        
        # Remove the directory we made
        cd ${curr_dir}
        rm -rf ${SCRATCH_DIR_LF}
    fi

    
    cd ${curr_dir}
	maxtime=36000
	SECONDS=0
    
    until [[ -f Outputs/0DOutputVector.dat || $SECONDS -ge $maxtime ]]
    do
      sleep 60s
    done

fi

# head -n 4 Outputs/0DOutputVector.dat > Outputs/0DOutputVector_head.dat
# tail -n 1 Outputs/0DOutputVector_head.dat > Outputs/0DOutputSingleQoI.dat
# rm Outputs/0DOutputVector_head.dat
