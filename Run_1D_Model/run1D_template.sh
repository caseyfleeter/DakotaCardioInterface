#!/bin/bash

# EXECUTED FROM INSIDE THE WORK DIRECTORY FOR THE LF MODEL
export curr_dir=$(pwd)
       
# Name of the executable you want to run
if [ ! -f Outputs/1DOutputVector.dat ]; then
        
    # pre-processing for the input format
    head -n 5 Inputs/MODELNAME_StochasticVector_DAKOTA.dat > Inputs/ao1.dat
    tail -n 4 Inputs/ao1.dat > Inputs/ao2.dat
    cut -c 23-100 Inputs/ao2.dat >  Inputs/MODELNAME_StochasticVector_Intermediary.dat
    rm Inputs/ao1.dat Inputs/ao2.dat

    # running the actual computation
    python run.py
    wait
    sbatch run_sim_1D.sh


    cd ${curr_dir}
    maxtime=36000
    SECONDS=0
    
    until [[ -f Outputs/1DOutputVector.dat || $SECONDS -ge $maxtime ]]
    do
      sleep 60s
    done

fi

# head -n 4 Outputs/0DOutputVector.dat > Outputs/0DOutputVector_head.dat
# tail -n 1 Outputs/0DOutputVector_head.dat > Outputs/0DOutputSingleQoI.dat
# rm Outputs/0DOutputVector_head.dat
