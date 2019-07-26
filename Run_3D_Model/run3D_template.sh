#!/bin/bash

# EXECUTED FROM INSIDE THE WORK DIRECTORY FOR THE HF MODEL
export curr_dir=$(pwd)

# pre-processing for the input format
if [ ! -f Outputs/3DOutputVector.dat ]; then

	# pre-processing for the input format
    head -n 5 Inputs/MODELNAME_StochasticVector_DAKOTA.dat > Inputs/ao1.dat
	tail -n 4 Inputs/ao1.dat > Inputs/ao2.dat
	cut -c 23-100 Inputs/ao2.dat >  Inputs/MODELNAME_StochasticVector_Intermediary.dat
	rm Inputs/ao1.dat Inputs/ao2.dat

	# running the actual computation
	python DAKOTA_driver.py -d3D MODELNAME > '3DJob.out' 2> '3DJob.err'


	cd ${curr_dir}
	maxtime=72000
	SECONDS=0
	
	until [[ -f Outputs/3DOutputVector.dat || $SECONDS -ge $maxtime ]]
    do
      sleep 600s
    done

fi

# head -n 4 Outputs/3DOutputVectorNew.dat > Outputs/3DOutputVector_head.dat
# tail -n 1 Outputs/3DOutputVector_head.dat > Outputs/3DOutputSingleQoI.dat
# rm Outputs/3DOutputVector_head.dat

