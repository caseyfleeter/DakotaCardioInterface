This project includes python scripts to interface between DAKOTA and the cardiovascular solvers.


#In the Inputs/ Folder:
	
	- [modelName]_StochasticVector.tmpl needs to be replaces with the file of stochastic parameter values (and flags) called [modelName]_StochasticVector.dat, where [modelName] is replaced with the name of the model (also used in the template file)
	- This file is a two column file with each line containing the stochastic value for a parameter followed by that parameters flag. This allows for the value to be inserted in the vector of default parameter values
	- For the 1D and 3D solvers, the [modelName]_StochasticVector.dat file contains the mesh level (-coarse, -medium or -fine) on the last line
	- And example file is provided:

		PARAMETER_STOCH_VALUE 	FLAG  # DO NOT INCLUDE THIS HEADER LINE IN YOUR FILE
		6058.52570 				-BC0
		24933.83724 			-BC1
		15911.40456 			-BC2
		15911.146007			-BC3
		20582.482622 			-BC4
		24933.94874 			-BC5
		7036.31991 				-BC6
		20581.48207 			-BC7
		11797.32087 			-BC8
		-coarse 				-MESHLEVEL

#In the Outputs/ Folder:
	
	- This empty folder with contain the output vector with all quantities of interest for the model. 
	- This vector is called:
		- 0DOutputVector.dat for a 0D simulation
		- 1DOutputVector.dat for a 1D simulation
		- 3DOutputVector.dat for a 3D simulation


	IN THE FUTURE, this file with contain a flag describing the QoI, but currently is just a plain vector.




#To run a 0D simulation:

##In the home folder for this simulation:
	
	- DAKOTA_driver.py needs to be included
	- DAKOTA_driver_0D.py needs to be included
	- 0D_Files/ folder needs to be included
	- Inputs/ folder needs to be included
	- Outputs/ folder needs to be included (an empty directory)

##In the Inputs/ Folder:
	
	- [modelName]_StochasticVector.tmpl needs to be replaces with the file of stochastic parameter values (and flags) called [modelName]_StochasticVector.dat, where [modelName] is replaced with the name of the model (also used in the template file)


##In the 0D_Files/ Folder:

	- Defaults0D.tmpl needs to replaced with the file of default parameter values called Defaults0D.dat
	- QoI-Branches.tmpl needs to replaced with a file containing the branch names called QoI-Branches.dat
	- QoI-WSSThresholds.tmpl needs to replaced with a file containing the thresholds for computing areas of low WSS called QoI-WSSThresholds.dat
	- rk4.py included but NOT modifed
	- runLPN.py included but NOT modified
	- Templates/ folder needs to be provided


##In the OD_Files/Templates/ Folder:
	
	- [modelName].tmpl needs to replaces with the template file for this model called [modelName].tmpl, where [modelName] is replaces with the name of the model (also used in the stochastic vector file)


## To run:

From the home folder for this simulation, the command line argument is

	python DAKOTA_driver.py -d0D [modelName]


## Files generated

After running the simulation successfully, the following files are added:

	- In Inputs/ the file 0DInputVector.dat containing the merged input vector (replacing default values with the stochastic values) is saved.
	- In Inputs/ the file LPN_Full.py is saved. This is the filled in template file for this simulation.
	- In Outputs/ the file 0DOutputVector.dat containing the output quantities of interest is saved.
	- In the home directory, a folder Run_[#]/ is added. This contains four files with the simulation results. These files contain headers describing each column.
		- all_results_flows.dat
		- all_results_pressures.dat
		- all_results_wss_branches.dat
		- all_results_wss_model.dat




#To run a 1D simulation:

##In the home folder for this simulation:
	
	- DAKOTA_driver.py needs to be included
	- DAKOTA_driver_1D.py needs to be included
	- 1D_Files/ folder needs to be included
	- Inputs/ folder needs to be included
	- Outputs/ folder needs to be included (an empty directory)

##In the Inputs/ Folder:
	
	- [modelName]_StochasticVector.tmpl needs to be replaces with the file of stochastic parameter values (and flags) called [modelName]_StochasticVector.dat, where [modelName] is replaced with the name of the model (also used in the template file). The desired mesh level (-coarse, -medium or -fine) needs to be provided on the last line.


##In the 1D_Files/ Folder:

	- Defaults1D.tmpl needs to replaced with the file of default parameter values called Defaults1D.dat
	- QoI-Flow.tmpl needs to replaced with a file containing the branch names where flow is to be computed, called QoI-Flow.dat
	- QoI-Pressure.tmpl needs to replaced with a file containing the branch names where pressure is to be computed, called QoI-Pressure.dat
	- QoI-WSS.tmpl needs to replaced with a file containing the names and number of segements in each branch, called QoI-WSS.dat
	- QoI-WSSThresholds.tmpl needs to replaced with a file containing the thresholds for computing areas of low WSS called QoI-WSSThresholds.dat
	- OneDSolver needs to be replaced with the executable for the 1DSolver (or a link to the executable)
	- Templates/ folder needs to be provided


##In the 1D_Files/Templates/ Folder:
	
	- [modelName]-coarse.tmpl needs to replaced with the template file for the COARSE MESH of this model called [modelName]-coarse.tmpl, where [modelName] is replaces with the name of the model (also used in the stochastic vector file)
	- [modelName]-medium.tmpl needs to replaced with the template file for the MEDIUM MESH of this model called [modelName]-medium.tmpl, where [modelName] is replaces with the name of the model (also used in the stochastic vector file)
	- [modelName]-fine.tmpl needs to replaced with the template file for the FINE MESH of this model called [modelName]-fine.tmpl, where [modelName] is replaces with the name of the model (also used in the stochastic vector file)
	- The code currently only support three mesh levels, but this can be adapted in the future


## To run:

From the home folder for this simulation, the command line argument is

	python DAKOTA_driver.py -d1D [modelName] 


## Files generated

After running the simulation successfully, the following files are added:

	- In Inputs/ the file 1DInputVector.dat containing the merged input vector (replacing default values with the stochastic values) is saved.
	- In Inputs/ the file [modelName]-[mesh].in is saved. This is the filled-in input file for this simulation.
	- In Outputs/ the file 1DOutputVector.dat containing the output quantities of interest is saved.
	- In the home directory, a folder Run_[#]/ is added. This contains files with the simulation results. These files are described in more detail in the OneDSolver documentation but there are: 
		- solver.out contains the STDOUT output from the OneDSolver 
		- solver.err contains any error messages sent to STDERR from the OneDSolver
		- echo.out contains a reformatted version of the OneDSolver input file 
		- *_area.dat: files containing the vessel area at each (saved) timestep for each segment of the model
		- *_flow.dat: files containing the flow rate at each (saved) timestep for each segment of the model
		- *_pressure.dat: files containing the pressure at each (saved) timestep for each segment of the model
		- *_Re.dat: files containing the Reynold's number at each (saved) timestep for each segment of the model
		- *_wss.dat: files containing the WSS at each (saved) timestep for each segment of the model

