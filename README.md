# Documentation for Dakota-Simvascular Cardiovascular Interface	

## General Notes
This project contains a framework for the Dakota-Simvascular coupled uncertainty quantification workflow. This project is designed to be customized to fit specific user desires. The framework couples Dakota to the Simvascular cardiovascular flow solvers, and can be adapted to utilize any of Dakota’s many features.

The files provided here were adapted from files for MLMF uncertainty quantification using three mesh levels each of high- and low-fidelity models. This framework is flexible in the sense that any number of mesh levels for high and low fidelity models can be added to the framework

Additionally, by changing the `dakota_MLMF_[3D-1D, 3D-1D-0D]_template.in` file, additional uncertainty quantification techniques, optimization, or parameter estimation techniques implemented in DAKOTA can be used with the Simvascular flow solvers.

## Dependencies not included here
To utilize this framework, the user needs to install the DAKOTA framework from Sandia National Labs (found at [Dakota | Explore and predict with confidence](https://dakota.sandia.gov)) and the Simvascular framework (found at [SimVascular](https://simvascular.github.io) and https://github.com/SimVascular/oneDSolver ). 

It is assumed the user will follow the installation instructions from both of those sites to install the needed executables for this coupled framework. 

## Further Details
More detailed documentation can be found in the `Documentation` directory. 

Below are details for how to specifically edit each provided template file to run simulations.

- - - -
## Outer Directory
### dakota_MLMF_[3D-1D, 3D-1D-0D]_template.in
* Replace line 11 with the desired number of `pilot_samples` per level
* Replace line 13 with the maximum number of iterations for convergence
* Replace line 14 with the desired convergence tolerance (i.e. the improvement factor in the variance of the quantity of interest).
* Replace line 30 with approximate costs of each low-fidelity (LF) model level
* Replace line 39 with approximate costs of each high-fidelity (HF) model level
* Replace lines 43-46 with the uncertain parameters, distributions, and bounds. Uniform uncertainties are shown.
* Add additional mesh levels in lines 47-52
* Define LF simulations in lines 55-65
	* Can change number of concurrent jobs launched by Dakota
* Define HF simulations in lines 68-78
	* Can change number of concurrent jobs launched by Dakota
	* On line 61, replace  `MODELNAME` with the name of your model
* Define output QoIs in lines 80-85
	* Add additional QoIs
	* On line 74, replace  `MODELNAME` with the name of your model
* After editing the file, remove the `3D-1D` or `3D-1D-0D` from the filename (i.e. rename the file `dakota_MLMF_template.in`)

### MakeExecutables_template.sh
* Update lines 6, 10-12, and 17 to include the correct paths to the `OneDSolver`, `svsolver`, `spree`, `svpost`, and `dakota` executables.
* After editing the file, remove the `_template` suffix from the filename

### runDakota_template.sh
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Add your email address to line 26.
* Switch between line 40 and line 41 depending on if you want a 3D-1D-0D or 3D-1D job.
* Switch from line 51 to line 52 if you are restarting your simulation.
* After editing the file, remove the `_template` suffix from the filename


## Run_0D_Model Directory
### DAKOTA_driver_0D_template.py
* 	On line 109, replace `i_parm == #` with the correct line number from `0D_Files/Defaults0DRun-*.dat`. This is the line with the parameter `-NoUncertain` for both the coarse and medium meshes.
* On lines  146 and 154, replace the first index into the `runInputVec` with the index (i.e. line number minus 1) containing the inflow waveform data in `0D_Files/Defaults0DRun-*.dat`  for both the coarse and medium meshes. This is the line with the parameter `-INLETWAVEFORM`.
* After editing the file, remove the `_template` suffix from the filename

### DAKOTA_driver_1D_template.py
* Replace filenames for model inlet on lines 19-20
* Adjust for your quantities of interest (QoIs)
	* Adjust lines 210-224 for flow QoIs
	* Adjust lines 227-262 for pressure QoIs
	* Adjust lines 266-610 for TAWSS and TABP QoIs
* After editing the file, remove the `_template` suffix from the filename

### DAKOTA_driver_template.py
* Adjust lines 26-44 for your specific uncertain parameters
* Change indices of `stochasticVectorIntermediary` throughout this file for your specific uncertain parameters
* Change area percentages on line 34 for your model (to evenly distribute boundary conditions).
* Change line 69 to the names of your branches.
* Change lines 79, 93, and 108 for the time steps of your inlet flow waveform
	* Be careful that line 108 needs the additional comma at the end of the timestep list.
* After editing the file, remove the `_template` suffix from the filename

### run_sim_0D_template.sh
*  This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change model name on line 36
* After editing the file, remove the `_template` suffix from the filename

### run_sim_1D_template.sh
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change model name on line 37
* After editing the file, remove the `_template` suffix from the filename

### run_template.py
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change maximum number of jobs appropriate to your computing environment on line 3
* Add your username to lines 6 and 12
* After editing the file, remove the `_template` suffix from the filename

### run0D_template.sh
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change model name on lines 10, 12, 17, and 37
* Change `5` and `4` on lines 10 and 11 to the number of uncertain parameters plus 2 and 1, respectively.
* Uncomment lines 60-62 to extract a single QoI (if desired). Change line 60 to the specific QoI.
* After editing the file, remove the `_template` suffix from the filename

### Inputs/ Directory
* This is an empty directory which DAKOTA automatically populated with the uncertain parameter values (and flags) in a file called `MODELNAME_StochasticVector.dat`
* Remove the empty file `MODELNAME_StochasticVector_template.dat` before running your simulation.

### Outputs/ Directory
* This is an empty directory which the pipeline automatically populated with the output vector with all quantities of interest for the model. 
* This vector is called `0DOutputVector.dat` for a 0D simulation
* Remove the empty  `0DOutputVector_template.dat`  file being running your simulation.

### 0DFiles/Defaults0DLPN-[coarse, medium]_template.dat
* Fill this file in with the default values of the parameters in the LPN file for your model
* Format is `VALUE:TYPE:-FLAG`
	* `TYPE` can be `INT`, `FLOAT`, `STRING,` or `LIST`. 
	* For uncertain parameters,  `-FLAG` is the as is used in the `dakota_MLMF_3D-1D_template.in` or `dakota_MLMF_3D-1D-0D_template.in` files.
* After editing the file, remove the `_template` suffix from the filename

### 0DFiles/Defaults0DRun-[coarse, medium]_template.dat
* Fill this file in with the default values of the parameters in the Run file for your model
* Format is `VALUE:TYPE:-FLAG`
	* `TYPE` can be `INT`, `FLOAT`, `STRING,` or `LIST`. 
	* For uncertain parameters,  `-FLAG` is the as is used in the `dakota_MLMF_3D-1D_template.in` or `dakota_MLMF_3D-1D-0D_template.in` files.
* The parameter with the flag `-NoUncertain` is used to determine how many uncertain parameters need to be passed to the 0D solver. In this example, there are 30 uncertain parameters (a combination of boundary conditions and material parameters). 
* After editing the file, remove the `_template` suffix from the filename

### 0DFiles/rk4.py
* Runge-Kutta 4th Order differential equation solver
* Can be substituted with the differential equation solver of your choice

### 0DFiles/unscaled_flow_template.dat
* Replace with the default flow waveform for your model.
* Line 1 are time steps; line 2 are flow values.
* After editing the file, remove the `_template` suffix from the filename

### 0DFiles/Templates/mesh-coarse/MODELNAME_LPN_template.tmpl
* Replaces lines 33-58 with the correct parameters for your model
* Replaces lines 64-98 with the correct linear system for your model (in the resistor only case with no differential equations)
* Replaces lines 107-145 with the quantities of interest you need for your model.
* Replace lines 169 and 189 with the index of the blood viscosity in the `QI` vector defined in `MODELNAME_run_template.tmpl` (if this is an uncertain parameter - otherwise can be hardcoded).
* Replace line 187 with the index of the Young’s modulus in the `QI` vector defined in `MODELNAME_run_template.tmpl` (if this is an uncertain parameter - otherwise can be hardcoded).
* Replace line 188 with the index of the blood density in the `QI` vector defined in `MODELNAME_run_template.tmpl` (if this is an uncertain parameter - otherwise can be hardcoded).
* Replace lines 191-225 with the geometry information for your model.
* Replace lines 252-272 with the parameters need for your model.
* Replace `MODELNAME` in filename with the name of your model
* After editing the file, remove the `_template` suffix from the filename

### 0DFiles/Templates/mesh-coarse/MODELNAME_run_template.tmpl
* Replace MODELNAME on line 7 with the name of your model.
* Adjust for your quantities of interest (QoIs)
	* Adjust lines 56-65 for flow QoIs
	* Adjust lines 72-136 for pressure QoIs
	* Adjust lines 184-416 for TAWSS and TABP QoIs
* Remove line 146 if not using WSS thresholds
* Replace lines 148-182 with the geometry information for your model.
* Replace `MODELNAME` in filename with the name of your model
* After editing the file, remove the `_template` suffix from the filename

### 0DFiles/Templates/mesh-medium/MODELNAME_LPN_template.tmpl
* Replaces lines 42-109 with the quantities of interest you need for your model.
* Replace lines 132  and 223 with the index of the blood viscosity in the `QI` vector defined in `MODELNAME_run_template.tmpl` (if this is an uncertain parameter - otherwise can be hardcoded).
* Replace lines 144-209 with the system of equations for your model.
* Replace line 221 with the index of the Young’s modulus in the `QI` vector defined in `MODELNAME_run_template.tmpl` (if this is an uncertain parameter - otherwise can be hardcoded).
* Replace line 222 with the index of the blood density in the `QI` vector defined in `MODELNAME_run_template.tmpl` (if this is an uncertain parameter - otherwise can be hardcoded).
* Replace lines 225-259 with the geometry information for your model.
* Replace lines 285-358 with the parameters need for your model.
* Replace `MODELNAME` in filename with the name of your model
* After editing the file, remove the `_template` suffix from the filename

### 0DFiles/Templates/mesh-medium/MODELNAME_run_template.tmpl
* Replace MODELNAME on line 7 with the name of your model.
* Adjust for your quantities of interest (QoIs)
	* Adjust lines 54-64 for flow QoIs
	* Adjust lines 72-135 for pressure QoIs
	* Adjust lines 184-416 for TAWSS and TABP QoIs
* Remove line 146 if not using WSS thresholds
* Replace lines 148-182 with the geometry information for your model.
* Replace `MODELNAME` in filename with the name of your model
* After editing the file, remove the `_template` suffix from the filename


## Run_1D_Model Directory
### DAKOTA_driver_1D_template.py
* Replace filenames for model inlet on lines 19-20
* Adjust for your quantities of interest (QoIs)
	* Adjust lines 210-224 for flow QoIs
	* Adjust lines 227-262 for pressure QoIs
	* Adjust lines 266-610 for TAWSS and TABP QoIs
* After editing the file, remove the `_template` suffix from the filename

### DAKOTA_driver_template.py
* Adjust lines 26-44 for your specific uncertain parameters
* Change indices of `stochasticVectorIntermediary` throughout this file for your specific uncertain parameters
* Change area percentages on line 34 for your model (to evenly distribute boundary conditions).
* Change line 69 to the names of your branches.
* Change lines 79, 93, and 108 for the time steps of your inlet flow waveform
	* Be careful that line 108 needs the additional comma at the end of the timestep list.
* After editing the file, remove the `_template` suffix from the filename

### run_sim_1D_template.sh
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change model name on line 37
* After editing the file, remove the `_template` suffix from the filename

### run_template.py
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change maximum number of jobs appropriate to your computing environment on line 3
* Add your username to lines 6 and 12
* After editing the file, remove the `_template` suffix from the filename

### run1D_template.sh
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change model name on lines 10 and 12
* Change `5` and `4` on lines 10 and 11 to the number of uncertain parameters plus 2 and 1, respectively.
* Uncomment lines 32-34 to extract a single QoI (if desired). Change line 32 to the specific QoI.
* After editing the file, remove the `_template` suffix from the filename

### Inputs/ Directory
* This is an empty directory which DAKOTA automatically populated with the uncertain parameter values (and flags) in a file called `MODELNAME_StochasticVector.dat`
* Remove the empty file `MODELNAME_StochasticVector_template.dat` before running your simulation.

### Outputs/ Directory
* This is an empty directory which the pipeline automatically populated with the output vector with all quantities of interest for the model. 
* This vector is called `1DOutputVector.dat` for a 1D simulation
* Remove the empty  `1DOutputVector_template.dat`  file being running your simulation.

### 1DFiles/Defaults1D_template.dat
* Fill this file in with the default values of the parameters in your runfile for your model.
* Format is `VALUE:TYPE:-FLAG`
	* `TYPE` can be `INT`, `FLOAT`, `STRING,` or `LIST`. 
	* For uncertain parameters,  `-FLAG` is the as is used in the `dakota_MLMF_3D-1D_template.in` or `dakota_MLMF_3D-1D-0D_template.in` files.
* After editing the file, remove the `_template` suffix from the filename

### 1DFiles/QoI-[Flow, Pressure, TA, TABranches, TAStrips, WSSThresholds]_template.dat
* Fill these file with specific geometry information for your model
	* `QoI-Flow` lists the file name  (branch names and segment numbers) on which to calculate inlet and outlet flows.
	* `QoI-Pressure` lists the file name  (branch names and segment numbers) on which to calculate inlet and outlet pressures.
	* `QoI-TA` lists the branch names and total number of segments for calculating time-averaged quantities (pressures and WSS).
	*  `QoI-TABranches` lists the branch names and the segment range for calculating time-averaged quantities (pressures and WSS) on a segment of the branch.
	*  `QoI-TAStrips` lists the branch names and the segment range for calculating time-averaged quantities (pressures and WSS) on strips of interest.
	*  `QoI-WSSThresholds` lists the thresholds used for calculating areas of low TAWSS in the model.
* After editing the files, remove the `_template` suffix from the filename

### 1DFiles/unscaled_flow_template.dat
* Replace with the default flow waveform for your model.
* Line 1 are time steps; line 2 are flow values.
* After editing the file, remove the `_template` suffix from the filename

### 1DFiles/Templates/MODELNAME-[coarse, fine, medium]_template.dat
* Replace these files with the template file for each mesh level of your model.
* These templates can be automatically generated from the SimVascular 3D model.
* After editing the file, remove the `_template` suffix from the filename


## Run_3D_Model Directory
### DAKOTA_driver_3D_template.py
* Replace filename of the model’s inlet cap on line 12
* After editing the file, remove the `_template` suffix from the filename

### DAKOTA_driver_template.py
* Adjust lines 26-44 for your specific uncertain parameters
* Change indices of `stochasticVectorIntermediary` throughout this file for your specific uncertain parameters
* Change area percentages on line 34 for your model (to evenly distribute boundary conditions).
* Change line 69 to the names of your branches.
* Change lines 79, 93, and 108 for the time steps of your inlet flow waveform
	* Be careful that line 108 needs the additional comma at the end of the timestep list.
* After editing the file, remove the `_template` suffix from the filename

### run3D_template.sh
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Change model name on lines 10, 12, and 16
* Change `5` and `4` on lines 10 and 11 to the number of uncertain parameters plus 2 and 1, respectively.
* Uncomment lines 30-32 to extract a single QoI (if desired). Change line 30 to the specific QoI.
* After editing the file, remove the `_template` suffix from the filename

### Inputs/ Directory
* This is an empty directory which DAKOTA automatically populated with the uncertain parameter values (and flags) in a file called `MODELNAME_StochasticVector.dat`
* Remove the empty file `MODELNAME_StochasticVector_template.dat` before running your simulation.

### Outputs/ Directory
* This is an empty directory which the pipeline automatically populated with the output vector with all quantities of interest for the model. 
* This vector is called `3DOutputVector.dat` for a 1D simulation
* Remove the empty  `3DOutputVector_template.dat`  file being running your simulation.

### 3DFiles/Defaults3DBct_template.dat
* Fill this file in with the default value of the `-INLETWAVEFORM` for your model.
* Keep in mind that the 3D svsolver requires flow values to be inverted compared to the 0D and 1D model.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Defaults3DRcrt_template.dat
* Fill this file in with the default value of the outlet boundary conditions for your model.
* RCR boundary conditions are shown here, but any boundary condition compatible with Simvascular can be used instead; simply supply default values for the needed parameters here.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Defaults3DSolver_template.dat
* Fill this file in with the default values of the parameters for the `solver.inp` file for your model.
* If resistance boundary conditions are used, this is where the default values will be specified (as shown).
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Defaults3DSvpre_template.dat
* Fill this file in with the default values of the parameters for the `MODELNAME.svpre` file for your model..
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Postsolve_master_template.py
* Replace line 16 with the name of your model.
* Replace lines 19-25 with the details for your quantities of interest.
* Replace lines 28-30 with the number of time steps for the restart files of your simulation.
* Edit lines 113-136 for the quantities of interest needed for your simulation.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Postsolve_template.py
* Replace lines 13-15 with the number of time steps for the restart files of your simulation.
* Replace line 18 with the names of the geometry files for your model.
* Replace line 19 with the name of the inlet cap for your model.
* Replace lines 22-39 with the point and normal vectors defining clipping planes of the full model to extract branches and strips of interest (if desired). This information can be determined by clipping the model in Paraview.
* Replace line 42 with the thresholds for areas of low TAWSS in the model, if desired.
* Edit lines 44-65 for the output files for the quantities of interest needed for your simulation.
* Edit lines 940-1113 to remove any QoIs not in use for your simulation.
* Remove any sections in the remainder of the file for QoIs not used in your simulation.  
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/unscaled_flow_template.dat
* Replace with the default flow waveform for your model.
* Line 1 are time steps; line 2 are flow values.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/mesh-[coarse, fine, medium]/mesh-complete/ Directory
* Replace this directory with the `mesh-complete` directory for each mesh level of your 3D model.
* Remove the empty file `mesh-complete/mesh-complete.exterior_template.vtp` before running your simulation.

### 3DFiles/mesh-[coarse, fine, medium]/Solver_Files/ Directory
* Replace the files in this directory with the proper files for each mesh level of your 3D model needed to run the simulations.
* These following files are only needed if your 3D simulations will have deformable walls:
	* `geombc_dat.1_template` should be replaced with the `geombc.dat.1` file generated from an initial rigid simulation with default uncertain parameter values .
	* `restart.0.1_template` should be replaced with the `restart.0.1` file generated from the final timestep of an initial rigid simulation with default uncertain parameter values .
	* After replacing the files, make sure to remove the `_template` suffix from the filenames
* These files are needed for rigid and deformable wall simulations:
	* `numstart.dat` can remain unchanged

### 3DFiles/mesh-[coarse, fine, medium]/Solver_Files/run.job.sh_template
* This job script is for a SLURM cluster - please change as needed for other infrastructures. Also, update the partition names and node configuration as appropriate for your cluster.
* Add your email address to line 24.
* Change model name on line 43.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Templates/MODELNAME-bct_template.tmpl
* Edit this file for additional parameters for a non-inflow waveform inlet boundary condition.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Defaults3DRcrt_template.dat
* Fill this file in with a templated value of the outlet boundary condition file for your model.
* RCR boundary conditions are shown here, but any boundary condition compatible with Simvascular can be used instead; simply supply the templated file here.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Defaults3DSolver_template.dat
* Adjust this file to match the `solver.inp` file for your model.
* After editing the file, remove the `_template` suffix from the filename

### 3DFiles/Defaults3DSvpre_template.dat
* Adjust this file to match the `MODELNAME.svpre` file for your model.
* After editing the file, remove the `_template` suffix from the filename

- - - - 
Please cite this work with [![DOI](https://zenodo.org/badge/170943673.svg)](https://zenodo.org/badge/latestdoi/170943673)
