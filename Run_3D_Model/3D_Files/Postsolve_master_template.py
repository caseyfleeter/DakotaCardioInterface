# --- VTK-PYTHON SCRIPT FOR POST-PROCESSING
# --- NOTE: THIS SCRIPT ASSUMES A TRIANGULAR SURFACE MESH

import sys
import os
import re
import vtk
import numpy as np
import time
import glob
import pandas as pd

# PATHS TO RELEVENT FILES
SV_POST_PATH = 'svpost'
RESTARTS_PATH = '*[0-9]-procs_case'
results_prefix = 'MODELNAME'

# NUMBER OF TRIALS AND NUMBER OF QOI/OUTPUTS
num_trials = 1
num_outputs = 133    # total QoIs
num_flow = 10        # flow QoIs
num_pressure_1 = 19  # pressure QoIs
num_pressure_2 = 10  # pressure amplitude QoIs
num_tabp = 45        # time-averaged pressure QoIs
num_tawss = 49       # TAWSS QoIs

# SPECIFY TIME POINTS TO DETERMINE QOI, ALONG WITH INCREMENT IN RESTART FILES
START_TIME = 500
END_TIME = 600
INCREMENT = 10

if __name__ == "__main__":
  
  outputVec = np.zeros(num_outputs)

  # check if simulation started
  sim_start_check = glob.glob(RESTARTS_PATH)
  while(len(sim_start_check) == 0):
    # check for error
    errorFile = glob.glob('*.e*')
    if len(errorFile) > 0:  
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 3D solution terminated. See ' + errorFile[0] + ' for details')
  
    # keep waiting
    time.sleep(600)
    sim_start_check = glob.glob(RESTARTS_PATH)

  
  # check if simulation completed
  results_check = glob.glob('all_results.*')
  if(len(results_check) == 0):  
    
    os.chdir(sim_start_check[0])
  
    # Check the numstart.dat to see if this simulation is complete
    numstart_file = open('numstart.dat', 'r')
    last_step = int(numstart_file.readline())
    print('Last completed step: ' + str(last_step))
    numstart_file.close()
    while(last_step < END_TIME):
      # check for error
      errorFile = glob.glob('../*.e*')
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 3D solution terminated. See ' + errorFile[0] + ' for details')

      # Remove restart files 
      for i_restart in xrange(0,START_TIME,INCREMENT):
        command_string = 'rm restart.' + str(i_restart) + "*"
        os.system(command_string)
      
      # keep waiting
      time.sleep(600)
      numstart_file = open('numstart.dat', 'r')
      last_step = int(numstart_file.readline())
      print('Last completed step: ' + str(last_step))
      numstart_file.close()
  
    # Run post solver
    command_string = '../' + SV_POST_PATH + ' -start ' + str(START_TIME) + ' -stop ' + str(END_TIME) + ' -incr ' + \
      str(INCREMENT) + ' -vtp all_results.vtp -vtu all_results.vtu -vtkcombo -all'
    os.system(command_string)

    results_check = glob.glob('all_results.*')
    if(len(results_check) == 0):
      sys.exit('ERROR: 3D solution could not be post-processed')
  
    # Move the all_results files out
    command_string = 'mv all_results* ..'
    os.system(command_string)
  
    os.chdir('..')

    # Remove results folder
    command_string = 'rm -rf ' + sim_start_check[0]
    os.system(command_string)

  postpro_check = glob.glob('all_results_tawss.dat')
  if(len(postpro_check) == 0):
    command_string = 'python Postsolve.py'
    os.system(command_string)

  
  # get average flow at outlets over final cycle
  flows_result = pd.read_csv("all_results_flows.dat",delimiter=' ')
  for i in xrange(0,num_flow):
    col = flows_result[flows_result.columns[i+1]]
    outputVec[i] = col.sum()/col.shape[0] 
  

  # get average pressure at outlets over final cycle
  pressures_result = pd.read_csv("all_results_pressures.dat",delimiter=' ')
  for i in xrange(num_flow,num_flow+num_pressure_1):
    col = pressures_result[pressures_result.columns[i-num_flow+1]]
    outputVec[i] = col.mean() 

  pressures_2_result = pd.read_csv("all_results_pressures_2.dat",delimiter=' ')
  for i in xrange(num_flow+num_pressure_1,num_flow+num_pressure_1+num_pressure_2):
    col = pressures_2_result[pressures_2_result.columns[i-(num_flow+num_pressure_1)]]
    outputVec[i] = col.mean() 


  # get the time-averaged pressures over final cycle
  tabp_result = pd.read_csv("all_results_tabp.dat",delimiter=' ')
  for i in xrange(num_flow+num_pressure_1+num_pressure_2, num_flow+num_pressure_1+num_pressure_2+num_tabp):
    col = tabp_result[tabp_result.columns[i-(num_flow+num_pressure_1+num_pressure_2)]]
    outputVec[i] = col.mean() 


  # get the time-averaged WSS over the final cycle
  tawss_result = pd.read_csv("all_results_tawss.dat",delimiter=' ')
  for i in xrange(num_flow+num_pressure_1+num_pressure_2+num_tabp, num_flow+num_pressure_1+num_pressure_2+num_tabp+num_tawss):
    col = tawss_result[tawss_result.columns[i-(num_flow+num_pressure_1+num_pressure_2+num_tabp)]]
    outputVec[i] = col.mean()


  # Remove copied over scripts
  command_string = 'rm Postsolve.py svpre svpost svsolver'
  os.system(command_string)
  
  # Return to home directory
  os.chdir('../')
    
    
  # Write out the results to file
  os.chdir('Outputs/')
  np.savetxt('3DOutputVector.dat',outputVec,delimiter=',')
  
   # Return to home directory
  os.chdir('../')
    
    
