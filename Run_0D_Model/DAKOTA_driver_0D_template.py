import sys,os
import re
import numpy as np
import math
import argparse
import glob
import time
import datetime
import ast



def Driver_0D_get_inputs(defaultRunVec,defaultLPNVec,stochasticVec):

  # Code executes while in OD_Files folder

  # iterate through all stochastic values and replace the default value
  mergedRunVec = defaultRunVec
  mergedLPNVec = defaultLPNVec

  numStoch = len(stochasticVec)
  stochasticVec = stochasticVec[0:numStoch-1]
  for i_stoch in stochasticVec:
    i_stoch = i_stoch.rsplit(' ')
    i_data = i_stoch[0]
    i_flag = i_stoch[1]

    # get index of parameter to be replaced
    flag_found_run = [i_flag in m for m in mergedRunVec]
    i_replace_run = [i for i, found in enumerate(flag_found_run) if found == True]

    flag_found_LPN = [i_flag in m for m in mergedLPNVec]
    i_replace_LPN = [i for i, found in enumerate(flag_found_LPN) if found == True]

    # error check
    num_found = len(i_replace_run) + len(i_replace_LPN)
    if num_found == 0:
      error_message = 'ERROR: invalid flag ' + i_flag + ' provided in stochastic vector'
      sys.exit(error_message)
    if num_found > 1:
      error_message = 'ERROR: flag ' + i_flag + ' appears multiple times in defaults vector'
      sys.exit(error_message)

    # replace default values in run file
    if len(i_replace_run) == 1:
      if (mergedRunVec[i_replace_run[0]][1] == 'STRING') or (mergedRunVec[i_replace_run[0]][1] == 'INT') or (mergedRunVec[i_replace_run[0]][1] == 'FLOAT') or (mergedRunVec[i_replace_run[0]][1] == 'LIST'):
        mergedRunVec[i_replace_run[0]][0] = i_data
    
    # replace default values in LPN file
    elif len(i_replace_LPN) == 1:
      if (mergedLPNVec[i_replace_LPN[0]][1] == 'STRING') or (mergedLPNVec[i_replace_LPN[0]][1] == 'INT') or (mergedLPNVec[i_replace_LPN[0]][1] == 'FLOAT') or (mergedLPNVec[i_replace_LPN[0]][1] == 'LIST'):
        mergedLPNVec[i_replace_LPN[0]][0] = i_data
    
  # save the input vector
  os.chdir('../Inputs/')
  saveFile = open('0DRunVector.dat','w')
  saveFile.write('\n'.join([':'.join(m) for m in mergedRunVec]))
  saveFile.close()
  saveFile2 = open('0DLPNVector.dat','w')
  saveFile2.write('\n'.join([':'.join(m) for m in mergedLPNVec]))
  saveFile2.close()

  # return to 0D_Files folder
  os.chdir('../0D_Files/')

  return mergedRunVec, mergedLPNVec


def Driver_0D_fill_template(runInputVec,LPNInputVec,templateFileName,mesh):

  # Code executes while in OD_Files folder

  # Generate Realization File: create 0D solver input file for this parameter vector
  # Open template file
  if '.tmpl' not in templateFileName:
    templateRunName = templateFileName + '_run.tmpl'
    templateLPNName = templateFileName + '_LPN.tmpl'
  else:
    templateFileNameShort = templateFileName.split('.')[0]
    templateRunName = templateFileNameShort + '_run.tmpl'
    templateLPNName = templateFileNameShort + '_LPN.tmpl'
  templateRun = open('Templates/mesh' + mesh + '/' + templateRunName,'r').read()
  templateLPN = open('Templates/mesh' + mesh + '/' + templateLPNName,'r').read()
  
  # Get dimensions of template files
  parIntRun = re.findall(r'<params,(.*?)>', templateRun)
  maxIntRun = np.max([int(p) for p in parIntRun]) + 1
  parIntLPN = re.findall(r'<params,(.*?)>', templateLPN)
  maxIntLPN = np.max([int(p) for p in parIntLPN]) + 1

  # Check if the dimensions are compatible with parameter vectors
  if (maxIntRun != len(runInputVec)):
      print('Dimensions in Template: ' + str(maxIntRun))
      print('Dimensions in Vector: ' + str(len(runInputVec)))
      print('Run template file dimensions not compatible. Terminating.')
      sys.exit(-1)
  if (maxIntLPN != len(LPNInputVec)):
      print('Dimensions in Template: ' + str(maxIntLPN))
      print('Dimensions in Vector: ' + str(len(LPNInputVec)))
      print('LPN template file dimensions not compatible. Terminating.')
      sys.exit(-1)

  # Replace parameters in run file
  for i_param in parIntRun:
    i_param = int(i_param)

    param = runInputVec[i_param]
    # generate uncertian parameter vector
    if (mesh == '-medium' and i_param == 5) or (mesh == '-coarse' and i_param == 4):
      numBCParam = int(param[0])
      # adding multiple items
      replace_string = '['
      for i_data in xrange(i_param + 1, i_param + numBCParam + 1):
        new_param = runInputVec[i_data]
        replace_string = replace_string + str(new_param[0]) + ',' 
      replace_string = replace_string + ']'
      templateRun = templateRun.replace('<params,' + str(i_param) + '>', replace_string )
    else:
      templateRun = templateRun.replace('<params,' + str(i_param) + '>', param[0] )


  # Replace parameters in LPN file
  for i_param in parIntLPN:
    i_param = int(i_param)

    param = LPNInputVec[i_param]
    templateLPN = templateLPN.replace('<params,' + str(i_param) + '>', param[0] ) 

  # Save the realization file in the Inputs Folder
  os.chdir('../Inputs')

  if '.tmpl' in templateFileName:
    templateFileName = templateFileName.split('.')[0]
  runFileName = templateFileName + mesh + '_run.py'
  runFile = open(runFileName,"w+")
  runFile.write(templateRun)
  runFile.close()

  LPNFileName = templateFileName + '_LPN.py'
  LPNFile = open(LPNFileName,"w+")
  LPNFile.write(templateLPN)
  LPNFile.close()

  # save inflow file
  if (mesh == '-medium'):
    string = runInputVec[36][0].split('[')
    times = string[2].split(',')[:-2]
    times = [float(t) for t in times]
    values = string[3].split(',')[:-1]
    values = [float(v) for v in values]
    flow = np.asarray(times + values).reshape(2,len(times))
    np.savetxt('../Inputs/0DInflow.flow',flow,fmt='%.8e',delimiter=',') 
  elif (mesh == '-coarse'):
    string = runInputVec[35][0].split('[')
    times = string[2].split(',')[:-2]
    times = [float(t) for t in times]
    values = string[3].split(',')[:-1]
    values = [float(v) for v in values]
    flow = np.asarray(times + values).reshape(2,len(times))
    np.savetxt('../Inputs/0DInflow.flow',flow,fmt='%.8e',delimiter=',') 
    

  # return to home directory
  os.chdir('../')



def Driver_0D_launch_sim(fileName,mesh):

  # Code executes from within /this_run folder

  # make new directory for this run
  now = datetime.datetime.utcnow()
  dirName = 'Run_' + now.strftime("%s")
  os.mkdir(dirName)

  # Move into new folder
  os.chdir(dirName)

  if '.tmpl' in fileName:
    fileName = fileName.split('.')[0]
  runFileName = fileName + mesh + '_run.py'
  LPNFileName = fileName + '_LPN.py'

  # Move necessary files into this folder
  command_string = 'cp ../Inputs/' + runFileName + ' ../Inputs/' + LPNFileName + ' ../Inputs/0DInflow.flow ../0D_Files/rk4.py .'
  os.system(command_string)

  # Run the simulations
  command_string = 'python ' + runFileName + ' > 0DModel.out 2> 0DModel.err'
  os.system(command_string)


  # Check if simulation has started
  outputFile = glob.glob('0DModel.o*')
  while(len(outputFile) == 0):
    # check for error
    errorFile = glob.glob('0DModel.e*')
    if len(errorFile) > 0:
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 0D solution terminated. See ' + errorFile[0] + ' for details')

    # keep waiting
    time.sleep(30)
    outputFile = glob.glob('0DModel.o*')

  # Check if simulation is complete
  completed_check = open(outputFile[0]).readlines()
  while(len(completed_check) == 0):
    # check for error
    errorFile = glob.glob('0DModel.e*')
    if len(errorFile) > 0:
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 0D solution terminated. See ' + errorFile[0] + ' for details')

    # keep waiting
    time.sleep(10)
    completed_check = open(outputFile[0]).readlines()

  # Get final line
  completed_check = open(outputFile[0]).readlines()[-1]
  while(completed_check != 'Completed!\n'):
    # check for error
    errorFile = glob.glob('0DModel.e*')
    if len(errorFile) > 0:
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 0D solution terminated. See ' + errorFile[0] + ' for details')

    # keep waiting
    time.sleep(10)
    completed_check = open(outputFile[0]).readlines()[-1]

  print("0D SLURM job complete")

  # return to home directory
  os.chdir('../')

  return dirName


# 0D Driver
def Driver_0D(stochasticVector,templateFileName):
  
  print("Made it to Driver_0D")
  os.chdir('0D_Files')

  mesh = stochasticVector[-1].strip().split(' ')[0]

  # get default run vector values
  defaultRunVector = open('Defaults0DRun' + mesh + '.dat').read().splitlines()
  defaultRunVector = [d.split(':') for d in defaultRunVector]

  # get default LPN vector values
  defaultLPNVector = open('Defaults0DLPN' + mesh + '.dat').read().splitlines()
  defaultLPNVector = [d.split(':') for d in defaultLPNVector]


  # replace stochastic values
  RunInputVec, LPNInputVec = Driver_0D_get_inputs(defaultRunVector,defaultLPNVector,stochasticVector)
  print("0D Inputs obtained")

  # fill template file
  Driver_0D_fill_template(RunInputVec,LPNInputVec,templateFileName,mesh)
  print("0D Templates filled")

  # Run the simulation
  runDirName = Driver_0D_launch_sim(templateFileName,mesh)
  print("0D Simulation complete")

  print("0D Results obtained") 



