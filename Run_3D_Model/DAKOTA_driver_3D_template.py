import sys,os
import re
import numpy as np
import math
import argparse
import glob
import time
import datetime
import ast
import subprocess

INLET_FLOW_CAP = "cap_aorta.flow"


def Driver_3D_get_inputs(defaultSvpreVec,defaultSolverVec,defaultRcrtVec,defaultBctVec,stochasticVec):

  # Code executes while in 3D_Files folder

  # iterate through all stochastic values and replace the default value
  mergedSvpreVec = defaultSvpreVec
  mergedSolverVec = defaultSolverVec
  mergedRcrtVec = defaultRcrtVec
  mergedBctVec = defaultBctVec

  numStoch = len(stochasticVec)
  stochasticVec = stochasticVec[0:numStoch-1]
  for i_stoch in stochasticVec:
    i_stoch = i_stoch.rsplit(' ')
    i_data = i_stoch[0]
    i_flag = i_stoch[1]

    # get index of parameter to be replaced
    flag_found_svpre = [i_flag in m for m in mergedSvpreVec]
    i_replace_svpre = [i for i, found in enumerate(flag_found_svpre) if found == True]

    flag_found_solver = [i_flag in m for m in mergedSolverVec]
    i_replace_solver = [i for i, found in enumerate(flag_found_solver) if found == True]

    flag_found_rcrt = [i_flag in m for m in mergedRcrtVec]
    i_replace_rcrt = [i for i, found in enumerate(flag_found_rcrt) if found == True]

    flag_found_bct = [i_flag in m for m in mergedBctVec]
    i_replace_bct = [i for i, found in enumerate(flag_found_bct) if found == True]

    # error check
    num_found = len(i_replace_solver) + len(i_replace_svpre) + len(i_replace_rcrt) + len(i_replace_bct)
    if num_found == 0:
      error_message = 'ERROR: invalid flag ' + i_flag + ' provided in stochastic vector'
      sys.exit(error_message)
    if num_found > 1:
      error_message = 'ERROR: flag ' + i_flag + ' appears multiple times in default vectors'
      sys.exit(error_message)

    # replace default values in svpre
    if len(i_replace_svpre) == 1:
      if (mergedSvpreVec[i_replace_svpre[0]][1] == 'STRING') or (mergedSvpreVec[i_replace_svpre[0]][1] == 'INT') or (mergedSvpreVec[i_replace_svpre[0]][1] == 'FLOAT') or (mergedSvpreVec[i_replace_svpre[0]][1] == 'LIST'):
        mergedSvpreVec[i_replace_svpre[0]][0] = i_data
    
    # replace default values in solver
    elif len(i_replace_solver) == 1:
      if (mergedSolverVec[i_replace_solver[0]][1] == 'STRING') or (mergedSolverVec[i_replace_solver[0]][1] == 'INT') or (mergedSolverVec[i_replace_solver[0]][1] == 'FLOAT') or (mergedSolverVec[i_replace_solver[0]][1] == 'LIST'):
        mergedSolverVec[i_replace_solver[0]][0] = i_data

    # replace default values in rcrt
    elif len(i_replace_rcrt) == 1:
      if (mergedRcrtVec[i_replace_rcrt[0]][1] == 'STRING') or (mergedRcrtVec[i_replace_rcrt[0]][1] == 'INT') or (mergedRcrtVec[i_replace_rcrt[0]][1] == 'FLOAT') or (mergedRcrtVec[i_replace_rcrt[0]][1] == 'LIST'):
        mergedRcrtVec[i_replace_rcrt[0]][0] = i_data

    # replace default values in bct
    elif len(i_replace_bct) == 1:
      if (mergedBctVec[i_replace_bct[0]][1] == 'STRING') or (mergedBctVec[i_replace_bct[0]][1] == 'INT') or (mergedBctVec[i_replace_bct[0]][1] == 'FLOAT') or (mergedBctVec[i_replace_bct[0]][1] == 'LIST'):
        mergedBctVec[i_replace_bct[0]][0] = i_data

  # save the input vectors
  os.chdir('../Inputs/')
  saveFile = open('3DSvpreVector.dat','w')
  saveFile.write('\n'.join([':'.join(m) for m in mergedSvpreVec]))
  saveFile.close()
  saveFile2 = open('3DSolverVector.dat','w')
  saveFile2.write('\n'.join([':'.join(m) for m in mergedSolverVec]))
  saveFile2.close()
  saveFile3 = open('3DRcrtVector.dat','w')
  saveFile3.write('\n'.join([':'.join(m) for m in mergedRcrtVec]))
  saveFile3.close()
  saveFile4 = open('3DBctVector.dat','w')
  saveFile4.write('\n'.join([':'.join(m) for m in mergedBctVec]))
  saveFile4.close()

  # return to 3D_Files folder
  os.chdir('../3D_Files/')

  return mergedSvpreVec, mergedSolverVec, mergedRcrtVec, mergedBctVec 



def Driver_3D_fill_template(svpreInputVec,solverInputVec,rcrtInputVec,bctInputVec,templateFileName,mesh):

  # Code executes while in 3D_Files folder

  # Generate Realization File: create 3D solver input files for this parameter vector
  # Open template files
  if '.tmpl' not in templateFileName:
    templateSvpreName  = templateFileName + '-svpre' + '.tmpl'
    templateSolverName = templateFileName + '-solver' + '.tmpl'
    templateRcrtName   = templateFileName + '-rcrt' + '.tmpl'
    templateBctName    = templateFileName + '-bct' + '.tmpl'
  else:
    templateFileNameShort = templateFileName.split('.')[0]
    templateSvpreName  = templateFileNameShort + '-svpre' + '.tmpl'
    templateSolverName = templateFileNameShort + '-solver' + '.tmpl'
    templateRcrtName   = templateFileNameShort + '-rcrt' + '.tmpl'
    templateBctName    = templateFileNameShort + '-bct' + '.tmpl'
  templateSvpre  = open('Templates/' + templateSvpreName,'r').read()
  templateSolver = open('Templates/' + templateSolverName,'r').read()
  templateRcrt   = open('Templates/' + templateRcrtName,'r').read()
  templateBct    = open('Templates/' + templateBctName,'r').read()
  
  # Get dimensions of template files
  parIntSvpre  = re.findall(r'<params,(.*?)>', templateSvpre)
  dimsSvpre    = len(np.unique(parIntSvpre))
  parIntSolver = re.findall(r'<params,(.*?)>', templateSolver)
  dimsSolver   = len(np.unique(parIntSolver))
  parIntRcrt   = re.findall(r'<params,(.*?)>', templateRcrt)
  dimsRcrt     = len(np.unique(parIntRcrt))
  parIntBct    = re.findall(r'<params,(.*?)>', templateBct)
  dimsBct      = len(np.unique(parIntBct))

  # Check if the dimensions are compatible with parameter vectors
  if (dimsSvpre != len(svpreInputVec)):
      print 'Dimensions in Template: ' + str(dimsSvpre)
      print 'Dimensions in Vector: ' + str(len(svpreInputVec))
      print 'Svpre template file dimensions not compatible. Terminating.'
      sys.exit(-1)
  if (dimsSolver != len(solverInputVec)):
      print 'Dimensions in Template: ' + str(dimsSolver)
      print 'Dimensions in Vector: ' + str(len(solverInputVec))
      print 'Solver.inp template file dimensions not compatible. Terminating.'
      sys.exit(-1)
  if (dimsRcrt != len(rcrtInputVec)):
      print 'Dimensions in Template: ' + str(dimsRcrt)
      print 'Dimensions in Vector: ' + str(len(rcrtInputVec))
      print 'rcrt.dat template file dimensions not compatible. Terminating.'
      sys.exit(-1)
  if (dimsBct != len(bctInputVec)):
      print 'Dimensions in Template: ' + str(dimsBct)
      print 'Dimensions in Vector: ' + str(len(bctInputVec))
      print 'Boundary condition (' + INLET_FLOW_CAP + ') template file dimensions not compatible. Terminating.'
      sys.exit(-1)
  
  # Replace parameters in svpre file
  for i_param in parIntSvpre:
    i_param = int(i_param)

    param = svpreInputVec[i_param]
    # adding one line
    if (param[1] == 'STRING') or (param[1] == 'INT') or (param[1] == 'FLOAT'):
      templateSvpre = templateSvpre.replace('<params,' + str(i_param) + '>', param[0] )
    # adding multiple lines
    elif (param[1] == 'LIST'):
      data = ast.literal_eval(param[0])
      replace_string = ''
      listoflists = isinstance(data[0],list)
      if(listoflists):
        for i_data in xrange(0,len(data[0])):
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries][i_data]) + ' ' 
          replace_string = replace_string + '\n'
      else:
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries]) + '\n' 
      templateSvpre = templateSvpre.replace('<params,' + str(i_param) + '>', replace_string )
    else:
      sys.exit('ERROR: Invalid data type provided. Only INT, FLOAT, STRING, or LIST allowed\nNote: check capitalization, type is case senstive.')
  
  # Replace parameters in solver.inp file
  for i_param in parIntSolver:
    i_param = int(i_param)

    param = solverInputVec[i_param]
    # adding one line
    if (param[1] == 'STRING') or (param[1] == 'INT') or (param[1] == 'FLOAT'):
      templateSolver = templateSolver.replace('<params,' + str(i_param) + '>', param[0] )
    # adding multiple lines
    elif (param[1] == 'LIST'):
      data = ast.literal_eval(param[0])
      replace_string = ''
      listoflists = isinstance(data[0],list)
      if(listoflists):
        for i_data in xrange(0,len(data[0])):
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries][i_data]) + ' ' 
          replace_string = replace_string + '\n'
      else:
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries]) + '\n' 
      templateSolver = templateSolver.replace('<params,' + str(i_param) + '>', replace_string )
    else:
      sys.exit('ERROR: Invalid data type provided. Only INT, FLOAT, STRING, or LIST allowed\nNote: check capitalization, type is case senstive.')
  
  # Replace parameters in rcrt.dat file
  for i_param in parIntRcrt:
    i_param = int(i_param)

    param = rcrtInputVec[i_param]
    # adding one line
    if (param[1] == 'STRING') or (param[1] == 'INT') or (param[1] == 'FLOAT'):
      templateRcrt = templateRcrt.replace('<params,' + str(i_param) + '>', param[0] )
    # adding multiple lines
    elif (param[1] == 'LIST'):
      data = ast.literal_eval(param[0])
      replace_string = ''
      listoflists = isinstance(data[0],list)
      if(listoflists):
        for i_data in xrange(0,len(data[0])):
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries][i_data]) + ' ' 
          replace_string = replace_string + '\n'
      else:
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries]) + '\n' 
      templateRcrt = templateRcrt.replace('<params,' + str(i_param) + '>', replace_string )
    else:
      sys.exit('ERROR: Invalid data type provided. Only INT, FLOAT, STRING, or LIST allowed\nNote: check capitalization, type is case senstive.')

  # Replace parameters in INLET_FLOW_CAP file
  for i_param in parIntBct:
    i_param = int(i_param)

    param = bctInputVec[i_param]
    # adding one line
    if (param[1] == 'STRING') or (param[1] == 'INT') or (param[1] == 'FLOAT'):
      templateBct = templateBct.replace('<params,' + str(i_param) + '>', param[0] )
    # adding multiple lines
    elif (param[1] == 'LIST'):
      data = ast.literal_eval(param[0])
      replace_string = ''
      listoflists = isinstance(data[0],list)
      if(listoflists):
        for i_data in xrange(0,len(data[0])):
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries][i_data]) + ' ' 
          replace_string = replace_string + '\n'
      else:
          for i_entries in xrange(0,len(data)):
            replace_string = replace_string + str(data[i_entries]) + '\n' 
      templateBct = templateBct.replace('<params,' + str(i_param) + '>', replace_string )
    else:
      sys.exit('ERROR: Invalid data type provided. Only INT, FLOAT, STRING, or LIST allowed\nNote: check capitalization, type is case senstive.')
  

  # Save the realization files in the Inputs Folder
  os.chdir('../Inputs')

  if '.tmpl' in templateFileName:
    templateFileName = templateFileName.split('.')[0]
  svpreFileName = templateFileName  + '.svpre'
  svpreFile = open(svpreFileName,"w+")
  svpreFile.write(templateSvpre)
  svpreFile.close()

  solverFile = open('solver.inp',"w+")
  solverFile.write(templateSolver)
  solverFile.close()

  rcrtFile = open('rcrt.dat',"w+")
  rcrtFile.write(templateRcrt)
  rcrtFile.close()

  bctFile = open(INLET_FLOW_CAP,"w+")
  bctFile.write(templateBct)
  bctFile.close()

  # return to home directory
  os.chdir('../')



def Driver_3D_launch_sim(templateFileName,mesh):

  # Code executes from within /this_run folder

  # make new directory for this run
  now = datetime.datetime.utcnow()
  dirName = 'Run_' + now.strftime("%s")
  os.mkdir(dirName)

  # Move into new folder
  os.chdir(dirName)

  # Determine which mesh file
  meshDir = 'mesh' + mesh

  # Determine svpre file name
  if '.tmpl' in templateFileName:
    templateFileName = templateFileName.split('.')[0]
  svpreFileName = templateFileName  + '.svpre'

  # Move necessary files into this folder
  command_string = 'cp -r ../3D_Files/' +  meshDir + '/Solver_Files/* ../3D_Files/' +  meshDir + \
      '/mesh-complete ../Inputs/solver.inp ../Inputs/rcrt.dat ../Inputs/' + INLET_FLOW_CAP + ' ../Inputs/' + svpreFileName + \
      ' ../3D_Files/svpre ../3D_Files/svsolver ../3D_Files/svpost ../3D_Files/Postsolve* .'
  os.system(command_string)

  
  # Run the simulation
  command_string = "sbatch run.job.sh"
  os.system(command_string)
  print("Submitted 3D SLURM job")
  
  
  # wait for job to complete
  job_finished_check = glob.glob('../Outputs/3DOutputVectorNew.dat')
  while(len(job_finished_check) == 0):
    # keep waiting
    time.sleep(300)
    job_finished_check = glob.glob('../Outputs/3DOutputVectorNew.dat')
  print("3D SLURM job complete")

  

# 3D Driver
def Driver_3D(stochasticVector,templateFileName):
  
  #print('cd 3D_Files')
  os.chdir('3D_Files')

  # get default input vector values for .svpre file
  defaultSvpreVector = open('Defaults3DSvpre.dat').read().splitlines()
  defaultSvpreVector = [d.split(':') for d in defaultSvpreVector]

  # get default input vector values for solver.inp file
  defaultSolverVector = open('Defaults3DSolver.dat').read().splitlines()
  defaultSolverVector = [d.split(':') for d in defaultSolverVector]

  # get default input vector values for rcrt.dat file
  defaultRcrtVector = open('Defaults3DRcrt.dat').read().splitlines()
  defaultRcrtVector = [d.split(':') for d in defaultRcrtVector]

  # get default input vector values for INLET_FLOW_CAP file
  defaultBctVector = open('Defaults3DBct.dat').read().splitlines()
  defaultBctVector = [d.split(':') for d in defaultBctVector]


  # replace stochastic values
  mesh = stochasticVector[-1].strip().split(' ')[0]
  SvpreInputVec, SolverInputVec, RcrtInputVec, BctInputVec = Driver_3D_get_inputs(defaultSvpreVector,defaultSolverVector,defaultRcrtVector,defaultBctVector,stochasticVector)
  print("3D Inputs obtained")

  # fill template file
  Driver_3D_fill_template(SvpreInputVec,SolverInputVec,RcrtInputVec,BctInputVec,templateFileName,mesh)
  print("3D Templates filled")

  # Run the simulation and get simulation results
  Driver_3D_launch_sim(templateFileName,mesh)
  print("3D Simulation complete")

  
