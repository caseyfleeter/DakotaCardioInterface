import sys,os
import re
import numpy as np
import math
import argparse
import glob
import time
import datetime
import ast
import pandas as pd

# CONSTANTS
START_TIME = 5000
END_TIME = 6000
INCREMENT = 10
CYCLES = 6
CYCLE_STEPS = END_TIME/6/10

INLET_FLOW_FILE = "aorta_0_flow.dat"
INLET_PRESS_FILE = "aorta_0_pressure.dat"


def Driver_1D_get_inputs(defaultInputVec,stochasticVec):

  # Code executes while in 1D_Files folder

  # iterate through all stochastic values and replace the default value
  mergedInputVec = defaultInputVec
  numStoch = len(stochasticVec)
  stochasticVec = stochasticVec[0:numStoch-1]
  for i_stoch in stochasticVec:
    i_stoch = i_stoch.rsplit(' ')
    i_data = i_stoch[0]
    i_flag = i_stoch[1]

    # get index of parameter to be replaced
    flag_found = [i_flag in m for m in mergedInputVec]
    i_replace = [i for i, found in enumerate(flag_found) if found == True]

    # error check
    if len(i_replace) == 0:
      error_message = 'ERROR: invalid flag ' + i_flag + ' provided in stochastic vector'
      sys.exit(error_message)
    if len(i_replace) > 1:
      error_message = 'ERROR: flag ' + i_flag + ' appears multiple times in defaults vector'
      sys.exit(error_message)

    # replace default value
    if (mergedInputVec[i_replace[0]][1] == 'STRING') or (mergedInputVec[i_replace[0]][1] == 'INT') or (mergedInputVec[i_replace[0]][1] == 'FLOAT') or (mergedInputVec[i_replace[0]][1] == 'LIST'):
      mergedInputVec[i_replace[0]][0] = i_data

  # save the input vector
  os.chdir('../Inputs/')
  saveFile = open('1DInputVector.dat','w')
  saveFile.write('\n'.join([':'.join(m) for m in mergedInputVec]))
  saveFile.close()

  # return to 1D_Files folder
  os.chdir('../1D_Files/')

  return mergedInputVec



def Driver_1D_fill_template(inputVec,templateFileName,mesh):

  # Code executes while in 1D_Files folder

  # Generate Realization File: create 1D solver input file for this parameter vector
  # Open template file
  if '.tmpl' not in templateFileName:
    templateFileNameFull = templateFileName + mesh + '.tmpl'
  else:
    templateFileNameFull = templateFileName.split('.')[0]
    templateFileNameFull = templateFileNameFull + mesh + '.tmpl'
  template = open('Templates/' + templateFileNameFull,'r').read()
  parInt = re.findall(r'<params,(.*?)>', template)
  dims = len(np.unique(parInt))

  # Check if the dimensions are compatible with parameter vector
  if (dims != len(inputVec)):
      print('Dimensions in Template: ' + str(dims))
      print('Dimensions in Vector: ' + str(len(inputVec)))
      print('Template file dimensions not compatible. Terminating.')
      sys.exit(-1)
  
  # Replace parameters 
  for i_param in parInt:
    i_param = int(i_param)

    param = inputVec[i_param]
    # adding one line
    if (param[1] == 'STRING') or (param[1] == 'INT') or (param[1] == 'FLOAT'):
      template = template.replace('<params,' + str(i_param) + '>', param[0] )
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
      template = template.replace('<params,' + str(i_param) + '>', replace_string )
    else:
      sys.exit('ERROR: Invalid data type provided. Only INT, FLOAT, STRING, or LIST allowed\nNote: check capitalization, type is case senstive.')
    
  # Save the realization files in the Inputs Folder
  os.chdir('../Inputs')

  if '.tmpl' in templateFileName:
    templateFileName = templateFileName.split('.')[0]
  realizationFileName = templateFileName + mesh + '.in'
  realizationFile = open(realizationFileName,"w+")
  realizationFile.write(template)
  realizationFile.close()


  # return to home directory
  os.chdir('../')


def Driver_1D_launch_sim(fileName,mesh):

  # Code executes from within /this_run folder

  # make new directory for this run
  now = datetime.datetime.utcnow()
  dirName = 'Run_' + now.strftime("%s")
  os.mkdir(dirName)

  # Move into new folder
  os.chdir(dirName)

  if '.tmpl' in fileName:
    fileName = fileName.split('.')[0]
  realizationFileName = fileName + mesh + '.in'

  # Move necessary files into this folder
  command_string = 'cp ../1D_Files/OneDSolver ../Inputs/' + realizationFileName + ' .'
  os.system(command_string)

  # Run the simulation
  command_string = './OneDSolver ' + realizationFileName + ' > 1DModel.out 2> 1DModel.err &'
  os.system(command_string)
  print("Submitted 1D SLURM job")

  # Check if simulation has started
  outputFile = glob.glob('1DModel.o*')
  while(len(outputFile) == 0):
    # check for error
    errorFile = glob.glob('1DModel.e*')
    if len(errorFile) > 0:
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 1D solution terminated. See ' + errorFile[0] + ' for details')

    # keep waiting
    time.sleep(30)
    outputFile = glob.glob('1DModel.o*')
  
  # Check if simulation has completed
  completed_check = open(outputFile[0]).readlines()
  while(len(completed_check) == 0):
    # check for error
    errorFile = glob.glob('1DModel.e*')
    if len(errorFile) > 0:
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 1D solution terminated. See ' + errorFile[0] + ' for details')

    # keep waiting
    time.sleep(30)
    completed_check = open(outputFile[0]).readlines()

  # Get final line
  completed_check = open(outputFile[0]).readlines()[-1]
  while(completed_check != 'Completed!\n'):
    # check for error
    errorFile = glob.glob('1DModel.e*')
    if len(errorFile) > 0:
      if os.stat(errorFile[0]).st_size > 0:
        sys.exit('ERROR: 1D solution terminated. See ' + errorFile[0] + ' for details')

    # keep waiting
    time.sleep(30)
    completed_check = open(outputFile[0]).readlines()[-1]

  print("1D SLURM job complete")

  # return to home directory
  os.chdir('../')

  return dirName



def Driver_1D_get_results(fileName,runDirName,mesh):

  # Code executes from home directory 

  if '.tmpl' in fileName:
    fileName = fileName.split('.')[0]
  realizationFileName = fileName + mesh + '.in'

  outputVec = []
  # Extract quantities of interest
  flow_check = glob.glob('1D_Files/QoI-Flow.dat')
  if (flow_check):
    flow_files = open('1D_Files/QoI-Flow.dat').read().splitlines()
    for file in flow_files:
      # get flow results at the final FEM in this segment, except for inlet segment
      flows_result = pd.read_csv(runDirName + fileName + mesh + file, delimiter=' ',header=None)
      
      flows_result = flows_result.iloc[:,range(-CYCLE_STEPS-1,-1)]
      if(file != INLET_FLOW_FILE):  
        row = flows_result.tail(1)
      else:    
        row = flows_result.head(1)

      # get average flow at outlets over final cycle only
      outputVec.append(row.mean(axis=1).values[0])


  pressure_check = glob.glob('1D_Files/QoI-Pressure.dat')
  if(pressure_check):
    # save inlet pressures
    inlet_pressures = []
    # save pressure amplitudes
    amplitudes = []

    pressure_files = open('1D_Files/QoI-Pressure.dat').read().splitlines()
    for file in pressure_files:
      # get pressure results at the final FEM in this segment, except for inlet segment
      pressures_result = pd.read_csv(runDirName + fileName + mesh + file, delimiter=' ',header=None)
      
      pressures_result = pressures_result.iloc[:,range(-CYCLE_STEPS-1,-1)]
      if(file != INLET_PRESS_FILE):  
        row = pressures_result.tail(1)
      else:    
        row = pressures_result.head(1)
        inlet_pressures = row
      
      # get average pressure at outlets over final cycle only
      outputVec.append(row.mean(axis=1).values[0]/1333.2239)

      # save amplitude over final cycle only
      tmp_amp_press = row.max(axis=1).values[0] - row.min(axis=1).values[0]
      amplitudes.append(tmp_amp_press/1333.2239)

      # get pressure drop at outlets (NOT inlet) over final cycle only
      if(file != INLET_PRESS_FILE):
        tmp_drop_press = inlet_pressures.loc[inlet_pressures.index[0]] - row.loc[row.index[0]]
        # get average pressure drop at outlets over final cycle only
        outputVec.append(tmp_drop_press.mean()/1333.2239)
    

    # Save pressure amplitudes
    for amp in amplitudes:
      outputVec.append(amp)


  # Time-averaged quantities (pressure and WSS)
  ta_check = glob.glob('1D_Files/QoI-TA.dat')
  if(ta_check):

    wssThresholds = open('1D_Files/QoI-WSSThresholds.dat').read().strip('\n').split(',')
    modelLowTAWSS = np.zeros(len(wssThresholds))
    
    taNames = open('1D_Files/QoI-TA.dat').read().splitlines()
    for i_name in xrange(0,len(taNames)):
      taNames[i_name] = taNames[i_name].split(',')

    taStrips = open('1D_Files/QoI-TAStrips.dat').read().splitlines()
    for i_strip in xrange(0,len(taStrips)):
      taStrips[i_strip] = taStrips[i_strip].split(',')

    taBranches = open('1D_Files/QoI-TABranches.dat').read().splitlines()
    for i_branch in xrange(0,len(taBranches)):
      taBranches[i_branch] = taBranches[i_branch].split(',')

    # Now, get TABP and TAWSS quantities across entire model
    modelMaxTABP = -sys.float_info.max
    modelMinTABP = sys.float_info.max
    modelMeanTABP = 0.0

    branchesMaxTABP = [-sys.float_info.max for i in xrange(0,len(taNames))]
    branchesMinTABP = [sys.float_info.max for i in xrange(0,len(taNames))]
    branchesMeanTABP = [0.0 for i in xrange(0,len(taNames))]

    stripsMaxTABP = [-sys.float_info.max for i in xrange(0,len(taStrips))]
    stripsMinTABP = [sys.float_info.max for i in xrange(0,len(taStrips))]
    stripsMeanTABP = [0.0 for i in xrange(0,len(taStrips))]
    
    modelMaxTAWSS = -sys.float_info.max
    modelMinTAWSS = sys.float_info.max
    modelMeanTAWSS = 0.0
   
    branchesMaxTAWSS = [-sys.float_info.max for i in xrange(0,len(taNames))]
    branchesMinTAWSS = [sys.float_info.max for i in xrange(0,len(taNames))]
    branchesMeanTAWSS = [0.0 for i in xrange(0,len(taNames))]

    stripsMaxTAWSS = [-sys.float_info.max for i in xrange(0,len(taStrips))]
    stripsMinTAWSS = [sys.float_info.max for i in xrange(0,len(taStrips))]
    stripsMeanTAWSS = [0.0 for i in xrange(0,len(taStrips))]

    modelTotalArea = 0.0
    branchesTotalArea = [0.0 for i in xrange(0,len(taNames))]
    stripsTotalArea = [0.0 for i in xrange(0,len(taStrips))]

    
    # Loop over all the result files for MODEL and BRANCH QoIs
    for i_branch in xrange(0,len(taNames)):

      # Getting starting and ending segments for this branch
      branch_start = taBranches[i_branch][1]
      branch_end = taBranches[i_branch][2]

      # Loop over segments of current branch
      for i_num in xrange(0,int(taNames[i_branch][1])+1):

        # read model input file for extracting geometry
        model_file = open('Inputs/' + realizationFileName,'r')

        # Get length and number of FEM of current segment
        findString = 'SEGMENT ' + str(taNames[i_branch][0]) + '_' + str(i_num)
        match = False
        for line in model_file:
          if re.search(findString, line):
            match = True
            
            result = line.split(' ')
            l = float(result[3])
            numFE = int(result[4])
            break

        if match == False:
          print('ERROR: Match for branch name not found in TAModelExtractor.')
          print('String = ' + findString)
          print('Terminating.')
          model_file.close()
          sys.exit(0)

        # Get areas of all FEM of current segment
        # Get current Area Result File
        areaFileName = str(fileName) + str(mesh) + str(taNames[i_branch][0]) + '_' + str(i_num) + '_area.dat'
        
        # Open the file with pandas
        currSegmentArea = pd.read_csv(runDirName + areaFileName, delimiter=' ',header=None)
        currSegmentArea = currSegmentArea.iloc[:,range(-CYCLE_STEPS-1,-1)]

        # Time average these values at all FEM of the segment
        currSegmentArea = currSegmentArea.mean(axis=1)
        
        # Split into start areas and end areas of each FEM
        area1 = currSegmentArea[:-1]
        area2 = currSegmentArea[1:]

        # Get length of reach FEM
        length = l / numFE

        # Get surface error of current segment (Isosceles trapezoid)
        currSurfaceArea = [0.25 * np.sqrt( ((a1+a2)**2) * np.abs(a1 - a2 + 2*length) * np.abs(a2 - a1 + 2*length) ) for a1,a2 in zip(area1, area2)]

        # Update total surface area
        for csa in currSurfaceArea:
          modelTotalArea = modelTotalArea + csa

        # Update total branch area if applicable
        if(i_num >= int(branch_start) and i_num <= int(branch_end)):
          for csa in currSurfaceArea:
            branchesTotalArea[i_branch] = branchesTotalArea[i_branch] + csa

        # Get current Pressure Result File
        tabpFileName = str(fileName) + str(mesh) + str(taNames[i_branch][0]) + '_' + str(i_num) + '_pressure.dat'
        # Get current WSS Result File
        tawssFileName = str(fileName) + str(mesh) + str(taNames[i_branch][0]) + '_' + str(i_num) + '_wss.dat'
        
        # Open the file with pandas
        currSegmentTABP_init = pd.read_csv(runDirName + tabpFileName, delimiter=' ',header=None)
        currSegmentTABP_init = currSegmentTABP_init.iloc[:,range(-CYCLE_STEPS-1,-1)]
        
        currSegmentTAWSS_init = pd.read_csv(runDirName + tawssFileName, delimiter=' ',header=None)      
        currSegmentTAWSS_init = currSegmentTAWSS_init.iloc[:,range(-CYCLE_STEPS-1,-1)]

        # Time average these values at all FEM of the segment
        currSegmentTABP_init = currSegmentTABP_init.mean(axis=1)
        currSegmentTAWSS_init = currSegmentTAWSS_init.mean(axis=1)     

        # Get values on each FEM of section
        currSegmentTABP = np.empty(numFE)
        currSegmentTAWSS = np.empty(numFE)
        for i_fem in xrange(0,numFE):
          currSegmentTABP[i_fem] = (currSegmentTABP_init[i_fem] + currSegmentTABP_init[i_fem+1])/2
          currSegmentTAWSS[i_fem] = (currSegmentTAWSS_init[i_fem] + currSegmentTAWSS_init[i_fem+1])/2

        # Update MODEL mean TABP and TAWSS and thresholded TAWSS
        for bp,wss,csa in zip(currSegmentTABP,currSegmentTAWSS,currSurfaceArea):
          modelMeanTABP = modelMeanTABP + bp*csa
          modelMeanTAWSS = modelMeanTAWSS + wss*csa
          
          # update area of low TAWSS regions
          for i_th in xrange(0,len(wssThresholds)):
            if wss < float(wssThresholds[i_th]):
              modelLowTAWSS[i_th] = modelLowTAWSS[i_th] + csa

        # Update MODEL min/max TABP and TAWSS
        maxSegmentTABP = currSegmentTABP.max()
        minSegmentTABP = currSegmentTABP.min()
        # Store if needed
        if(maxSegmentTABP > modelMaxTABP):
          modelMaxTABP = maxSegmentTABP
        if(minSegmentTABP < modelMinTABP):
          modelMinTABP = minSegmentTABP

        maxSegmentTAWSS = currSegmentTAWSS.max()
        minSegmentTAWSS = currSegmentTAWSS.min()
        # Store if needed
        if(maxSegmentTAWSS > modelMaxTAWSS):
          modelMaxTAWSS = maxSegmentTAWSS
        if(minSegmentTAWSS < modelMinTAWSS):
          modelMinTAWSS = minSegmentTAWSS


        # Update BRANCH QoIs if applicable
        if(i_num >= int(branch_start) and i_num <= int(branch_end)):

          # Update BRANCH mean TABP and TAWSS and thresholded TAWSS
          for bp,wss,csa in zip(currSegmentTABP,currSegmentTAWSS,currSurfaceArea):
            branchesMeanTABP[i_branch] = branchesMeanTABP[i_branch] + bp*csa
            branchesMeanTAWSS[i_branch] = branchesMeanTAWSS[i_branch] + wss*csa


          # Update BRANCH min/max TABP and TAWSS
          maxSegmentTABP = currSegmentTABP.max()
          minSegmentTABP = currSegmentTABP.min()
          # Store if needed
          if(maxSegmentTABP > branchesMaxTABP[i_branch]):
            branchesMaxTABP[i_branch] = maxSegmentTABP
          if(minSegmentTABP < branchesMinTABP[i_branch]):
            branchesMinTABP[i_branch] = minSegmentTABP

          maxSegmentTAWSS = currSegmentTAWSS.max()
          minSegmentTAWSS = currSegmentTAWSS.min()
          # Store if needed
          if(maxSegmentTAWSS > branchesMaxTAWSS[i_branch]):
            branchesMaxTAWSS[i_branch] = maxSegmentTAWSS
          if(minSegmentTAWSS < branchesMinTAWSS[i_branch]):
            branchesMinTAWSS[i_branch] = minSegmentTAWSS


    # spatially average MODEL mean TABP and TAWSS
    modelMeanTABP = modelMeanTABP / modelTotalArea
    modelMeanTAWSS = modelMeanTAWSS / modelTotalArea

    # spatially average BRANCHES mean TABP and TAWSS
    branchesMeanTABP = [bp / a for bp,a in zip(branchesMeanTABP, branchesTotalArea)]

    branchesMeanTAWSS = [wss / a for wss,a in zip(branchesMeanTAWSS, branchesTotalArea)]

    # make lowWSS a percentage of total area
    modelLowTAWSS = [l / modelTotalArea for l in modelLowTAWSS]


    # Loop over all the result files for STRIP QoIs
    for i_strip in xrange(0,len(taStrips)):

      # Get branch name and all segments for this strip
      strip_name = taStrips[i_strip][0]
      strip_segments = taStrips[i_strip][1:]

      # Loop over segments of current strip
      for i_num in strip_segments:

        # read model input file for extracting geometry
        model_file = open('Inputs/' + realizationFileName,'r')

        # Get surface area of current segment
        findString = 'SEGMENT ' + str(strip_name) + '_' + str(i_num)
        match = False
        for line in model_file:
          if re.search(findString, line):
            match = True
            
            result = line.split(' ')
            l = float(result[3])
            numFE = int(result[4])
            break

        if match == False:
          print('ERROR: Match for branch name not found in TAModelExtractor.')
          print('String = ' + findString)
          print('Terminating.')
          model_file.close()
          sys.exit(0)

        # Get areas of all FEM of current segment
        # Get current Area Result File
        areaFileName = str(fileName) + str(mesh) + str(strip_name) + '_' + str(i_num) + '_area.dat'
        
        # Open the file with pandas
        currSegmentArea = pd.read_csv(runDirName + areaFileName, delimiter=' ',header=None)
        currSegmentArea = currSegmentArea.iloc[:,range(-CYCLE_STEPS-1,-1)]

        # Time average these values at all FEM of the segment
        currSegmentArea = currSegmentArea.mean(axis=1)
        
        # Split into start areas and end areas of each FEM
        area1 = currSegmentArea[:-1]
        area2 = currSegmentArea[1:]

        # Get length of reach FEM
        length = l / numFE

        # Get surface error of current segment (Isosceles trapezoid)
        currSurfaceArea = [0.25 * np.sqrt( ((a1+a2)**2) * (a1 - a2 + 2*length) * (a2 - a1 + 2*length) ) for a1,a2 in zip(area1, area2)]
        
        # Update total surface area
        for csa in currSurfaceArea:
          stripsTotalArea[i_strip] = stripsTotalArea[i_strip] + csa


        # Get current BP Result File
        tabpFileName = str(fileName) + str(mesh) + str(strip_name) + '_' + str(i_num) + '_pressure.dat'
        # Get current WSS Result File
        tawssFileName = str(fileName) + str(mesh) + str(strip_name) + '_' + str(i_num) + '_wss.dat'
        
        # Open the file with pandas
        currSegmentTABP_init = pd.read_csv(runDirName + tabpFileName, delimiter=' ',header=None)
        currSegmentTABP_init = currSegmentTABP_init.iloc[:,range(-CYCLE_STEPS-1,-1)]
        
        currSegmentTAWSS_init = pd.read_csv(runDirName + tawssFileName, delimiter=' ',header=None)      
        currSegmentTAWSS_init = currSegmentTAWSS_init.iloc[:,range(-CYCLE_STEPS-1,-1)]

        # Time average these values at all FEM of the segment
        currSegmentTABP_init = currSegmentTABP_init.mean(axis=1)
        currSegmentTAWSS_init = currSegmentTAWSS_init.mean(axis=1) 

        # Get values on each FEM of section
        currSegmentTABP = np.empty(numFE)
        currSegmentTAWSS = np.empty(numFE)
        for i_fem in xrange(0,numFE):
          currSegmentTABP[i_fem] = (currSegmentTABP_init[i_fem] + currSegmentTABP_init[i_fem+1])/2
          currSegmentTAWSS[i_fem] = (currSegmentTAWSS_init[i_fem] + currSegmentTAWSS_init[i_fem+1])/2


        # Update STRIP QoIs 
        # Update STRIP mean TABP and TAWSS and thresholded TAWSS
        for bp,wss,csa in zip(currSegmentTABP,currSegmentTAWSS,currSurfaceArea):
          stripsMeanTABP[i_strip] = stripsMeanTABP[i_strip] + bp*csa
          stripsMeanTAWSS[i_strip] = stripsMeanTAWSS[i_strip] + wss*csa
          

        # Update STRIP min/max TABP and TAWSS
        maxSegmentTABP = currSegmentTABP.max()
        minSegmentTABP = currSegmentTABP.min()
        # Store if needed
        if(maxSegmentTABP > stripsMaxTABP[i_strip]):
          stripsMaxTABP[i_strip] = maxSegmentTABP
        if(minSegmentTABP < stripsMinTABP[i_strip]):
          stripsMinTABP[i_strip] = minSegmentTABP

        maxSegmentTAWSS = currSegmentTAWSS.max()
        minSegmentTAWSS = currSegmentTAWSS.min()
        # Store if needed
        if(maxSegmentTAWSS > stripsMaxTAWSS[i_strip]):
          stripsMaxTAWSS[i_strip] = maxSegmentTAWSS
        if(minSegmentTAWSS < stripsMinTAWSS[i_strip]):
          stripsMinTAWSS[i_strip] = minSegmentTAWSS

    # spatially average STRIPS mean TABP and TAWSS
    stripsMeanTABP = [bp / a for bp,a in zip(stripsMeanTABP, stripsTotalArea)]
    stripsMeanTAWSS = [wss / a for wss,a in zip(stripsMeanTAWSS, stripsTotalArea)]

    # Save all in the outputVec

    # TABP results first: branches, model, strips (in order: mean, min, max)
    for i_branch in xrange(0,len(taNames)):
      outputVec.append(branchesMeanTABP[i_branch]/1333.2933)
      outputVec.append(branchesMinTABP[i_branch]/1333.2933)
      outputVec.append(branchesMaxTABP[i_branch]/1333.2933)

    outputVec.append(modelMeanTABP/1333.2933)
    outputVec.append(modelMinTABP/1333.2933)
    outputVec.append(modelMaxTABP/1333.2933)

    for i_strip in xrange(0,len(taStrips)):
      outputVec.append(stripsMeanTABP[i_strip]/1333.2933)
      outputVec.append(stripsMinTABP[i_strip]/1333.2933)
      outputVec.append(stripsMaxTABP[i_strip]/1333.2933)

    # TAWSS results second: branches, model, strips (in order: mean, min, max, thresholds (when applicable))
    for i_branch in xrange(0,len(taNames)):
      outputVec.append(branchesMeanTAWSS[i_branch])
      outputVec.append(branchesMinTAWSS[i_branch])
      outputVec.append(branchesMaxTAWSS[i_branch])

    outputVec.append(modelMeanTAWSS)
    outputVec.append(modelMinTAWSS)
    outputVec.append(modelMaxTAWSS)

    for i_wss in xrange(0,len(modelLowTAWSS)):
      outputVec.append(modelLowTAWSS[i_wss])

    for i_strip in xrange(0,len(taStrips)):
      outputVec.append(stripsMeanTAWSS[i_strip])
      outputVec.append(stripsMinTAWSS[i_strip])
      outputVec.append(stripsMaxTAWSS[i_strip])


  # Return QoI vector
  return outputVec




# 1D Driver
def Driver_1D(stochasticVector,templateFileName):
  
  os.chdir('1D_Files')

  # get default input vector values
  defaultInputVector = open('Defaults1D.dat').read().splitlines()
  defaultInputVector = [d.split(':') for d in defaultInputVector]


  # replace stochastic values
  mesh = stochasticVector[-1].strip().split(' ')[0]
  inputVec = Driver_1D_get_inputs(defaultInputVector,stochasticVector)
  sys.stdout.write("1D Inputs obtained\n")
  sys.stdout.flush()


  # fill template file
  Driver_1D_fill_template(inputVec,templateFileName,mesh)
  sys.stdout.write("1D Templates filled\n")
  sys.stdout.flush()

  # Run the simulation
  runDirName = Driver_1D_launch_sim(templateFileName,mesh)
  sys.stdout.write("1D Simulation complete\n")
  sys.stdout.flush()

  # get simulation results
  dirName = runDirName + '/'
  outputVec = Driver_1D_get_results(templateFileName,dirName,mesh)
  
  # Write output to file
  os.chdir('Outputs/')
  np.savetxt('0DOutputVector.dat',outputVec,delimiter=',',fmt='%s')

  # Return to home directory
  os.chdir('../')

  sys.stdout.write("1D Results obtained\n")
  sys.stdout.flush()


