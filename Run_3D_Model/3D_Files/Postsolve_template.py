# --- VTK-PYTHON SCRIPT FOR READING VTP, INTEGRATING VARIABLES OF INTEREST USING TRAPEZOIDAL RULE
# --- CALCULATES FLOW AND PRESSURE AT OUTLETS OF THE MODEL
# --- NOTE: TO BE RUN FROM DIRECTORY CONTAINING all_results.vtu FILE AND mesh-complete/ DIRECTORY
# --- NOTE: THIS SCRIPT ASSUMES A TRIANGULAR SURFACE MESH

import sys
import os
import vtk
import numpy as np
import glob

# SPECIFIY TIME POINTS TO DETERMINE QoIs
START_TIME = 500
END_TIME = 600
INCREMENT = 10

# NAMES OF CAP AND WALL GEOMETRY FILES OF MODEL, WITHOUT CAP_/WALL_ PREFIX
input_filenames = ['aorta_2','celiac_branch','right_iliac']
inlet_file = 'aorta'

# POINT AND NORMAL DEFINING PLANES FOR CLIPPING FULL MODEL TO EACH BRANCH AND STRIP OF INTEREST
branch_clip_planes = [  [[[-0.708118482547202, -0.903954731420905, -7.55794848970956], [-3.68450951964806, 0.874679925552799, -5.86868534441482], [-6.61416460962624, 2.47637539848582, -5.83113193746924]],
                        [[-0.905949297123049, -0.0803344775715136, -0.415694891423464], [-0.455898769227772, 0.871951474909039, 0.178485118764996], [0.409480528765804, 0.166464849522705, -0.897003428329603]]
                        ],
                        [[[-1.94298892480308, 1.71454311555629, 5.65017670683522], [-1.30287036931842, 5.83344677863165, 4.55818706302802],[-0.819141864496289, 3.97939951146734, 4.22053052732145]],
                         [[-0.517817224783634, 0.851546634307298, 0.0820588283933473], [0.979729347213665, 0.191351715685473, -0.0592868207320274], [-0.00489053848000058, -0.107074250159097, 0.994238999228175]]
                        ],
                        [[[-1.62704205614731, 5.47782376675798, -5.60856532837656], [-6.69572377298467, -0.341179982123416, 3.54104649308491],[1.97967518206377, 4.01611977721367, 4.11511862973538],[-0.440964081252454, 4.15578703295313, 2.83330182336447]],
                         [[-0.28955051714008, 0.902939993223789, 0.31758442446212], [0.459716134775713, -0.47936478691392, -0.747576401777003], [-0.435062826240822, 0.883935498956326, 0.171401198677779],[0.198293801576317, -0.812399924239667, -0.548348366781369]]
                        ]
                     ]

strip_clip_planes = [ [[[-1.66715865760017, -0.988124292679531, 6.13460968029932], [-1.88242297393781, 0.0303006465428083, 8.49511915378237]],
                       [[-0.0653674041139646, -0.245497018496443, 0.967190941018762], [-0.0153293827308922, 0.300741602714008, -0.953582454967634]]
                      ],
                      [[[-1.58476920351838, 0.0583473352559409, 1.46456681294995], [-2.24495184966637, 3.18965162505489, 0.508682254830264],[0.279690756527583, 2.08205290731091, 1.32131888619727]],
                       [[0.0196896475531627, 0.382818817176909, -0.923613593985333], [-0.0223457613321743, -0.331978819560708, 0.943022126099682], [-0.514366991413212, -0.738593498846872, -0.435782332828738]]
                      ]
                    ]

# LOW WSS THRESHHOLD VALUES
thresholds = [1.0, 5.0, 10.0, 20.0]

# PATH OF OUTPUT FILE FOR QoIs
output_filename_flow = 'all_results_flows.dat'
output_collection_flow = []    # This holds the flows averaged on each outlet cap at each time step

output_filename_pressure = 'all_results_pressures.dat'
output_collection_pressure = [] # This holds the pressures averaged on each outlet cap and 
                                # pressure drop from inlet to outlet at each time step
output_filename_pressure_2 = 'all_results_pressures_2.dat'                                
output_collection_pressure_2 = [] # This holds the pressure amplitude at each outlet cap (over entire time cycle)                    

branch_vols = [] # This holds the volumes of each branch, used to volumetrically average
output_filename_tabp = 'all_results_tabp.dat'
output_collection_tabp = [] # This holds the mean, min and max TABP volumetrically averaged in various strips,
                            # the mean, min and max TABP volumetrically averaged in each branch,
                            # and the mean, min and max TABP volumetrically averaged in the whole model (over entire time cycle)  

branch_areas = [] # This holds the surface areas of each branch, used to spatially average
output_filename_tawss = 'all_results_tawss.dat'
output_collection_tawss = [] # This holds the mean, min and max TAWSS spatially averaged in each branch,
                             # the mean, min and max TAWSS spatially averaged in the whole model,
                             # the mean, min and max TAWSS spatially averaged in various strips, 
                             # and the % area in the model with TAWSS less than various thresholds (over entire time cycle)


if __name__ == "__main__":

  # how many timesteps for post-processing
  timestep_count = len(xrange(START_TIME,END_TIME+INCREMENT,INCREMENT)) 

  #############################################################
  #  Read in the .vtp file containing quantities of interest  #
  #############################################################
  all_results_vtp_reader = vtk.vtkXMLPolyDataReader()             # Create vtk instance
  all_results_vtp_reader.SetFileName('all_results.vtp')           # Open file
  all_results_vtp_reader.Update()
  all_results_vtp_model = vtk.vtkPolyData() 
  all_results_vtp_model = all_results_vtp_reader.GetOutput()          # Read file into new variable for manipulation
  all_results_vtp_numPts = all_results_vtp_model.GetNumberOfPoints()  # Determine number of points in the mesh of the entire model
  all_results_vtp_IDs = all_results_vtp_model.GetPointData().GetArray('GlobalNodeID')   # Extract node IDs of full model solution
  all_results_vtp_numCells = all_results_vtp_model.GetNumberOfCells() 


  #############################################################
  #  Read in the .vtu file containing quantities of interest  #
  #############################################################
  all_results_vtu_reader = vtk.vtkXMLUnstructuredGridReader()             # Create vtk instance
  all_results_vtu_reader.SetFileName('all_results.vtu')           # Open file
  all_results_vtu_reader.Update()
  all_results_vtu_model = all_results_vtu_reader.GetOutput()          # Read file into new variable for manipulation
  all_results_vtu_numPts = all_results_vtu_model.GetNumberOfPoints()  # Determine number of points in the mesh of the entire model
  all_results_vtu_IDs = all_results_vtu_model.GetPointData().GetArray('GlobalNodeID')   # Extract node IDs of full model solution
  all_results_vtu_numCells = all_results_vtu_model.GetNumberOfCells() 


  #####################################################
  #  Load data for caps, walls and interior of model  #
  #####################################################

  # For velocity arrays from the all_results.vtp
  velocity_vectors = []

  # For pressure arrays from the all_results.vtp
  pressure_vectors = []

  # For pressure arrays from the all_results.vtu
  pressure_vectors_vtu = []

  # keep track of how many timesteps in solution
  t_count = 0
  for i_time in xrange(START_TIME, END_TIME+INCREMENT, INCREMENT):
    velocity_vectors.append(vtk.vtkDoubleArray())
    pressure_vectors.append(vtk.vtkDoubleArray())
    pressure_vectors_vtu.append(vtk.vtkDoubleArray())
    
    if i_time < 10:
      velocity_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('velocity_' + '0000' + str(i_time))
      pressure_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('pressure_' + '0000' + str(i_time))
      pressure_vectors_vtu[t_count] = all_results_vtu_model.GetPointData().GetArray('pressure_' + '0000' + str(i_time))
    elif i_time < 100:
      velocity_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('velocity_' + '000' + str(i_time))
      pressure_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('pressure_' + '000' + str(i_time))
      pressure_vectors_vtu[t_count] = all_results_vtu_model.GetPointData().GetArray('pressure_' + '000' + str(i_time))
    elif i_time < 1000:
      velocity_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('velocity_' + '00' + str(i_time))
      pressure_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('pressure_' + '00' + str(i_time))
      pressure_vectors_vtu[t_count] = all_results_vtu_model.GetPointData().GetArray('pressure_' + '00' + str(i_time))
    elif i_time < 10000:
      velocity_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('velocity_' + '0' + str(i_time))
      pressure_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('pressure_' + '0' + str(i_time))
      pressure_vectors_vtu[t_count] = all_results_vtu_model.GetPointData().GetArray('pressure_' + '0' + str(i_time))
    else:
      velocity_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('velocity_' + str(i_time)) 
      pressure_vectors[t_count] = all_results_vtp_model.GetPointData().GetArray('pressure_' + str(i_time)) 
      pressure_vectors_vtu[t_count] = all_results_vtu_model.GetPointData().GetArray('pressure_' + str(i_time))
    
    t_count += 1

  # For WSS arrays from the all_results.vtp
  wss_vectors = vtk.vtkDoubleArray()
  wss_vectors = all_results_vtp_model.GetPointData().GetArray('vTAWSS')


  #############################
  #  Calculate TABP in model  #
  #############################

  tabp_vectors = pressure_vectors_vtu[0]
  
  # sum pressures at all mesh points over all times
  for i_pt in xrange(0, all_results_vtu_numPts):
    for i_time in xrange(1, t_count):
      update_val = [t + u for t,u in zip(tabp_vectors.GetTuple(i_pt),pressure_vectors_vtu[i_time].GetTuple(i_pt))]
      tabp_vectors.SetTuple(i_pt,update_val)

    # average pressure at each mesh point in time
    update_val = [float(t/t_count) for t in tabp_vectors.GetTuple(i_pt)]
    tabp_vectors.SetTuple(i_pt,update_val)


  ######################################################
  #  Extract QoIs from each wall and cap of the model  #
  ######################################################
  num_file = -1;
  for file in input_filenames:
    print('Post-processing ' + file)

    ##########################################
    #  Find node IDs for cap of the model    #
    ##########################################

    # Load in the mesh file for the cap you want flows and pressures. This can be found in the mesh-surfaces folder.
    os.chdir('mesh-complete/mesh-surfaces')

    # Read geometry (mesh) information from this cap
    cap_reader = vtk.vtkXMLPolyDataReader()          # Create vtk instance
    cap_reader.SetFileName('cap_' + file + '.vtp')   # Open file
    cap_reader.Update()
    cap_model = vtk.vtkPolyData()
    cap_model = cap_reader.GetOutput()            # Read file into new variable for manipulation
    cap_numPts = cap_model.GetNumberOfPoints()           # Determine number of points in the mesh at this outlet
    cap_IDs = cap_model.GetPointData().GetArray("GlobalNodeID")  # Extract node IDs to match with full model solution
    
    # Compute the normals to each mesh cell on the cap
    normalGenerator = vtk.vtkPolyDataNormals()          # Create vtk instance
    normalGenerator.SetInputData(cap_model)          # Open file
    normalGenerator.ComputePointNormalsOff()
    normalGenerator.ComputeCellNormalsOn()              # normals to each mesh cell
    normalGenerator.Update()
    normals_test = normalGenerator.GetOutput()          # Read file into new variable for manipulation
    cap_normal = normals_test.GetCellData().GetArray("Normals").GetTuple3(1) # one normal for the entire cap, assumed flat


    ##########################################
    #  Find node IDs for wall of the model   #
    ##########################################

    # Read geometry (mesh) information from this branch wall
    wall_check = glob.glob('wall_' + file + '.vtp')
    if(len(wall_check) != 0):
      wall_reader = vtk.vtkXMLPolyDataReader()          # Create vtk instance
      wall_reader.SetFileName('wall_' + file + '.vtp')   # Open file
      wall_reader.Update()
      wall_model = vtk.vtkPolyData()
      wall_model = wall_reader.GetOutput()            # Read file into new variable for manipulation
      wall_numPts = wall_model.GetNumberOfPoints()           # Determine number of points in the mesh at this outlet
      wall_IDs = wall_model.GetPointData().GetArray("GlobalNodeID")  # Extract node IDs to match with full model solution


    ##############################################
    #  Find node IDs for branches of the model   #
    ##############################################

    # check if this is a valid branch
    branch_check = glob.glob('wall_' + file + '.vtp')

    # return to the home directory
    os.chdir('../..')

    # get branch if valid
    if(len(branch_check) != 0):
      # increment file count
      num_file = num_file + 1

      # Get geometry (mesh) information for this branch wall
      branch_model = all_results_vtu_reader.GetOutput()
      
      # Iterate over all planes needed for clipping model
      for i_plane in xrange(0,len(branch_clip_planes[num_file][0])):  
        curr_point = branch_clip_planes[num_file][0][i_plane]
        curr_norm = branch_clip_planes[num_file][1][i_plane]
        
        clipPlane = vtk.vtkPlane()
        clipPlane.SetOrigin(curr_point)
        clipPlane.SetNormal(curr_norm)

        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(branch_model)
        clipper.SetClipFunction(clipPlane) 
        clipper.Update()

        triangulator = vtk.vtkDataSetTriangleFilter()
        triangulator.SetInputData(clipper.GetOutput())
        triangulator.Update()

        branch_model = triangulator.GetOutput()

      # Get info for branch model
      branch_numPts = branch_model.GetNumberOfPoints()           # Determine number of points in the mesh at this outlet
      branch_IDs = branch_model.GetPointData().GetArray("GlobalNodeID")  # Extract node IDs to match with full model solution

      # Write out the branch to ensure they all look fine
      # branch_writer = vtk.vtkXMLUnstructuredGridWriter()
      # branch_writer.SetInputData(branch_model)
      # branch_writer.SetFileName('branch_' + file + '.vtu')
      # branch_writer.Write()
    


    ################################################################################################
    #  Find IDs of corresponding nodes (for caps/walls/branches) in all_results.vt(p/u) for model  #
    ################################################################################################

    # Find the nodes on all_results.vtp that correspond to the cap of interest
    cap_nodes = []
    # iterate through all nodes in model
    for i_node in xrange(0, cap_numPts):
      cap_ID = cap_IDs.GetTuple1(i_node)

      # iterate through all nodes in model
      for i_full in xrange(0, all_results_vtp_numPts):
        full_ID = all_results_vtp_IDs.GetTuple1(i_full)

        if(full_ID == cap_ID):
          cap_nodes.append(i_full)
          break
          
    # Just to make sure we found all the cap nodes in all_results
    assert(len(cap_nodes) == cap_numPts)


    # Find the nodes on all_results.vtp that correspond to the wall of interest
    wall_nodes = []
    # iterate through all nodes in model
    for i_node in xrange(0, wall_numPts):
      wall_ID = wall_IDs.GetTuple1(i_node)

      # iterate through all nodes in model
      for i_full in xrange(0, all_results_vtp_numPts):
        full_ID = all_results_vtp_IDs.GetTuple1(i_full)

        if(full_ID == wall_ID):
          wall_nodes.append(i_full)
          break
          
    # Just to make sure we found all the wall nodes in all_results
    assert(len(wall_nodes) == wall_numPts)


    # Find the nodes on all_results.vtp that correspond to the branch of interest
    branch_nodes = []
    # iterate through all nodes in model
    for i_node in xrange(0, branch_numPts):
      branch_ID = branch_IDs.GetTuple1(i_node)

      # iterate through all nodes in model
      for i_full in xrange(0, all_results_vtu_numPts):
        full_ID = all_results_vtu_IDs.GetTuple1(i_full)

        if(full_ID == branch_ID):
          branch_nodes.append(i_full)
          break
          
    # Just to make sure we found all the branch nodes in all_results
    assert(len(branch_nodes) == branch_numPts)
    

    ############################################
    #  Extract normal velocities for each cap  #
    ############################################  
    
    # Create a matrix to hold all the normal velocities for every node at every time on the outlet face
    normal_velocities = np.zeros((cap_numPts, t_count))
    
    # Fill out the normal_velocities matrix by first looping over the nodes in the outlet,
    # computing the normal velocity, then storing it for each time point
    for i_node in xrange(0, len(cap_nodes)):
      
      node_ind = cap_nodes[i_node]
      for i_time in xrange(0, t_count):
          
        nodal_velocities = velocity_vectors[i_time]
        x_vel = nodal_velocities.GetTuple3(node_ind)[0]
        y_vel = nodal_velocities.GetTuple3(node_ind)[1]
        z_vel = nodal_velocities.GetTuple3(node_ind)[2]
          
        norm_vel_temp = cap_normal[0]*x_vel + \
                        cap_normal[1]*y_vel + \
                        cap_normal[2]*z_vel
                          
        normal_velocities[i_node][i_time] = norm_vel_temp
        

    ######################################################
    #  Calculate velocity and pressure QoIs at each cap  #
    ######################################################  

    # Now that we know the normal velocities at every node and time, integrate
    # them over the surface of the outlet face to get the flow on this face at every time
    temp_mean_flow = np.zeros(t_count)
    
    # Integrate pressures over the surface of the outlet face to get the pressure on this face at each time
    temp_mean_press = np.zeros(t_count)
    temp_amp_press  = 0.0  # also the amplitude at the cap (over all time steps)
    if(file != inlet_file):
      temp_drop_press = np.zeros(t_count) # if not the inlet cap, then also get the pressure drop inlet to outlet

    for i_time in xrange(0, t_count):
    
      # Compute the integral using trapezoidal rule
      total_area = 0.0
      
      # store velocity information for the entire outlet face at this time step
      curr_flow = 0.0 

      # store pressure information for the entire outlet face at this time step
      curr_press = 0.0

      # iterate over all mesh cells on outlet face
      for i_cell in xrange(0, cap_model.GetNumberOfCells()):
        
        # extract information about cell vertices
        temp_cell = cap_model.GetCell(i_cell)
        pts_cell = temp_cell.GetPointIds()
        cell_pts = temp_cell.GetPoints()
        p0 = cell_pts.GetPoint(0)
        p1 = cell_pts.GetPoint(1)
        p2 = cell_pts.GetPoint(2)
          
        # compute area of mesh cell (triangular mesh assumed)  
        local_area = temp_cell.TriangleArea(p0, p1, p2)
        total_area = total_area + local_area
        
        local_temp_flow = 0.0
        local_temp_press = 0.0
        # add contributions from each vertex of cell
        for ipt in xrange(0, pts_cell.GetNumberOfIds()):
          
          iid = pts_cell.GetId(ipt)     # get node number of this point on cap
          node_ind = cap_nodes[iid]     # get node number on full model
     
          norm_vel_temp = normal_velocities[iid][i_time]  # get normal velocity of this node
          # add normal velocity contribution of this point to the cell velocity
          local_temp_flow = local_temp_flow + norm_vel_temp 

          temp_press_vec = float(pressure_vectors[i_time].GetTuple(node_ind)[0]) # get pressure at this point
          # add pressure contribution of this point to the total cell pressure
          local_temp_press = local_temp_press + temp_press_vec
          
        # Complete the trapezoidal rule integration for this cell by multiplying the sum of velocites
        # by the local area and dividing by the number of vertices
        # Add the contribution of this cell to the curr_flow for the entire outlet face
        curr_flow = curr_flow + local_temp_flow*local_area/3.0

        # Add the contribution of this cell to the curr_press for the entire outlet face
        curr_press = curr_press + local_temp_press*local_area/3.0
      
      # save flow information at the outlet face for the current timestep  
      temp_mean_flow[i_time] = curr_flow

      # save pressure information at the outlet face (with pressure normalized by the total area 
      # of the outlet face) for the current timestep   
      temp_mean_press[i_time] = curr_press/ total_area

      # save pressure drop from inlet to current outlet face at current timestep
      if(file != inlet_file):
        temp_drop_press[i_time] = output_collection_pressure[0][i_time] - (temp_mean_press[i_time]/1333.2239)

    # save pressure waveform amplitude (max press - min press) at current outlet face over current time cycle (all timesteps)
    temp_amp_press = temp_mean_press.max() - temp_mean_press.min()
    
    # save flow information for all timesteps for the current outlet face
    if(file != inlet_file):
      output_collection_flow.append(temp_mean_flow)
    else:
      # invert flow values to make positive
      temp_mean_flow_neg = [-1*f for f in temp_mean_flow]
      output_collection_flow.append(temp_mean_flow_neg)

    # save pressure information for all timesteps for the current outlet face
    if(file != inlet_file):
      # convert from Barye to mmHg
      temp_mean_press_mmHg = [p/1333.2239 for p in temp_mean_press]

      output_collection_pressure.append([temp_mean_press_mmHg, temp_drop_press])
    else:
      # convert from Barye to mmHg
      temp_mean_press_mmHg = [p/1333.2239 for p in temp_mean_press]

      output_collection_pressure.append(temp_mean_press_mmHg)

    # convert from Barye to mmHg
    temp_amp_press_mmHg = temp_amp_press/1333.2239
    
    output_collection_pressure_2.append(temp_amp_press_mmHg)



    #########################################
    #  Calculate TAWSS QoIs on branch wall  #
    #########################################

    # Check for valid branch
    os.chdir('mesh-complete/mesh-surfaces')

    wall_check = glob.glob('wall_' + file + '.vtp')
    
    # return to the home directory
    os.chdir('../..')


    if(len(wall_check) != 0):
      
      # Integrate TAWSS over the surface of the branch to get the mean TAWSS on this branch
      temp_mean_tawss = 0.0

      # Also keep track of max and min TAWSS on the branch
      temp_max_tawss = -sys.maxint
      temp_min_tawss = sys.maxint

      # Also keep track of areas of regions of low TAWSS on the branch
      temp_low_wss_area = np.zeros(len(thresholds))

      # Compute the integral using trapezoidal rule
      total_area = 0.0
      # store WSS information for the entire branch at final time step
      curr_wss = 0.0
          
      # iterate over all mesh cells on branch wall
      for i_cell in xrange(0,wall_model.GetNumberOfCells()):
            
        # extract information about cell vertices
        temp_cell = wall_model.GetCell(i_cell)
        pts_cell = temp_cell.GetPointIds()
        cell_pts = temp_cell.GetPoints()
        p0 = cell_pts.GetPoint(0)
        p1 = cell_pts.GetPoint(1)
        p2 = cell_pts.GetPoint(2)

        # compute area of mesh cell (triangular mesh assumed)
        local_area = vtk.vtkTriangle().TriangleArea(p0,p1,p2)
        total_area = total_area + local_area
            
        local_temp_wss = 0.0
        flag_wss = np.zeros(len(thresholds))

        # add contributions from each vertex of cell
        for ipt in xrange(0, pts_cell.GetNumberOfIds()):
              
          iid = pts_cell.GetId(ipt)     # get node number of this point on wall
          node_ind = wall_nodes[iid]  # get node number on full model

          # get TAWSS magnitude at this point
          temp_wss_vec = float(wss_vectors.GetTuple(node_ind)[0])

          # check for region of low WSS
          flag_wss = flag_wss + [temp_wss_vec > th for th in thresholds]
              
          # add TAWSS contribution of this point to the total cell TAWSS
          local_temp_wss = local_temp_wss + temp_wss_vec

          # update min_wss and max_wss values for this branch
          if temp_wss_vec < temp_min_tawss:
            temp_min_tawss = temp_wss_vec
          if temp_wss_vec > temp_max_tawss:
            temp_max_tawss = temp_wss_vec
                
        # To complete the trapezoidal rule integration, multiply each summed quantity
        # by the area of the cell, then divide by the number of vertices
        # Complete the trapezoidal rule integration for this cell by multiplying the sum of
        # the TAWSS by the local area and dividing by the number of vertices
        # Add the contribution of this cell to the curr_wss for the entire branch wall
        curr_wss = curr_wss + local_temp_wss*local_area/3.0

        # update area of low TAWSS regions for this branch
        for i_th in xrange(0,len(thresholds)):
          if flag_wss[i_th] == 0:
            temp_low_wss_area[i_th] = temp_low_wss_area[i_th] + local_area
        
      # save TAWSS information at the branch wall (with TAWSS normalized by the total area 
      # of the branch wall) for the final timestep
      temp_mean_tawss = float(curr_wss / total_area)

      # precentage of area of low WSS for the final timestep
      temp_low_wss_area = [th / total_area for th in temp_low_wss_area]

      branch_areas.append(total_area)

      # save WSS information for the current branch
      output_collection_tawss.append([temp_mean_tawss, temp_min_tawss, temp_max_tawss, temp_low_wss_area])


    ###################################
    #  Calculate TABP QoIs in branch  #
    ###################################

    # Check if this is a valid branch
    os.chdir('mesh-complete/mesh-surfaces')

    branch_check = glob.glob('wall_' + file + '.vtp')

    # return to the home directory
    os.chdir('../..')

    if(len(branch_check) != 0):

      # Integrate TABP over the volume of the branch to get the mean TABP in this branch
      temp_mean_tabp = 0.0

      # Also keep track of max and min TABP in the branch
      temp_max_tabp = -sys.maxint
      temp_min_tabp = sys.maxint
      
      # Compute the integral using trapezoidal rule
      total_vol = 0.0
      # store TABP information for the entire branch
      curr_tabp = 0.0
          
      # iterate over all mesh cells in branch
      for i_cell in xrange(0,branch_model.GetNumberOfCells()):
            
        # extract information about cell vertices
        temp_cell = branch_model.GetCell(i_cell)
        pts_cell = temp_cell.GetPointIds()
        cell_pts = temp_cell.GetPoints()
        p0 = cell_pts.GetPoint(0)
        p1 = cell_pts.GetPoint(1)
        p2 = cell_pts.GetPoint(2)
        p3 = cell_pts.GetPoint(3)

        # compute volume of mesh cell (tetrahedral mesh assumed)
        local_vol = vtk.vtkTetra().ComputeVolume(p0,p1,p2,p3)
        total_vol = total_vol + local_vol
            
        local_temp_tabp = 0.0
        
        # add contributions from each vertex of cell
        for ipt in xrange(0, pts_cell.GetNumberOfIds()):
              
          iid = pts_cell.GetId(ipt)     # get node number of this point on branch
          node_ind = branch_nodes[iid]  # get node number on full model

          # get TABP at this point
          temp_tabp_vec = float(tabp_vectors.GetTuple(node_ind)[0])
     
          # add TABP contribution of this point to the total cell TABP
          local_temp_tabp = local_temp_tabp + temp_tabp_vec

          # update min_wss and max_wss values for this branch
          if (temp_tabp_vec/1333.2239) < temp_min_tabp:
            temp_min_tabp = (temp_tabp_vec/1333.2239)
          if (temp_tabp_vec/1333.2239) > temp_max_tabp:
            temp_max_tabp = (temp_tabp_vec/1333.2239)
                
        # To complete the trapezoidal rule integration, multiply each summed quantity
        # by the vol of the cell, then divide by the number of vertices
        # Add the contribution of this cell to the curr_tabp for the entire branch
        curr_tabp = curr_tabp + local_temp_tabp*local_vol/4.0
        
      # save TABP information at the branch (with TABP normalized by the total vol of the branch)
      temp_mean_tabp = float(curr_tabp / total_vol)

      branch_vols.append(total_vol)

      # convert from Barye to mmHg
      temp_mean_tabp_mmHg = temp_mean_tabp/1333.2239
      
      # save TABP information for the current branch
      output_collection_tabp.append([temp_mean_tabp_mmHg, temp_min_tabp, temp_max_tabp])
  


  print('Post-processing model QoIs')

  #########################################
  #  Calculate TAWSS QoIs on whole model  #
  #########################################
  model_mean_wss = 0.0
  model_min_wss = sys.maxint
  model_max_wss = -sys.maxint
  model_low_wss_area = np.zeros(len(thresholds))
  model_total_area = 0.0

  for i_file in xrange(0,len(input_filenames)-1):
    # add contribution of this branch to mean TAWSS of model
    model_mean_wss = model_mean_wss + output_collection_tawss[i_file][0]*branch_areas[i_file]
    model_total_area = model_total_area + branch_areas[i_file]

    # update min_wss and max_wss values for this branch at this timestep
    if output_collection_tawss[i_file][1] < model_min_wss:
      model_min_wss = output_collection_tawss[i_file][1]
    if output_collection_tawss[i_file][2] > model_max_wss:
      model_max_wss = output_collection_tawss[i_file][2]

    # add contribution of this branch to area of low TAWSS
    temp_low_wss_area = output_collection_tawss[i_file][3]
    model_low_wss_area = [ml + th*branch_areas[i_file] for ml, th in zip(model_low_wss_area, temp_low_wss_area)]

  # spatially average mean TAWSS in model 
  model_mean_wss = model_mean_wss / model_total_area

  # precentage of area of low TAWSS in model
  model_low_wss_area = [ml / model_total_area for ml in model_low_wss_area]
  
  # save WSS information for the model
  output_collection_tawss.append([model_mean_wss, model_min_wss, model_max_wss, model_low_wss_area])


  ########################################
  #  Calculate TABP QoIs on whole model  #
  ########################################
  model_mean_tabp = 0.0
  model_min_tabp = sys.maxint
  model_max_tabp = -sys.maxint
  model_total_vol = 0.0

  for i_file in xrange(0,len(input_filenames)-1):
    # add contribution of this branch to mean TABP of model
    model_mean_tabp = model_mean_tabp + output_collection_tabp[i_file][0]*branch_vols[i_file]
    model_total_vol = model_total_vol + branch_vols[i_file]

    # update min_wss and max_wss values for this branch at this timestep
    if output_collection_tabp[i_file][1] < model_min_tabp:
      model_min_tabp = output_collection_tabp[i_file][1]
    if output_collection_tabp[i_file][2] > model_max_tabp:
      model_max_tabp = output_collection_tabp[i_file][2]

  # volumetrically average mean TABP in model 
  model_mean_tabp = model_mean_tabp / model_total_vol

  # Convert from Barye to mmHg
  #model_mean_tabp_mmHg = model_mean_tabp/1333.2239

  # save TABP information for the model
  output_collection_tabp.append([model_mean_tabp, model_min_tabp, model_max_tabp]) 


  print('Post-processing strip QoIs')
  #############################################
  #  Calculate TAWSS and TABP QoIs in strips  #
  #############################################

  ##############################################
  #  Find node IDs for strips of the model   #
  ##############################################
  for i_strip in xrange(0, len(strip_clip_planes)):
    # Get geometry (mesh) information for this strip 
    strip_model_vtu = all_results_vtu_reader.GetOutput()
    strip_model_vtp = all_results_vtp_reader.GetOutput()
    
    # Iterate over all planes needed for clipping model
    for i_plane in xrange(0,len(strip_clip_planes[i_strip][0])):  
      curr_point = strip_clip_planes[i_strip][0][i_plane]
      curr_norm = strip_clip_planes[i_strip][1][i_plane]

      clipPlane = vtk.vtkPlane()
      clipPlane.SetOrigin(curr_point)
      clipPlane.SetNormal(curr_norm)

      clipper = vtk.vtkClipDataSet()
      clipper.SetInputData(strip_model_vtu)
      clipper.SetClipFunction(clipPlane)
      clipper.Update()

      clipper2 = vtk.vtkClipPolyData()
      clipper2.SetInputData(strip_model_vtp)
      clipper2.SetClipFunction(clipPlane)
      clipper2.Update()

      triangulator = vtk.vtkDataSetTriangleFilter()
      triangulator.SetInputData(clipper.GetOutput())
      triangulator.Update()

      triangulator2 = vtk.vtkTriangleFilter()
      triangulator2.SetInputData(clipper2.GetOutput())
      triangulator2.Update()

      strip_model_vtu = triangulator.GetOutput()
      strip_model_vtp = triangulator2.GetOutput()


    # Write out the strips to ensure they all look fine
    # strip_writer = vtk.vtkXMLUnstructuredGridWriter()
    # strip_writer.SetInputData(strip_model_vtu)
    # strip_writer.SetFileName('strip_' + str(i_strip) + '.vtu')
    # strip_writer.Write()

    # strip_writer2 = vtk.vtkXMLPolyDataWriter()
    # strip_writer2.SetInputData(strip_model_vtp)
    # strip_writer2.SetFileName('strip_' + str(i_strip) + '.vtp')
    # strip_writer2.Write()

    # Get info for strips
    strip_vtu_numPts = strip_model_vtu.GetNumberOfPoints()           # Determine number of points in the mesh at this outlet
    strip_vtu_IDs = strip_model_vtu.GetPointData().GetArray("GlobalNodeID")  # Extract node IDs to match with full model solution

    strip_vtp_numPts = strip_model_vtp.GetNumberOfPoints()           # Determine number of points in the mesh at this outlet
    strip_vtp_IDs = strip_model_vtp.GetPointData().GetArray("GlobalNodeID")  # Extract node IDs to match with full model solution
    

    # Find the nodes on all_results that correspond to the vtu_strip
    strip_vtu_nodes = []
    for i_node in xrange(0, strip_vtu_numPts):
      strip_vtu_ID = strip_vtu_IDs.GetTuple1(i_node)
      
      # iterate through all nodes in model
      for i_full in xrange(0, all_results_vtu_numPts):
        full_ID = all_results_vtu_IDs.GetTuple1(i_full)

        if(full_ID == strip_vtu_ID):
          strip_vtu_nodes.append(i_full)
          break
    
    # Just to make sure we found all the outlet nodes in all_results
    assert(len(strip_vtu_nodes) == strip_vtu_numPts)


    # Find the nodes on all_results that correspond to the vtp_strip
    strip_vtp_nodes = []
    for i_node in xrange(0, strip_vtp_numPts):
      strip_vtp_ID = strip_vtp_IDs.GetTuple1(i_node)
      
      # iterate through all nodes in model
      for i_full in xrange(0, all_results_vtp_numPts):
        full_ID = all_results_vtp_IDs.GetTuple1(i_full)

        if(full_ID == strip_vtp_ID):
          strip_vtp_nodes.append(i_full)
          break

    # Just to make sure we found all the outlet nodes in all_results
    assert(len(strip_vtp_nodes) == strip_vtp_numPts)


    ########################################
    #  Calculate TAWSS QoIs on strip wall  #
    ########################################

    # Integrate TAWSS over the surface of the branch to get the mean TAWSS on this branch
    temp_mean_tawss = 0.0

    # Also keep track of max and min TAWSS on the branch
    temp_max_tawss = -sys.maxint
    temp_min_tawss = sys.maxint

    # Compute the integral using trapezoidal rule
    total_area = 0.0
    # store WSS information for the entire branch at final time step
    curr_wss = 0.0
        
    # iterate over all mesh cells on branch wall
    for i_cell in xrange(0,strip_model_vtp.GetNumberOfCells()):
          
      # extract information about cell vertices
      temp_cell = strip_model_vtp.GetCell(i_cell)
      pts_cell = temp_cell.GetPointIds()
      cell_pts = temp_cell.GetPoints()
      p0 = cell_pts.GetPoint(0)
      p1 = cell_pts.GetPoint(1)
      p2 = cell_pts.GetPoint(2)

      # compute area of mesh cell (triangular mesh assumed)
      local_area = vtk.vtkTriangle().TriangleArea(p0,p1,p2)
      total_area = total_area + local_area
          
      local_temp_wss = 0.0
      flag_wss = np.zeros(len(thresholds))

      # add contributions from each vertex of cell
      for ipt in xrange(0, pts_cell.GetNumberOfIds()):
            
        iid = pts_cell.GetId(ipt)     # get node number of this point on outlet face
        node_ind = strip_vtp_nodes[iid]  # get node number on full model

        # get TAWSS magnitude at this point
        temp_wss_vec = float(wss_vectors.GetTuple(node_ind)[0])

        # check for region of low WSS
        flag_wss = flag_wss + [temp_wss_vec > th for th in thresholds]
            
        # add TAWSS contribution of this point to the total cell TAWSS
        local_temp_wss = local_temp_wss + temp_wss_vec

        # update min_wss and max_wss values for this branch
        if temp_wss_vec < temp_min_tawss:
          temp_min_tawss = temp_wss_vec
        if temp_wss_vec > temp_max_tawss:
          temp_max_tawss = temp_wss_vec
              
      # To complete the trapezoidal rule integration, multiply each summed quantity
      # by the area of the cell, then divide by the number of vertices
      # Complete the trapezoidal rule integration for this cell by multiplying the sum of
      # the TAWSS by the local area and dividing by the number of vertices
      # Add the contribution of this cell to the curr_wss for the entire branch wall
      curr_wss = curr_wss + local_temp_wss*local_area/3.0
      
    # save TAWSS information at the branch wall (with TAWSS normalized by the total area 
    # of the branch wall) for the final timestep
    temp_mean_tawss = float(curr_wss / total_area)

    # save WSS information for the current branch
    output_collection_tawss.append([temp_mean_tawss, temp_min_tawss, temp_max_tawss])



    ##################################
    #  Calculate TABP QoIs in strip  #
    ##################################

    ## tabp_vectors defined above...

    # Integrate TABP over the volume of the branch to get the mean TABP in this branch
    temp_mean_tabp = 0.0

    # Also keep track of max and min TABP in the branch
    temp_max_tabp = -sys.maxint
    temp_min_tabp = sys.maxint
    
    # Compute the integral using trapezoidal rule
    total_vol = 0.0
    # store TABP information for the entire branch
    curr_tabp = 0.0
        
    # iterate over all mesh cells in branch
    for i_cell in xrange(0,strip_model_vtu.GetNumberOfCells()):
          
      # extract information about cell vertices
      temp_cell = strip_model_vtu.GetCell(i_cell)
      pts_cell = temp_cell.GetPointIds()
      cell_pts = temp_cell.GetPoints()
      p0 = cell_pts.GetPoint(0)
      p1 = cell_pts.GetPoint(1)
      p2 = cell_pts.GetPoint(2)
      p3 = cell_pts.GetPoint(3)

      # compute volume of mesh cell (tetrahedral mesh assumed)
      local_vol = vtk.vtkTetra().ComputeVolume(p0,p1,p2,p3)
      total_vol = total_vol + local_vol
          
      local_temp_tabp = 0.0
      
      # add contributions from each vertex of cell
      for ipt in xrange(0, pts_cell.GetNumberOfIds()):
            
        iid = pts_cell.GetId(ipt)     # get node number of this point on outlet face
        node_ind = strip_vtu_nodes[iid]  # get node number on full model

        # get TABP at this point
        temp_tabp_vec = float(tabp_vectors.GetTuple(node_ind)[0])
   
        # add TABP contribution of this point to the total cell TABP
        local_temp_tabp = local_temp_tabp + temp_tabp_vec

        # update min_wss and max_wss values for this branch
        if (temp_tabp_vec/1333.2239) < temp_min_tabp:
          temp_min_tabp = (temp_tabp_vec/1333.2239)
        if (temp_tabp_vec/1333.2239) > temp_max_tabp:
          temp_max_tabp = (temp_tabp_vec/1333.2239)
              
      # To complete the trapezoidal rule integration, multiply each summed quantity
      # by the vol of the cell, then divide by the number of vertices
      # Add the contribution of this cell to the curr_tabp for the entire branch
      curr_tabp = curr_tabp + local_temp_tabp*local_vol/4.0
      
    # save TABP information at the branch (with TABP normalized by the total vol of the branch)
    temp_mean_tabp = float(curr_tabp / total_vol)


    # Convert from Barye to mmHg
    temp_mean_tabp_mmHg = temp_mean_tabp/1333.2239

    # save TABP information for the current branch
    output_collection_tabp.append([temp_mean_tabp_mmHg, temp_min_tabp, temp_max_tabp])



  print('Writing results to file')
  ####################################
  #  Write out all results to files  #
  #################################### 

  # Now that we have looped over all our files of interest and integrated
  # the variables, it is time to save them to the output file. 

  # Flow QoIs
  outfile = open(output_filename_flow, 'w')
  
  # First print a header that tells what each integrated quantity of interest is
  out_string = 'Time_Step'
  for i_file in xrange(0, len(input_filenames)):
    out_string = out_string + ' ' + input_filenames[i_file] + '_cap_flow'
  out_string = out_string + '\n'
  outfile.write(out_string)

  # Now print the data for each quantity of interest at each time step
  for i_time in xrange(0, timestep_count):
    # Print time step
    out_string = str(i_time)
    
    # Print each quantity of interest at that timestep
    for i_file in xrange(0,len(input_filenames)):  
      out_string = out_string + ' ' + str(output_collection_flow[i_file][i_time])
        
    out_string = out_string + '\n'
    outfile.write(out_string)

  outfile.close()


  # Pressure QoIs (individual timesteps) 
  outfile = open(output_filename_pressure, 'w')
  
  # First print a header that tells what each integrated quantity of interest is
  out_string = 'Time_Step'
  for i_file in xrange(0, len(input_filenames)):
    out_string = out_string + ' ' + input_filenames[i_file] + '_cap_pressure'
    if(input_filenames[i_file] != inlet_file):
      out_string = out_string + ' ' + input_filenames[i_file] + '_pressure_drop'
  out_string = out_string + '\n'
  outfile.write(out_string)

  # Now print the data for each quantity of interest at each time step
  for i_time in xrange(0, timestep_count):
      # Print time step
      out_string = str(i_time)
      
      # Print each quantity of interest at that timestep
      for i_file in xrange(0,len(input_filenames)):
        if(input_filenames[i_file] != inlet_file):  
          out_string = out_string + ' ' + str(output_collection_pressure[i_file][0][i_time]) \
                       + ' ' + str(output_collection_pressure[i_file][1][i_time])
        else:
          out_string = out_string + ' ' + str(output_collection_pressure[i_file][i_time])

      out_string = out_string + '\n'
      outfile.write(out_string)

  outfile.close() 


  # Pressure QoIs (all timesteps)
  outfile = open(output_filename_pressure_2, 'w')
  
  # First print a header that tells what each integrated quantity of interest is
  out_string = ''
  for i_file in xrange(0, len(input_filenames)):
    out_string = out_string + input_filenames[i_file] + '_cap_pressure_amplitude '
  out_string = out_string + '\n'
  outfile.write(out_string)

  # Now print the data for each quantity of interest
  out_string = ''
  for i_file in xrange(0,len(input_filenames)):
    out_string = out_string + str(output_collection_pressure_2[i_file]) + ' '

  out_string = out_string + '\n'
  outfile.write(out_string)

  outfile.close() 


  # TABP QoIs (all timesteps)
  outfile = open(output_filename_tabp, 'w')

  # First print a header that tells what each integrated quantity of interest is
  out_string = ''
  # branches
  for i_file in xrange(1, len(input_filenames)):
    out_string = out_string + input_filenames[i_file] + '_meanTABP ' + \
                 input_filenames[i_file] + '_minTABP ' + \
                 input_filenames[i_file] + '_maxTABP '
  # model
  out_string = out_string + 'model_meanTABP model_minTABP model_maxTABP '
  # strips
  for i_strip in xrange(0, len(strip_clip_planes)):
    out_string = out_string + 'strip_' + str(i_strip) + '_meanTABP ' + \
                 'strip_' + str(i_strip) + '_minTABP ' + \
                 'strip_' + str(i_strip) + '_maxTABP '
  out_string = out_string + '\n'
  outfile.write(out_string)

  # Now print the data for each quantity of interest
  out_string = ''
  qoi_count = 0

  # Print each quantity of interest (branches)
  for i_file in xrange(1,len(input_filenames)):
    out_string = out_string + str(output_collection_tabp[i_file][0]) \
                 + ' ' + str(output_collection_tabp[i_file][1]) \
                 + ' ' + str(output_collection_tabp[i_file][2]) + ' ' 
    qoi_count += 1
  
  # Print each quantity of interest (model)
  out_string = out_string + str(output_collection_tabp[qoi_count][0]) \
                 + ' ' + str(output_collection_tabp[qoi_count][1]) \
                 + ' ' + str(output_collection_tabp[qoi_count][2]) + ' ' 
  qoi_count += 1

  # Print each quantity of interest (strips)
  for i_ in xrange(0,len(strip_clip_planes)):
    out_string = out_string + str(output_collection_tabp[qoi_count][0]) \
                 + ' ' + str(output_collection_tabp[qoi_count][1]) \
                 + ' ' + str(output_collection_tabp[qoi_count][2]) + ' ' 
    qoi_count += 1
  out_string = out_string + '\n'
  outfile.write(out_string)

  outfile.close()


 # TAWSS QoIs (all timesteps)
  outfile = open(output_filename_tawss, 'w')

  # First print a header that tells what each integrated quantity of interest is
  out_string = ''
  # branches
  for i_file in xrange(1, len(input_filenames)):
    out_string = out_string + input_filenames[i_file] + '_meanTAWSS ' + \
                 input_filenames[i_file] + '_minTAWSS ' + \
                 input_filenames[i_file] + '_maxTAWSS '
  # model
  out_string = out_string + 'model_meanTAWSS model_minTAWSS model_maxTAWSS '
  for i_th in xrange(0,len(thresholds)):
    out_string = out_string + 'model_lowWSS_th=' + str(thresholds[i_th]) + ' '
  # strips
  for i_strip in xrange(0, len(strip_clip_planes)):
    out_string = out_string + 'strip_' + str(i_strip) + '_meanTAWSS ' + \
                 'strip_' + str(i_strip) + '_minTAWSS ' + \
                 'strip_' + str(i_strip) + '_maxTAWSS '
  out_string = out_string + '\n'
  outfile.write(out_string)

  # Now print the data for each quantity of interest
  out_string = ''
  qoi_count = 0

  # Print each quantity of interest (branches)
  for i_file in xrange(1,len(input_filenames)):
    out_string = out_string + str(output_collection_tawss[i_file][0]) \
                 + ' ' + str(output_collection_tawss[i_file][1]) \
                 + ' ' + str(output_collection_tawss[i_file][2]) + ' ' 
    qoi_count += 1
  
  # Print each quantity of interest (model)
  out_string = out_string + str(output_collection_tawss[qoi_count][0]) \
                 + ' ' + str(output_collection_tawss[qoi_count][1]) \
                 + ' ' + str(output_collection_tawss[qoi_count][2]) + ' ' 
  for i_th in xrange(0,len(thresholds)):
    out_string = out_string + str(output_collection_tawss[qoi_count][3][i_th]) + ' '
  qoi_count += 1

  # Print each quantity of interest (strips)
  for i_ in xrange(0,len(strip_clip_planes)):
    out_string = out_string + str(output_collection_tawss[qoi_count][0]) \
                 + ' ' + str(output_collection_tawss[qoi_count][1]) \
                 + ' ' + str(output_collection_tawss[qoi_count][2]) + ' ' 
    qoi_count += 1
  out_string = out_string + '\n'
  outfile.write(out_string)

  outfile.close()


