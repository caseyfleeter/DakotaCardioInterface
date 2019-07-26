#!/bin/bash

# chmod up folders
chmod 700 Run_0D_Model/*.sh

cp /path/to/OneDSolver Run_1D_Model/1D_Files/
chmod 700 Run_1D_Model/1D_Files/OneDSolver
chmod 700 Run_1D_Model/*.sh

cp /path/to/svpre-gcc-gfortran.exe Run_3D_Model/3D_Files/svpre
cp /path/to/svsolver-gcc-gfortran-mpich.exe Run_3D_Model/3D_Files/svsolver
cp /path/to/svpost-gcc-gfortran.exe Run_3D_Model/3D_Files/svpost
chmod 700 Run_3D_Model/3D_Files/sv*
chmod 700 Run_3D_Model/*.sh

chmod 700 runDakota.sh
cp /path/to/dakota .
chmod 700 dakota
