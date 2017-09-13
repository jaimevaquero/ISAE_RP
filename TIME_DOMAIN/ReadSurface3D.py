# ****************************************************************************
#
# Function for surface reading in the advanced time approach of the 
# FW-H method.
# Autor: Jaime Vaquero
# September 2017
#
# ****************************************************************************
import numpy as np

def ReadSurface3D(x,y,z,variables,Cartesian=True):
    if (Cartesian == True):
        print(' ReadSurface3D: using cartesian coordinates')
    else: 
        print(' ReadSurface3D: using cylindrical coordinates')
              
