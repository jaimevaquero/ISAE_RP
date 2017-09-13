# ****************************************************************************
#
# Function for surface reading in the advanced time approach of the 
# FW-H method.
# Author: Jaime Vaquero
# September 2017
#
# ****************************************************************************
import numpy as np

def ReadSurface3D(x,y,z,variables,Cartesian=True):
    if (Cartesian == True):
        print(' ReadSurface3D: using cartesian coordinates')
    else: 
        print(' ReadSurface3D: using cylindrical coordinates')
   

def SortData(variables,sortingVariable):
    sortindex  = np.argsort(sortingVariable)
    outputList = [] 
    for x in variables:
       x2 = np.zeros(length(x))
       i = -1
       for m in sortindex:
           i     = i + 1
           x2[i] = x[m] 
       outputList.append(x) 
    return outputList
