# ****************************************************************************
#
# Formulation 1A for time domain noise prediction through FW-H method.
# Autor: Jaime Vaquero
# July 2017
# Margnat 2010 algorithm
#
# ****************************************************************************
#import numpy as np
import numba

@numba.jit(nopython=True)
def Margnat2010(LengthTime,LengthTimeAdv,LengthSpaceX,LengthSpaceZ,ObserverPosition,\
                ObserverPressure,c0,dt,ObserverPressure_discret,ObserverTime):
     
    for k in range(0,LengthTime-1):
        for j in range(0,LengthSpaceX):
            for l in range(0,LengthSpaceZ):
                kadvM = k + 1 + int(ObserverPosition[j,l,k+1] / (c0*dt))
                wM = 1 - (ObserverPosition[j,l,k+1] / (c0*dt) - int(ObserverPosition[j,l,k+1] / (c0*dt)))
                pinter = (ObserverPressure[j,l,k+1] - ObserverPressure[j,l,k]) * wM + ObserverPressure[j,l,k]
                ObserverPressure_discret[j,l,kadvM] = pinter
                ObserverTime[kadvM] = kadvM * dt    
                
    return ObserverPressure_discret,ObserverTime