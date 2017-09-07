# ****************************************************************************
#
# Formulation 1A for time domain noise prediction through FW-H method.
# Autor: Jaime Vaquero
# May 2017
#
# ****************************************************************************
import numpy as np
#import numba
#from math import *
#from cmath import *


#@numba.jit(numba.types.Tuple((int,int,int,np.zeros(1),np.zeros([1,1]),np.zeros([1,1]),\
#np.zeros(1),np.zeros([1,1,1]),np.zeros([1,1]),np.zeros([1,1,1]),np.zeros([1,1,1]),\
#np.zeros([1,1,1]),np.zeros([1,1,1]),np.zeros([1,1,1]),np.zeros([1,1,1]),float,float,float,float))\
#(numba.typeof(np.zeros([1,1,1])),numba.typeof(np.zeros([1,1,1])),numba.typeof(np.zeros([1,1,1]))),\
#nopython = True)

def Formulation1A(LengthTime,LengthSpaceX,LengthSpaceExtZ,tsrc,SpaceX,SpaceY,SpaceZ,normal,elesurf,FunctionV1,FunctionV2, \
                  FunctionV3,FunctionP,FunctionRho,Xobs,dt,c0,p0,rho0):

    # INITIALISATION

    # Pressure
    p    = np.zeros(LengthTime)

    # Density
    rho    = np.zeros(LengthTime)

    # Source coordinates
    y    = np.zeros([3,LengthTime])
    
    x    = np.zeros([3,LengthTime])

    # Flow velocities at surface f=0
    unm4  = np.zeros([3,LengthTime-4])
    unm3  = np.zeros([3,LengthTime-4])
    unm2  = np.zeros([3,LengthTime-4])
    unm1  = np.zeros([3,LengthTime-4])
    u_int = np.zeros([3,LengthTime-4])
    u     = np.zeros([3,LengthTime])

    # Surface f=0 velocities
    v    = np.zeros([3,LengthTime])

    # Normal flow velocities at surface f=0
    un     = np.zeros(LengthTime)

    # Normal surface f=0 velocities at surface f=0
    vn     = np.zeros(LengthTime)

    # Normal vector at surface f=0
    n_nm4  = np.zeros([3,LengthTime-4]) # The direction is constant along the surface. If it is not, source position dependance must be added.
    n_nm3  = np.zeros([3,LengthTime-4])
    n_nm2  = np.zeros([3,LengthTime-4])
    n_nm1  = np.zeros([3,LengthTime-4])
    n_int  = np.zeros([3,LengthTime-4])
    n      = np.zeros([3,LengthTime])
    # Flow rate
    Unm4  = np.zeros([3,LengthTime-4]) # The direction is constant along the surface. If it is not, source position dependance must be added.
    Unm3  = np.zeros([3,LengthTime-4])
    Unm2  = np.zeros([3,LengthTime-4])
    Unm1  = np.zeros([3,LengthTime-4])
    U_int = np.zeros([3,LengthTime-4])
    U     = np.zeros([3,LengthTime])

    # Force
    Lnm4  = np.zeros([3,LengthTime-4]) # The direction is constant along the surface. If it is not, source position dependance must be added.
    Lnm3  = np.zeros([3,LengthTime-4])
    Lnm2  = np.zeros([3,LengthTime-4])
    Lnm1  = np.zeros([3,LengthTime-4])
    L_int = np.zeros([3,LengthTime-4])
    L     = np.zeros([3,LengthTime])

    # Output initialisation
    tobs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])
    robs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])
    pobs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime]) # Attention: la pression acoustique devrait dependre de la position obs! Oui,
                                              # En fait, on resout ce code pour un seul observateur, donc on va obtenir
                                              # p = f(sourcepos, t) qu'on devra correctement additioner plus tard
    ri      = np.zeros([3,LengthTime])
    r       = np.zeros(LengthTime)
    Mi      = np.zeros([3,LengthTime])
    M       = np.zeros(LengthTime)  
    Mr      = np.zeros(LengthTime)  
    Lr      = np.zeros(LengthTime)
    Lm      = np.zeros(LengthTime)  
    
    dUdt    = np.zeros([3,LengthTime])
    dLdt    = np.zeros([3,LengthTime])
    dndt    = np.zeros([3,LengthTime])
    dudt    = np.zeros([3,LengthTime])
    
    U_m1 = np.zeros(3)
    L_m1 = np.zeros(3)
    n_m1 = np.zeros(3)
    u_m1 = np.zeros(3)
    
    U_m2 = np.zeros(3)
    L_m2 = np.zeros(3)
    n_m2 = np.zeros(3)
    u_m2 = np.zeros(3)

    U_m3 = np.zeros(3)
    L_m3 = np.zeros(3)
    n_m3 = np.zeros(3)
    u_m3 = np.zeros(3)
    
    U_m4 = np.zeros(3)
    L_m4 = np.zeros(3)
    n_m4 = np.zeros(3)
    u_m4 = np.zeros(3)
    
    dUdt_m1 = np.zeros(3)
    dLdt_m1 = np.zeros(3)
    dndt_m1 = np.zeros(3)
    dudt_m1 = np.zeros(3)
    
    dUdt_m2 = np.zeros(3)
    dLdt_m2 = np.zeros(3)
    dndt_m2 = np.zeros(3)
    dudt_m2 = np.zeros(3)
    
    dUdt_m3 = np.zeros(3)
    dLdt_m3 = np.zeros(3)
    dndt_m3 = np.zeros(3)
    dudt_m3 = np.zeros(3)
    
    dUdt_m4 = np.zeros(3)
    dLdt_m4 = np.zeros(3)
    dndt_m4 = np.zeros(3)
    dudt_m4 = np.zeros(3)
    
    dMidt   = np.zeros([3,LengthTime])
    
    pQ1     = np.zeros(LengthTime)
    pQ2     = np.zeros(LengthTime)
    pL1     = np.zeros(LengthTime)
    pL2     = np.zeros(LengthTime)
    pL3     = np.zeros(LengthTime)
    
    for l in range(0,LengthSpaceExtZ):
        print('  Formulation 1A step: ' + str(l+1) + ' of ' + str(LengthSpaceExtZ))
        for j in range(0,LengthSpaceX):

            n[0,:] = normal[0,j,l]
            n[1,:] = normal[1,j,l]
            n[2,:] = normal[2,j,l]
            y[0,:] = SpaceX[j,l]
            y[1,:] = SpaceY[j,l]
            y[2,:] = SpaceZ[l]
            x[0,:] = Xobs[0]
            x[1,:] = Xobs[1]
            x[2,:] = Xobs[2]            

            u[0,:] = FunctionV1[j,l,:]
            u[1,:] = FunctionV2[j,l,:]
            u[2,:] = FunctionV3[j,l,:]
            v[0,:] = 0.
            v[1,:] = 0.
            v[2,:] = 0.
            un   = np.einsum('ij,ij->j',n,u)
            vn   = np.einsum('ij,ij->j',n,v)         

            p    = FunctionP[j,l,:]
            rho  = FunctionRho[j,l,:]

            
            ri    = x - y
            r[:]  = np.sqrt(ri[0,:]**2. + ri[1,:]**2. + ri[2,:]**2.)
            Mi    = u / c0
            M[:]  = np.sqrt(Mi[0,:]**2. + Mi[1,:]**2. + Mi[2,:]**2.)
            Mr    = np.einsum('ij,ij->j',Mi,ri/r)            
            U     = (1. - rho/rho0) * v + rho/rho0 * u
            L     = (p - p0) * n + rho * u * (un - vn) # Viscosity tensor is not considered
            Un    = np.einsum('ij,ij->j',U,n)
            Lr    = np.einsum('ij,ij->j',L,ri/r)
            Lm    = np.einsum('ij,ij->j',L,Mi)
            
            U_int = U[:,4:]
            L_int = L[:,4:]
            n_int = n[:,4:]
            u_int = u[:,4:]
            
            Unm1   = U[:,3:-1]
            Lnm1   = L[:,3:-1]
            n_nm1  = n[:,3:-1]
            unm1   = u[:,3:-1]      
            
            Unm2   = U[:,2:-2]
            Lnm2   = L[:,2:-2]
            n_nm2  = n[:,2:-2]
            unm2   = u[:,2:-2]
            
            Unm3   = U[:,1:-3]
            Lnm3   = L[:,1:-3]
            n_nm3  = n[:,1:-3]
            unm3   = u[:,1:-3]

            Unm4   = U[:,0:-4]
            Lnm4   = L[:,0:-4]
            n_nm4  = n[:,0:-4]
            unm4   = u[:,0:-4]            

            U_m1 = U[:,3]
            L_m1 = L[:,3]
            n_m1 = n[:,3]
            u_m1 = u[:,3]

            U_m2 = U[:,2]
            L_m2 = L[:,2]
            n_m2 = n[:,2]
            u_m2 = u[:,2]

            U_m3 = U[:,1]
            L_m3 = L[:,1]
            n_m3 = n[:,1]
            u_m3 = u[:,1]
            
            U_m4 = U[:,0]
            L_m4 = L[:,0]
            n_m4 = n[:,0]
            u_m4 = u[:,0]  
            
            
            # Time derivatives
            dUdt_m3 = (U_m3 - U_m4)/dt
            dLdt_m3 = (L_m3 - L_m4)/dt
            dndt_m3 = (n_m3 - n_m4)/dt
            dudt_m3 = (u_m3 - u_m4)/dt
            
            dUdt_m4 = dUdt_m3
            dLdt_m4 = dLdt_m3
            dndt_m4 = dndt_m3
            dudt_m4 = dudt_m3
            
            dUdt_m2 = (1.5*U_m2 -2.*U_m3 +0.5*U_m4)/dt
            dLdt_m2 = (1.5*L_m2 -2.*L_m3 +0.5*L_m4)/dt
            dndt_m2 = (1.5*n_m2 -2.*n_m3 +0.5*n_m4)/dt
            dudt_m2 = (1.5*u_m2 -2.*u_m3 +0.5*u_m4)/dt

            dUdt_m1 = (11./6.*U_m1 -3.*U_m2 + 1.5*U_m3  -1./3.*U_m4)/dt
            dLdt_m1 = (11./6.*L_m1 -3.*L_m2 + 1.5*L_m3  -1./3.*L_m4)/dt
            dndt_m1 = (11./6.*n_m1 -3.*n_m2 + 1.5*n_m3  -1./3.*n_m4)/dt
            dudt_m1 = (11./6.*u_m1 -3.*u_m2 + 1.5*u_m3  -1./3.*u_m4)/dt

            dUdt_int = (25./12.*U_int - 4.*Unm1  + 3.*Unm2  - 4./3.*Unm3  + 0.25*Unm4 )/dt
            dLdt_int = (25./12.*L_int - 4.*Lnm1  + 3.*Lnm2  - 4./3.*Lnm3  + 0.25*Lnm4 )/dt
            dndt_int = (25./12.*n_int - 4.*n_nm1 + 3.*n_nm2 - 4./3.*n_nm3 + 0.25*n_nm4)/dt
            dudt_int = (25./12.*u_int - 4.*unm1  + 3.*unm2  - 4./3.*unm3  + 0.25*unm4 )/dt
            
            for i in range(0,3):
                dUdt[i,:] = np.concatenate([[dUdt_m4[i],dUdt_m3[i],dUdt_m2[i],dUdt_m1[i]],dUdt_int[i,:]])
                dLdt[i,:] = np.concatenate([[dLdt_m4[i],dLdt_m3[i],dLdt_m2[i],dLdt_m1[i]],dLdt_int[i,:]])
                dndt[i,:] = np.concatenate([[dndt_m4[i],dndt_m3[i],dndt_m2[i],dndt_m1[i]],dndt_int[i,:]])
                dudt[i,:] = np.concatenate([[dudt_m4[i],dudt_m3[i],dudt_m2[i],dudt_m1[i]],dudt_int[i,:]])
            

            dMidt = dudt / c0
            dMdtr = np.einsum('ij,ij->j',dMidt,ri/r)
            Undt  = np.einsum('ij,ij->j',U,dndt)
            dUdtn = np.einsum('ij,ij->j',dUdt,n)
            dLdtr = np.einsum('ij,ij->j',dLdt,ri/r) 

            # Computation of pQ1, pQ2, pL1, pL2, and pL3
            pQ1[:] = 1./(4. * np.pi)      * (rho0 * (dUdtn + Undt)) /\
                                                (r * (1. - Mr)**2.)             * elesurf[j,l]
            pQ2[:] = 1./(4. * np.pi)      * (rho0 * Un * (r * dMdtr + c0*(Mr - M**2.))) /\
                                                (r**2. * (1. - Mr)**3.)         * elesurf[j,l]
            pL1[:] = 1./(4. * np.pi * c0) * (dLdtr) / (r * (1. - Mr)**2.)          * elesurf[j,l]
            pL2[:] = 1./(4. * np.pi)      * (Lr - Lm) / (r**2. * (1. - Mr)**2.)    * elesurf[j,l]
            pL3[:] = 1./(4. * np.pi * c0) * (Lr * (r * dMdtr + c0 * (Mr - M**2.))) / \
                                                (r**2. * (1. - Mr)**3.)         * elesurf[j,l]
                                            
            pobs[j,l,:] = pQ1 + pQ2 + pL1 + pL2 + pL3
            tobs[j,l,:] = tsrc[:] + r/c0
            robs[j,l,:] = r

    return pobs,tobs,robs
