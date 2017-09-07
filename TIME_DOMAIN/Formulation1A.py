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
                  FunctionV3,FunctionP,FunctionRho,x,dt,c0,p0,rho0):

    # INITIALISATION

    # Pressure
    pnm4 = 0.
    pnm3 = np.copy(pnm4)
    pnm2 = np.copy(pnm3)
    pnm1 = np.copy(pnm2)
    p    = np.copy(pnm1)

    # Density
    rhonm4 = 0.
    rhonm3 = np.copy(rhonm4)
    rhonm2 = np.copy(rhonm3)
    rhonm1 = np.copy(rhonm2)
    rho    = np.copy(rhonm1)

    # Source coordinates
    ynm4 = np.zeros(3)
    ynm3 = np.copy(ynm4)
    ynm2 = np.copy(ynm3)
    ynm1 = np.copy(ynm2)
    y    = np.copy(ynm1)

    # Flow velocities at surface f=0
    unm4 = np.zeros(3)
    unm3 = np.copy(unm4)
    unm2 = np.copy(unm3)
    unm1 = np.copy(unm2)
    u    = np.copy(unm1)

    # Surface f=0 velocities
    vnm4 = np.zeros(3)
    vnm3 = np.copy(vnm4)
    vnm2 = np.copy(vnm3)
    vnm1 = np.copy(vnm2)
    v    = np.copy(vnm1)

    # Normal flow velocities at surface f=0
    un_nm4 = 0.
    un_nm3 = np.copy(un_nm4)
    un_nm2 = np.copy(un_nm3)
    un_nm1 = np.copy(un_nm2)
    un     = np.copy(un_nm1)

    # Normal surface f=0 velocities at surface f=0
    vn_nm4 = 0.
    vn_nm3 = np.copy(vn_nm4)
    vn_nm2 = np.copy(vn_nm3)
    vn_nm1 = np.copy(vn_nm2)
    vn     = np.copy(vn_nm1)

    # Normal vector at surface f=0
    n_nm4 = np.zeros(3) # The direction is constant along the surface. If it is not, source position dependance must be added.
    n_nm3 = np.copy(n_nm4)
    n_nm2 = np.copy(n_nm3)
    n_nm1 = np.copy(n_nm2)
    n     = np.copy(n_nm1)

    # Flow rate
    Unm4 = np.zeros(3) # The direction is constant along the surface. If it is not, source position dependance must be added.
    Unm3 = np.copy(Unm4)
    Unm2 = np.copy(Unm3)
    Unm1 = np.copy(Unm2)
    U    = np.copy(Unm1)

    # Force
    Lnm4 = np.zeros(3) # The direction is constant along the surface. If it is not, source position dependance must be added.
    Lnm3 = np.copy(Lnm4)
    Lnm2 = np.copy(Lnm3)
    Lnm1 = np.copy(Lnm2)
    L    = np.copy(Lnm1)

    # Output initialisation
    tobs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])
    robs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])
    pobs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime]) # Attention: la pression acoustique devrait dependre de la position obs! Oui,
                                              # En fait, on resout ce code pour un seul observateur, donc on va obtenir
                                              # p = f(sourcepos, t) qu'on devra correctement additioner plus tard
    ri    = np.zeros(3)
    Mi    = np.zeros(3)
    dUdt  = np.zeros(3)
    dLdt  = np.zeros(3)
    dndt  = np.zeros(3)
    dudt  = np.zeros(3)
    dMidt = np.zeros(3)
    
    for l in range(0,LengthSpaceExtZ):
        print('  Formulation 1A step: ' + str(l+1) + ' of ' + str(LengthSpaceExtZ))
        for j in range(0,LengthSpaceX):
            for k in range(0,LengthTime):
                # Saving values for time derivatives
                pnm4 = np.copy(pnm3)
                pnm3 = np.copy(pnm2)
                pnm2 = np.copy(pnm1)
                pnm1 = np.copy(p)
        
                rhonm4 = np.copy(rhonm3)
                rhonm3 = np.copy(rhonm2)
                rhonm2 = np.copy(rhonm1)
                rhonm1 = np.copy(rho)
        
                ynm4 = np.copy(ynm3)
                ynm3 = np.copy(ynm2)
                ynm2 = np.copy(ynm1)
                ynm1 = np.copy(y)
        
                unm4 = np.copy(unm3)
                unm3 = np.copy(unm2)
                unm2 = np.copy(unm1)
                unm1 = np.copy(u)
        
                vnm4 = np.copy(vnm3)
                vnm3 = np.copy(vnm2)
                vnm2 = np.copy(vnm1)
                vnm1 = np.copy(v)
        
                un_nm4 = np.copy(un_nm3)
                un_nm3 = np.copy(un_nm2)
                un_nm2 = np.copy(un_nm1)
                un_nm1 = np.copy(un)
        
                vn_nm4 = np.copy(vn_nm3)
                vn_nm3 = np.copy(vn_nm2)
                vn_nm2 = np.copy(vn_nm1)
                vn_nm1 = np.copy(vn)
        
                n_nm4 = np.copy(n_nm3)
                n_nm3 = np.copy(n_nm2)
                n_nm2 = np.copy(n_nm1)
                n_nm1 = np.copy(n)
        
                Unm4 = np.copy(Unm3)
                Unm3 = np.copy(Unm2)
                Unm2 = np.copy(Unm1)
                Unm1 = np.copy(U)
        
                Lnm4 = np.copy(Lnm3)
                Lnm3 = np.copy(Lnm2)
                Lnm2 = np.copy(Lnm1)
                Lnm1 = np.copy(L)


                n[:] = normal[:,j,l]
                y[0] = SpaceX[j,l]
                y[1] = SpaceY[j,l]
                y[2] = SpaceZ[l]

                u[0] = FunctionV1[j,l,k]
                u[1] = FunctionV2[j,l,k]
                u[2] = FunctionV3[j,l,k]
                v[0] = 0.
                v[1] = 0.
                v[2] = 0.
                un   = np.dot(n,u[:])
                vn   = np.dot(n,v[:])

                p    = FunctionP[j,l,k]
                rho  = FunctionRho[j,l,k]


                ri[:] = x[:] - y[:]
                r     = np.sqrt(ri[0]**2. + ri[1]**2. + ri[2]**2.)
                Mi[:] = u[:] / c0
                M     = np.sqrt(Mi[0]**2. + Mi[1]**2. + Mi[2]**2.)
                Mr    = np.dot(Mi,ri/r)           
                U[:]  = (1. - rho/rho0) * v[:] + rho/rho0 * u[:]
                L[:]  = (p - p0) * n[:] + rho * u[:] * (un - vn) # Viscosity tensor is not considered
                Un    = np.dot(U,n)
                Lr    = np.dot(L,ri/r)
                Lm    = np.dot(L,Mi)
                

                # Time derivatives
                if k == 1:
                    dUdt[:] = (U[:] - Unm1[:])/dt
                    dLdt[:] = (L[:] - Lnm1[:])/dt
                    dndt[:] = (n[:] - n_nm1[:])/dt
                    dudt[:] = (u[:] - unm1[:])/dt
                elif k == 2:
                    dUdt[:] = (1.5*U[:] -2.*Unm1[:]  +0.5*Unm2[:] )/dt
                    dLdt[:] = (1.5*L[:] -2.*Lnm1[:]  +0.5*Lnm2[:] )/dt
                    dndt[:] = (1.5*n[:] -2.*n_nm1[:] +0.5*n_nm2[:])/dt
                    dudt[:] = (1.5*u[:] -2.*unm1[:]  +0.5*unm2[:] )/dt
                elif k == 3:
                    dUdt[:] = (11./6.*U[:] -3.*Unm1[:]  + 1.5*Unm2[:]  -1./3.*Unm3[:] )/dt
                    dLdt[:] = (11./6.*L[:] -3.*Lnm1[:]  + 1.5*Lnm2[:]  -1./3.*Lnm3[:] )/dt
                    dndt[:] = (11./6.*n[:] -3.*n_nm1[:] + 1.5*n_nm2[:] -1./3.*n_nm3[:])/dt
                    dudt[:] = (11./6.*u[:] -3.*unm1[:]  + 1.5*unm2[:]  -1./3.*unm3[:] )/dt
                elif k >= 4:
                    dUdt[:] = (25./12.*U[:] - 4.*Unm1[:]  + 3.*Unm2[:]  - 4./3.*Unm3[:]  + 0.25*Unm4[:] )/dt
                    dLdt[:] = (25./12.*L[:] - 4.*Lnm1[:]  + 3.*Lnm2[:]  - 4./3.*Lnm3[:]  + 0.25*Lnm4[:] )/dt
                    dndt[:] = (25./12.*n[:] - 4.*n_nm1[:] + 3.*n_nm2[:] - 4./3.*n_nm3[:] + 0.25*n_nm4[:])/dt
                    dudt[:] = (25./12.*u[:] - 4.*unm1[:]  + 3.*unm2[:]  - 4./3.*unm3[:]  + 0.25*unm4[:] )/dt

                dMidt[:] = dudt[:] / c0
                dMdtr    = np.dot(dMidt,ri/r)
                Undt     = np.dot(U,dndt)
                dUdtn    = np.dot(dUdt,n)
                dLdtr    = np.dot(dLdt,ri/r)  
                

                # Computation of pQ1, pQ2, pL1, pL2, and pL3
                pQ1 = 1./(4. * np.pi)      * (rho0 * (dUdtn + Undt)) /\
                                                    (r * (1. - Mr)**2.)             * elesurf[j,l]
                pQ2 = 1./(4. * np.pi)      * (rho0 * Un * (r * dMdtr + c0*(Mr - M**2.))) /\
                                                    (r**2. * (1. - Mr)**3.)         * elesurf[j,l]
                pL1 = 1./(4. * np.pi * c0) * (dLdtr) / (r * (1. - Mr)**2.)          * elesurf[j,l]
                pL2 = 1./(4. * np.pi)      * (Lr - Lm) / (r**2. * (1. - Mr)**2.)    * elesurf[j,l]
                pL3 = 1./(4. * np.pi * c0) * (Lr * (r * dMdtr + c0 * (Mr - M**2.))) / \
                                                    (r**2. * (1. - Mr)**3.)         * elesurf[j,l]

                pobs[j,l,k] = pQ1 + pQ2 + pL1 + pL2 + pL3
                tobs[j,l,k] = tsrc[k] + r/c0
                robs[j,l,k] = r

    return pobs,tobs,robs
