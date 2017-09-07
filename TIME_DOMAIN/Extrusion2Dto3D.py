# ****************************************************************************
#
# Function Extrusion: gives back data on a surface defined by a line extruded
# in Z direction, converting thus 2D data to 3D data.
#
# Autor: Jaime Vaquero
# June 2017
#
# ****************************************************************************
import numpy as np

def Extrusion2Dto3D(nxs,nt,FunctionV1,FunctionV2,FunctionV3,FunctionP,FunctionRho,xs,ys,zs,\
                    normal,elesurfXY,nzExtr,Lhalfz):
    nzs        = nzExtr
    xsExt      = np.zeros([nxs,nzExtr])
    ysExt      = np.zeros([nxs,nzExtr])
    zsExt      = np.zeros(nzExtr)
    Aext_u     = np.zeros([nxs,nzExtr,nt])
    Aext_v     = np.zeros([nxs,nzExtr,nt])
    Aext_w     = np.zeros([nxs,nzExtr,nt])
    Aext_p     = np.zeros([nxs,nzExtr,nt])
    Aext_rho   = np.zeros([nxs,nzExtr,nt])
    elesurfExt = np.zeros([nxs,nzExtr])
    normalExt  = np.zeros([3,nxs,nzExtr])

    for i in range(0,nxs):
        xsExt[i,:] = xs[i]
        ysExt[i,:] = ys[i]
        for m in range(0,3):
            normalExt[m,i,:] = normal[m,i]
        for k in range(0,nt):
            Aext_u[i,:,k]   = FunctionV1 [i,k]
            Aext_v[i,:,k]   = FunctionV2 [i,k]
            Aext_w[i,:,k]   = FunctionV3 [i,k]
            Aext_p[i,:,k]   = FunctionP  [i,k]
            Aext_rho[i,:,k] = FunctionRho[i,k]
    if nzExtr == 1:
        elesurfZ = np.ones(nzExtr)
    else:
        zsExt[0] = zs-Lhalfz  # SHOULD I USE NP.COPY?
        dzExt = 2.*Lhalfz / (nzExtr-1)
        for l in range(1,nzExtr):
            if (l==1):
                print('  Extrusion 2D to 3D: steps remaining ' + str(nxs+nzExtr-2+nzExtr-1-1) + ' of ' + str(nxs+nzExtr-1+nzExtr-2+nxs))            
            zsExt[l] = zsExt[l-1] + dzExt

        elesurfZ     = np.zeros(nzExtr)
        elesurfZ[0]  = 0.5 * (zsExt[1]  - zsExt[0] )
        elesurfZ[-1] = 0.5 * (zsExt[-1] - zsExt[-2])
        for i in range(1,nzExtr-1):
            if (i==1):
                print('  Extrusion 2D to 3D: steps remaining ' + str(nxs+nzExtr-2-1) + ' of ' + str(nxs+nzExtr-1+nzExtr-2+nxs))            
            elesurfZ[i] = 0.5 * (zsExt[i+1]-zsExt[i-1])

    for i in range(0,nxs):
        if (i==0):
            print('  Extrusion 2D to 3D: steps remaining ' + str(nxs-1) + ' of ' + str(nxs+nzExtr-1+nzExtr-2+nxs))
        for l in range(0,nzExtr):
            elesurfExt[i,l] = elesurfXY[i] * elesurfZ[l]

    return Aext_u,Aext_v,Aext_w,Aext_p,nzs,Aext_rho,xsExt,ysExt,zsExt,normalExt,elesurfExt
