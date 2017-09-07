# ****************************************************************************
#
# Advanced Time approach for pressure reconstruction of FW-H method.
#
# Autor: Jaime Vaquero
# June 2017
#
# ****************************************************************************
import numpy as np
#from math import *
#from cmath import *
from Extrusion2Dto3D import Extrusion2Dto3D
#from Formulation1A import Formulation1A
from Formulation1A_2 import Formulation1A  # Optimized version
from Margnat2010 import Margnat2010
import matplotlib.pyplot as plt



# ntf  = 77#nt of file
# nxs = 600
# NT  = 1#reproducing the file periodically
# nt = ntf*NT

#ntf = 123#nt of file SCCM
ntf = 481#nt of file ChX
#ntf = 180#nt of file ChX
NT  = 1#reproducing the file periodically
nt  = ntf*NT

#nxs = 600 # SCCM
nxs = 1480 # CharlesX
ScaleSurf = 15
if (nxs % ScaleSurf == 0):
    nxsFWH = int(nxs/ScaleSurf)
else:
    nxsFWH = int(nxs/ScaleSurf)+1

t = np.zeros(nt)

zs      = 0.
xs      = np.zeros(nxsFWH)
ys      = np.zeros(nxsFWH)
rs      = np.zeros(nxsFWH)
ths     = np.zeros(nxsFWH)
elesurf = np.zeros(nxsFWH)
drs     = np.zeros([3,nxsFWH])
kvector = np.zeros([3,nxsFWH])
normal  = np.zeros([3,nxsFWH])

ordIndex = np.zeros(nxs)

R_xs  = np.zeros(nxs)
R_ys  = np.zeros(nxs)
R_ux  = np.zeros([nxs,nt])
R_uy  = np.zeros([nxs,nt])
R_p   = np.zeros([nxs,nt])
R_rho = np.zeros([nxs,nt])

O_rs  = np.zeros(nxs)
O_ths = np.zeros(nxs)
O_xs  = np.zeros(nxs)
O_ys  = np.zeros(nxs)
O_ux  = np.zeros([nxs,nt])
O_uy  = np.zeros([nxs,nt])
O_p   = np.zeros([nxs,nt])
O_rho = np.zeros([nxs,nt])

A_ux  = np.zeros([nxsFWH,nt])
A_uy  = np.zeros([nxsFWH,nt])
A_uz  = np.zeros([nxsFWH,nt])
A_p   = np.zeros([nxsFWH,nt])
A_rho = np.zeros([nxsFWH,nt])




# tobs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])
# robs = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])
# pobs  = np.zeros([LengthSpaceX,LengthSpaceExtZ,LengthTime])


# OBSERVER
xo = 0.0
yo = 200.*0.019
zo = 0.0

# Problem PARAMETERS
p0   = 1.e5
#rho0 = 1.225 # SCCM
#rho0 = 1.1476 # CharlesX only diff mu
rho0 = 1.211 # Charles Good parameters
#c0   = 338.1 #SCCM
c0   = 340. # CharlesX Good parameters
Mach = 0.2
#dt   = 1e-4 # SCCM
#dt   = 5.17e-5 # CharlesX only diff mu
dt   = 5.11e-5 # CharlesX Good parameters

#Incude date and time
from datetime import datetime
print('\n')
print('Nous sommes le ' + datetime.now().strftime('%d-%m-%Y') + '.')
print('Il est ' + datetime.now().strftime('%H:%M') + ' et '
                + datetime.now().strftime('%S.%f') + ' secondes.')
print('\n')

print('***********************************************')
print('*  ADVANCED TIME APPROACH : FWH WEM           *')
print('*  MACH                   :',Mach,'              *')
print('***********************************************')

print('')
for kk in range(0,NT):
    for k in range(0,ntf):
        if(float(k+1+kk*ntf)/float(NT*ntf)*100. % 10 == 0 or float(k+1+kk*ntf)/float(NT*ntf)*100. % 10 < 0.2):
            print('  Reading data: ' + str(round(float(k+1+kk*ntf)/float(NT*ntf)*100)) + ' %')
        t[k]  = float(k+1+kk*ntf)*dt

#        file201 = np.loadtxt(open('../../CYLINDRE_W/BDD/13-06-SCCM/XYZ_Internal_Table_C0019_Re_150_cylinderwall_DT1e-4_' + str(k) + '.csv', 'rb'),delimiter = ',', skiprows = 1) # SCCM Stru 13-06
#        file201 = np.loadtxt(open('../../CYLINDRE_W/BDD/C0019_Re150_Mlower02_onlyDiffMU/C0019_Re150_Mlower02_onlydiffmu_c350_DTparaview5-17e-5.' + str(k) + '.csv', 'rb'),delimiter = ',', skiprows = 1) # ChX
        file201 = np.loadtxt(open('../../CYLINDRE_W/BDD/C0019_Re150_M02_ChX_BONpar/C0019_Re150_M02_ChX_BONpar/C0019_Re150_M02_ChX_bonspars_DTparav5-11e-5_1R.' + str(k) + '.csv', 'rb'),delimiter = ',', skiprows = 1) # CharlesX Good Parameters

#        R_xs[:]    = file201[:,5] #SCCM stru
#        R_ys[:]    = file201[:,6] # Maybe useless!
#        R_ux [:,k] = file201[:,3]
#        R_uy [:,k] = file201[:,4]
#        R_p  [:,k] = file201[:,0]
#        R_rho[:,k] = file201[:,1]

        R_xs[:]    = file201[:,5] #ChX
        R_ys[:]    = file201[:,6] # Maybe useless!
        R_ux [:,k] = file201[:,2]
        R_uy [:,k] = file201[:,3]
        R_p  [:,k] = file201[:,1]
        R_rho[:,k] = file201[:,0]

        O_rs[:]       = np.sqrt(R_xs[:]*R_xs[:] + R_ys[:]*R_ys[:])
        O_ths[:]      = np.arctan2(R_ys[:],R_xs[:])
        ordIndex[:]   = np.argsort(O_ths[:]) # reorganizes the indexes from lowest ths to highest ths


        iss = -1
        for m in ordIndex:
            iss = iss +1
            O_xs[iss]    = R_xs[int(m)]
            O_ys[iss]    = R_ys[int(m)]
            O_ux[iss,k + kk*ntf]  = R_ux[int(m),k]
            O_uy[iss,k + kk*ntf]  = R_uy[int(m),k]
            O_p[iss,k + kk*ntf]   = R_p[int(m),k]
            O_rho[iss,k + kk*ntf] = R_rho[int(m),k]

O_rs  = np.sqrt(O_xs*O_xs[:] + O_ys[:]*O_ys[:])
O_ths = np.arctan2(O_ys,O_xs)


iss = -1
for isc in range(0,nxs,ScaleSurf):
    iss = iss + 1
    xs[iss]      = O_xs[isc]
    ys[iss]      = O_ys[isc]
    A_ux[iss,:]  = O_ux[isc,:]
    A_uy[iss,:]  = O_uy[isc,:]
    A_p[iss,:]   = O_p[isc,:]
    A_rho[iss,:] = O_rho[isc,:]
    rs[iss]      = O_rs[isc]
    ths[iss]     = O_ths[isc]

# Tangential vector ----------------
drs[0,0]   = 0.5 * (xs[1] - xs[-1])
drs[1,0]   = 0.5 * (ys[1] - ys[-1])
drs[0,-1]  = 0.5 * (xs[0] - xs[-2])
drs[1,-1]  = 0.5 * (ys[0] - ys[-2])
for i in range(1,nxsFWH-1):
    drs[0,i] = 0.5 * (xs[i+1] - xs[i-1])
    drs[1,i] = 0.5 * (ys[i+1] - ys[i-1])


# Surface (line) element and normal vector. I must be sure of the direction
# of variaton of dths in order to have the normal pointing outwards
if (ths[0] - ths[-1] > 0):
    kvector[2,0] = 1.
else:
    kvector[2,0] = -1.
for iss in range(1,nxsFWH):
    if (ths[iss] - ths[iss-1] > 0):
        kvector[2,iss] = 1.
    else:
        kvector[2,iss] = -1.

for i in range(0,nxsFWH):
    elesurf[i]  = np.sqrt(drs[0,i]**2. + drs[1,i]**2.)
    normal[:,i] = np.cross(drs[:,i],kvector[:,i]) / np.linalg.norm(np.cross(drs[:,i],kvector[:,i]))

# Extrusion of data from 2D to 3D in z direction
Lhz = 10.*27.3*0.019 # It is half the extrusion length! Extrusion is symmetric
ScaleNzs = 100.
#nzs = int((2.*Lhz/(min(elesurf)))/ScaleNzs)+1
nzs = 23
print('')
print('  Data reading     :  OK')
print('  Extrusion points : ' , nzs)
print('  Total points     : ',nxsFWH*nzs)
print('')


Aext_ux,Aext_uy,Aext_uz,Aext_p,nzs,Aext_rho,xsExt,ysExt,zsExt,normalExt,elesurfExt = \
    Extrusion2Dto3D(nxsFWH,nt,A_ux,A_uy,A_uz,A_p,A_rho,xs,ys,zs,normal,elesurf,nzs,Lhz)

print('\n')
if nzs == 1:
    print('  Extrusion        : not required')
else:
    print('  Extrusion        : OK')

# FARASSAT FORMULATION

pobs,tobs,robs = Formulation1A(nt,nxsFWH,nzs,t,xsExt,ysExt,zsExt,normalExt,elesurfExt,Aext_ux,Aext_uy, \
                  Aext_uz,Aext_p,Aext_rho,np.array([xo,yo,zo]),dt,c0,p0,rho0)
print('')                
print('  Formulation1A    : OK')



nlM        = nt+1+int(robs.max()/(c0*dt))
tadvmin    = (1 + int(robs.min() / (c0*dt)))*dt
tadvmax    = nt-1+int(robs.max()/(c0*dt))

print('')
print('  First signal received at: ' , (1 + int(robs.min() / (c0*dt)) )*dt*Mach*c0/0.019)
print('  First valid time        : ' , (2 + int(robs.max() / (c0*dt)) )*dt*Mach*c0/0.019)
print('  Last valid time         : ' , (nt +int(robs.min() / (c0*dt)))*dt*Mach*c0/0.019)
print('')


pobs_discret = np.zeros([nxsFWH,nzs,nlM])
pac  = np.zeros(nlM)
tadv = np.zeros(nlM)



#################################### MARGNAT ##################################
pobs_discret,tadv = Margnat2010(nt,nlM,nxsFWH,nzs,robs,pobs,c0,dt,pobs_discret,tadv)
print('  Margnat 2010 Algorithm   : OK')
print('')


################################# CASALINO ####################################
#wij = np.zeros([nxsFWH,nzs,nlM])
#for j in range(0,nxsFWH):
#    for l in range(0,nzs):
#        for k in range(0,nt):
#            kadv = int(tobs[j,l,k]/dt)
#            w    = tobs[j,l,k]/dt-kadv
#            if pobs_discret[j,l,kadv] == 0:
#                pobs_discret[j,l,kadv] = pobs[j,l,k]   # k or kadv????
#                wij[j,l,kadv] = w # k or kadv????
#            else:
#                pobs_discret[j,l,kadv] = pobs[j,l,k] - (pobs_discret[j,l,kadv]-pobs[j,l,k])/(wij[j,l,kadv]-w)*w
#                wij[j,l,kadv] = 0.0

for k in range(0,nlM):
    pac[k]  = pobs_discret[:,:,k].sum()
            
                

prms  = 0.0
ntadv = nt   + int(robs.min() / (c0*dt))  - (2 + int(robs.max() / (c0*dt))) + 1
kmin  = 2    + int(robs.max() / (c0*dt))
kmax  = nt+1 + int(robs.min() / (c0*dt))

#for k in range(kmin,kmax):
#    prms = prms + (pac[k]-np.mean(pac[kmin:kmax-1]))**2.*dt
#prms = (prms/((ntadv)*dt))**0.5/(rho0*c0**2.)
prms = np.sqrt(np.mean((pac[kmin:kmax-1]-np.mean(pac[kmin:kmax-1]))**2.))/(rho0*c0**2.)
print('  Prms    : ', prms)
print('  Mean Pac: ', np.mean(pac[kmin:kmax-1])/(rho0*c0**2.*Mach**2.5))




################## PLOTTING #################################################
# Pressure square over time
plt.figure()
plt.plot(tadv[kmin:kmax-1]*68./0.019,(pac[kmin:kmax-1]-np.mean(pac[kmin:kmax-1]))**2./(rho0*c0**2.),'k-o')
plt.ylabel('(p-pmean)^2')
plt.xlabel('tU0/D')
plt.grid()


# Pressure over time
filePACCylFWH = open('../../CYLINDRE_W/FWH/pacAdim90rendT_Starccm_2w_8p.dat','r')
pacFromCyl = filePACCylFWH.read().split()
filePACCylFWH.close()


plt.figure()
plt.plot(tadv*68./0.019,(pac-np.mean(pac[kmin:kmax-1]))/(rho0*c0**2.*Mach**2.5),'r-')
plt.plot(np.linspace(0,36*0.35789,len(pacFromCyl))+70,pacFromCyl,'k-')
plt.xlabel('tU0/D')
plt.ylabel('pac')
plt.grid()
plt.show()


print('\n')
print('Nous sommes le ' + datetime.now().strftime('%d-%m-%Y') + '.')
print('Il est ' + datetime.now().strftime('%H:%M') + ' et '
                + datetime.now().strftime('%S.%f') + ' secondes.')
print('\n')
