#
#**********************************************************************
# Programme d'extrapolation d'ondes acoustiques
# Formalisme de Ffwocs-Williams & Hawkings
# Densité moyenne
# RS - 30/07/2014
#
# Monopole solution analytique These de Xavier Gloerfelt p.216
# Adaptee 3D
# RS - 26/05/2016
# Ecrit sur python
# JV 01/02017
#**********************************************************************
#

import numpy as np
from math import *
from cmath import *
from scipy.io import FortranFile
from scipy.special import hankel2

# Variables definition
ii   = 1j
nxo  = 201; nyo = nxo; nzo = nxo # 201
T    = 100.          # 100
nper = 4
nt   = floor(nper * T) + 728 - 300
nxs  = 81; nys = nxs; nzs = nxs  # 81
ds   = (nxs - 1) / 2

xo = np.zeros(nxo)
yo = np.zeros(nyo)
zo = np.zeros(nzo)

t      = np.zeros(nt)
k_onde = np.zeros(nt)
omega  = np.zeros(nt)

D2px = np.zeros([nxo,nzo])
D2ux = np.zeros([nxo,nzo])
D2vx = np.zeros([nxo,nzo])
D2wx = np.zeros([nxo,nzo])

phi      = np.zeros(nt, dtype = np.complex)
dphisdt  = np.zeros(nt, dtype = np.complex)
dphisdxo = np.zeros(nt, dtype = np.complex)
dphisdyo = np.zeros(nt, dtype = np.complex)
dphisdzo = np.zeros(nt, dtype = np.complex)
pprim    = np.zeros(nt, dtype = np.complex)


# D2pxt = np.zeros([nxo,nzo,nt])
# D2uxt = np.zeros([nxo,nzo,nt])
# D2vxt = np.zeros([nxo,nzo,nt])
# D2wxt = np.zeros([nxo,nzo,nt])

D2pxt = np.zeros([nxs*2 + nys*2,nt])
D2uxt = np.zeros([nxs*2 + nys*2,nt])
D2vxt = np.zeros([nxs*2 + nys*2,nt])
D2wxt = np.zeros([nxs*2 + nys*2,nt])



D3px = np.zeros([nxs,nys,nzs])
D3ux = np.zeros([nxs,nys,nzs])
D3vx = np.zeros([nxs,nys,nzs])
D3wx = np.zeros([nxs,nys,nzs])

D3pxt = np.zeros([nxs,nys,nzs,nt])
D3uxt = np.zeros([nxs,nys,nzs,nt])
D3vxt = np.zeros([nxs,nys,nzs,nt])
D3wxt = np.zeros([nxs,nys,nzs,nt])

#Include date and time
from datetime import datetime
print('\n')
print('Nous sommes le ' + datetime.now().strftime('%d-%m-%Y') + '.')
print('Il est ' + datetime.now().strftime('%H:%M') + ' et '
                + datetime.now().strftime('%S.%f') + ' secondes.')
print('\n')

# ============================================================================

def numcar3(num):
    if num >= 100:
        car = str(num)
    elif num >= 10:
        car = '0'  + str(num)
    else:
        car = '00' + str(num)

    return car


# PARAMETERS ==================================================================
Amplitude = 5.
lac       = 10.

Mach = 0.
rho0 = 1.22
p0   = 1.0e5
co   = 340.

b2   = 1. - Mach**2.
fact = - Amplitude / (4. * np.pi)

dt = lac / (100. * co)

# SOURCE ======================================================================
Lxo = 100.
Lyo = Lxo
Lzo = Lxo

dxo = Lxo / (nxo - 1)
xo[0] = -0.
for io in range(1, nxo):
    xo[io] = xo[io-1] + dxo

dyo = Lyo / (nyo - 1)
yo[0] = -0.
for jo in range(1, nyo):
    yo[jo] = yo[jo-1] + dyo

dzo = Lzo / (nzo - 1)
zo[0] = -50.
for lo in range(1, nxo):
    zo[lo] = zo[lo-1] + dzo

# POSITION DE LA SOURCE, on enlève 1 indice pour que ce soit
# cohérent avec le code fortran
XS  = 50.
YS  = 50.
ZS  = 0.
ism = floor((nxo+1)/2 - ds - 1)
isp = floor((nxo+1)/2 + ds - 1)
jsm = floor((nyo+1)/2 - ds - 1)
jsp = floor((nyo+1)/2 + ds - 1)
lsm = floor((nzo+1)/2 - ds - 1)
lsp = floor((nzo+1)/2 + ds - 1)
print(ism,isp,jsm,jsp,lsm,lsp)

file215 = open('../POSTRAITEMENT_py/xobs.dat', 'w')
file216 = open('../POSTRAITEMENT_py/yobs.dat', 'w')
file217 = open('../POSTRAITEMENT_py/zobs.dat', 'w')

for io in range(0, nxo):
	file215.write(str(xo[io])+'\n')
for jo in range(0, nyo):
	file216.write(str(yo[jo])+'\n')
for lo in range(0, nzo):
	file217.write(str(zo[lo])+'\n')

file215.close()
file216.close()
file217.close()

print(' Observer coordinates writing : OK')




#-------------------------------
for k in range(0,nt):
    t[k]      = float(k+1)*dt # Attention to indexes!
    k_onde[k] = 2*pi/lac
    omega[k]  = k_onde[k]*co



# CIRC---------------------------------------------------------------
# D2pxt_th = np.zeros([nth,nt])
# D2uxt_th = np.zeros([nth,nt])
# D2vxt_th = np.zeros([nth,nt])
# D2wxt_th = np.zeros([nth,nt])
#
# for iss in range(0,nth):
# # for iss in range(0,nxs2):
#     lo = int((nyo+1)/2) - 1  # Attention to indexes!
#     # for io in range(0,nxo):
#     # r1 = (xo2[iss] - XS)
#     # r2 = (yo2[iss] - YS)
#     # r3 = (zo2     - ZS)
#     r1 = (xsth[iss] - XS + XS)
#     r2 = (ysth[iss] - YS + YS)
#     r3 = (zo2     - ZS)
#     xb = np.sqrt( r1**2. + b2*r2**2. + b2*r3**2. )
#     if xb == 0:
#         phi      = 0
#         dphisdt  = 0
#         dphisdxo = 0
#         dphisdyo = 0
#         dphisdzo = 0
#     else:
#         # phi      = fact * \
#         #             np.exp( ii*omega*t - (ii*k_onde/b2)*(xb-Mach*r1))/xb
#         # dphisdt  = ii*omega*phi
#         # dphisdxo = fact*( np.exp( ii*omega*t)/xb**2)* \
#         #                 (-xb*ii*k_onde/b2 * ((r1/xb)-Mach) - (r1/xb))* \
#         #                 np.exp(-(ii*k_onde/b2)*(xb-Mach*r1))
#         # dphisdyo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
#         #                 (r2/xb)  - b2 * (r2/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
#         # dphisdzo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
#         #                 (r3/xb)  - b2 * (r3/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
#
#         phi      = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb)
#         dphisdt  = ii*omega*phi
#         dphisdxo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb) *ii*Mach*k_onde/b2 \
#                              + Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*\
#                              k_onde/b2 * r1/xb
#         dphisdyo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*k_onde * r2/xb
#         dphiszo = 0.0
#
#
#     Re = ( dphisdt + (Mach*co)*dphisdxo ).real
#     D2pxt_th[iss,:] = -rho0 * Re + p0
#     D2uxt_th[iss,:] = Mach*co + dphisdxo.real
#     D2vxt_th[iss,:] = dphisdyo.real
#     D2wxt_th[iss,:] = dphisdzo.real
#     if iss== 42:
#         print(-rho0 * Re + p0)
#         print(dphisdxo.real)
#         print(dphisdyo.real)
#
#
# prms = np.zeros(nth)
# for iss in range(0,nth):
#     sum1 = 0.
#     for k in range(0,nt):
#         sum1 = sum1 + (D2pxt_th[iss,k]-p0)**2.*dt
#     prms[iss] = (sum1 / (nt*dt))**0.5
#
# import matplotlib.pyplot as plt
# plt.figure
# plt.plot(prms*np.cos(theta),prms*np.sin(theta),'o')
# plt.show()
#
#
# # D2pxt_bis = np.zeros(nxs2*nt)
# # D2uxt_bis = np.zeros(nxs2*nt)
# # D2vxt_bis = np.zeros(nxs2*nt)
# # D2wxt_bis = np.zeros(nxs2*nt)
# D2pxt_bis_th = np.zeros(nth*nt)
# D2uxt_bis_th = np.zeros(nth*nt)
# D2vxt_bis_th = np.zeros(nth*nt)
# D2wxt_bis_th = np.zeros(nth*nt)
# for k in range(0,nt):
#     for iss in range(0,nth):
#     # for iss in range(0,nxs2):
#         # D2pxt_bis[iss + k*nxs2] = D2pxt[iss,k]
#         # D2uxt_bis[iss + k*nxs2] = D2uxt[iss,k]
#         # D2vxt_bis[iss + k*nxs2] = D2vxt[iss,k]
#         # D2wxt_bis[iss + k*nxs2] = D2wxt[iss,k]
#         D2pxt_bis_th[iss + k*nth] = D2pxt_th[iss,k]
#         D2uxt_bis_th[iss + k*nth] = D2uxt_th[iss,k]
#         D2vxt_bis_th[iss + k*nth] = D2vxt_th[iss,k]
#         D2wxt_bis_th[iss + k*nth] = D2wxt_th[iss,k]
#
# file1007 = open('../BDD_py/pp.ct.dat', 'w')
# file1008 = open('../BDD_py/ux.ct.dat', 'w')
# file1009 = open('../BDD_py/uy.ct.dat', 'w')
# file1010 = open('../BDD_py/uz.ct.dat', 'w')
# # for k in range(0,nt*nxs2):
# for k in range(0,nt*nth):
#     file1007.write(str(D2pxt_bis_th[k])+'\n')
#     file1008.write(str(D2uxt_bis_th[k])+'\n')
#     file1009.write(str(D2vxt_bis_th[k])+'\n')
#     file1010.write(str(D2wxt_bis_th[k])+'\n')
# file1007.close()
# file1008.close()
# file1009.close()
# file1010.close()
# print(' CIRC  surface writing:  OK')


#TOP---------------------------------------------------------------
iss = -1
for io in range(ism,isp+1):
    iss = iss + 1
    lo = int((nyo+1)/2) - 1  # Attention to indexes!
    # for io in range(0,nxo):
    r1 = (xo[io]  - XS)
    r2 = (yo[jsp] - YS)
    r3 = (0.0     - ZS)
    xb = np.sqrt( r1**2. + b2*r2**2. + b2*r3**2. )
    if xb == 0.:
        phi      = 0.
        dphisdt  = 0.
        dphisdxo = 0.
        dphisdyo = 0.
        dphisdzo = 0.
    else:
        # phi      = fact * \
        #             np.exp( ii*omega*t - (ii*k_onde/b2)*(xb-Mach*r1))/xb
        # dphisdt  = ii*omega*phi
        # dphisdxo = fact*( np.exp( ii*omega*t)/xb**2)* \
        #                 (-xb*ii*k_onde/b2 * ((r1/xb)-Mach) - (r1/xb))* \
        #                 np.exp(-(ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdyo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r2/xb)  - b2 * (r2/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdzo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r3/xb)  - b2 * (r3/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))

        phi      = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb)
        dphisdt  = ii*omega*phi
        dphisdxo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb) *ii*Mach*k_onde/b2 \
                             + Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*\
                             k_onde/b2 * r1/xb
        dphisdyo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*k_onde * r2/xb
        dphiszo = 0.0


    Re = ( dphisdt + (Mach*co)*dphisdxo ).real
    D2pxt[iss,:] = -rho0 * Re + p0
    D2uxt[iss,:] = Mach*co + dphisdxo.real
    D2vxt[iss,:] = dphisdyo.real
    D2wxt[iss,:] = dphisdzo.real

# D2pxt_bis = np.zeros(nxs2*nt)
# D2uxt_bis = np.zeros(nxs2*nt)
# D2vxt_bis = np.zeros(nxs2*nt)
# D2wxt_bis = np.zeros(nxs2*nt)
D2pxt_bis = np.zeros(nxs*nt)
D2uxt_bis = np.zeros(nxs*nt)
D2vxt_bis = np.zeros(nxs*nt)
D2wxt_bis = np.zeros(nxs*nt)
for k in range(0,nt):
    for iss in range(0,nxs):
    # for iss in range(0,nxs2):
        # D2pxt_bis[iss + k*nxs2] = D2pxt[iss,k]
        # D2uxt_bis[iss + k*nxs2] = D2uxt[iss,k]
        # D2vxt_bis[iss + k*nxs2] = D2vxt[iss,k]
        # D2wxt_bis[iss + k*nxs2] = D2wxt[iss,k]
        D2pxt_bis[iss + k*nxs] = D2pxt[iss,k]
        D2uxt_bis[iss + k*nxs] = D2uxt[iss,k]
        D2vxt_bis[iss + k*nxs] = D2vxt[iss,k]
        D2wxt_bis[iss + k*nxs] = D2wxt[iss,k]

file1007 = open('../BDD_py/pp.xtp.dat', 'w')
file1008 = open('../BDD_py/ux.xtp.dat', 'w')
file1009 = open('../BDD_py/uy.xtp.dat', 'w')
file1010 = open('../BDD_py/uz.xtp.dat', 'w')
# for k in range(0,nt*nxs2):
for k in range(0,nt*nxs):
    file1007.write(str(D2pxt_bis[k])+'\n')
    file1008.write(str(D2uxt_bis[k])+'\n')
    file1009.write(str(D2vxt_bis[k])+'\n')
    file1010.write(str(D2wxt_bis[k])+'\n')
file1007.close()
file1008.close()
file1009.close()
file1010.close()
print(' TOP    surface writing:  OK')

# BOTTOM ------------------------------------------------------
iss = -1
for io in range(ism,isp+1):
    iss = iss + 1
    lo = int((nyo+1)/2) - 1  # Attention to indexes!
    # for io in range(0,nxo):
    r1 = (xo[io]  - XS)
    r2 = (yo[jsp] - YS)
    r3 = (0.     - ZS)
    xb = np.sqrt( r1**2. + b2*r2**2. + b2*r3**2. )
    if xb == 0.:
        phi      = 0.
        dphisdt  = 0.
        dphisdxo = 0.
        dphisdyo = 0.
        dphisdzo = 0.
    else:
        # phi      = fact * \
        #             np.exp( ii*omega*t - (ii*k_onde/b2)*(xb-Mach*r1))/xb
        # dphisdt  = ii*omega*phi
        # dphisdxo = fact*( np.exp( ii*omega*t)/xb**2)* \
        #                 (-xb*ii*k_onde/b2 * ((r1/xb)-Mach) - (r1/xb))* \
        #                 np.exp(-(ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdyo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r2/xb)  - b2 * (r2/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdzo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r3/xb)  - b2 * (r3/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))

        phi      = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb)
        dphisdt  = ii*omega*phi
        dphisdxo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb) *ii*Mach*k_onde/b2 \
                             + Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*\
                             k_onde/b2 * r1/xb
        dphisdyo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*k_onde * r2/xb
        dphiszo = 0.0


    Re = ( dphisdt + (Mach*co)*dphisdxo ).real
    D2pxt[iss,:] = -rho0 * Re + p0
    D2uxt[iss,:] = Mach*co + dphisdxo.real
    D2vxt[iss,:] = dphisdyo.real
    D2wxt[iss,:] = dphisdzo.real

# D2pxt_bis = np.zeros(nxs2*nt)
# D2uxt_bis = np.zeros(nxs2*nt)
# D2vxt_bis = np.zeros(nxs2*nt)
# D2wxt_bis = np.zeros(nxs2*nt)
D2pxt_bis = np.zeros(nxs*nt)
D2uxt_bis = np.zeros(nxs*nt)
D2vxt_bis = np.zeros(nxs*nt)
D2wxt_bis = np.zeros(nxs*nt)
for k in range(0,nt):
    for iss in range(0,nxs):
    # for iss in range(0,nxs2):
        # D2pxt_bis[iss + k*nxs2] = D2pxt[iss,k]
        # D2uxt_bis[iss + k*nxs2] = D2uxt[iss,k]
        # D2vxt_bis[iss + k*nxs2] = D2vxt[iss,k]
        # D2wxt_bis[iss + k*nxs2] = D2wxt[iss,k]
        D2pxt_bis[iss + k*nxs] = D2pxt[iss,k]
        D2uxt_bis[iss + k*nxs] = D2uxt[iss,k]
        D2vxt_bis[iss + k*nxs] = D2vxt[iss,k]
        D2wxt_bis[iss + k*nxs] = D2wxt[iss,k]


file1007 = open('../BDD_py/pp.xtm.dat', 'w')
file1008 = open('../BDD_py/ux.xtm.dat', 'w')
file1009 = open('../BDD_py/uy.xtm.dat', 'w')
file1010 = open('../BDD_py/uz.xtm.dat', 'w')
# for k in range(0,nt*nxs2):
for k in range(0,nt*nxs):
    file1007.write(str(D2pxt_bis[k])+'\n')
    file1008.write(str(D2uxt_bis[k])+'\n')
    file1009.write(str(D2vxt_bis[k])+'\n')
    file1010.write(str(D2wxt_bis[k])+'\n')
file1007.close()
file1008.close()
file1009.close()
file1010.close()
print(' BOTTOM surface writing:  OK')


# RIGHT ------------------------------------------------------
iss = -1
for jo in range(jsm,jsp+1):
    iss = iss + 1
    lo = int((nyo+1)/2) - 1  # Attention to indexes!
    # for io in range(0,nxo):
    r1 = (xo[isp]  - XS)
    r2 = (yo[jo] - YS)
    r3 = (0.     - ZS)
    xb = np.sqrt( r1**2. + b2*r2**2. + b2*r3**2. )
    if xb == 0.:
        phi      = 0.
        dphisdt  = 0.
        dphisdxo = 0.
        dphisdyo = 0.
        dphisdzo = 0.
    else:
        # phi      = fact * \
        #             np.exp( ii*omega*t - (ii*k_onde/b2)*(xb-Mach*r1))/xb
        # dphisdt  = ii*omega*phi
        # dphisdxo = fact*( np.exp( ii*omega*t)/xb**2)* \
        #                 (-xb*ii*k_onde/b2 * ((r1/xb)-Mach) - (r1/xb))* \
        #                 np.exp(-(ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdyo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r2/xb)  - b2 * (r2/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdzo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r3/xb)  - b2 * (r3/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))

        phi      = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb)
        dphisdt  = ii*omega*phi
        dphisdxo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb) *ii*Mach*k_onde/b2 \
                             + Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*\
                             k_onde/b2 * r1/xb
        dphisdyo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*k_onde * r2/xb
        dphiszo = 0.0


    Re = ( dphisdt + (Mach*co)*dphisdxo ).real
    D2pxt[iss,:] = -rho0 * Re + p0
    D2uxt[iss,:] = Mach*co + dphisdxo.real
    D2vxt[iss,:] = dphisdyo.real
    D2wxt[iss,:] = dphisdzo.real

# D2pxt_bis = np.zeros(nxs2*nt)
# D2uxt_bis = np.zeros(nxs2*nt)
# D2vxt_bis = np.zeros(nxs2*nt)
# D2wxt_bis = np.zeros(nxs2*nt)
D2pxt_bis = np.zeros(nys*nt)
D2uxt_bis = np.zeros(nys*nt)
D2vxt_bis = np.zeros(nys*nt)
D2wxt_bis = np.zeros(nys*nt)
for k in range(0,nt):
    for iss in range(0,nys):
    # for iss in range(0,nxs2):
        # D2pxt_bis[iss + k*nxs2] = D2pxt[iss,k]
        # D2uxt_bis[iss + k*nxs2] = D2uxt[iss,k]
        # D2vxt_bis[iss + k*nxs2] = D2vxt[iss,k]
        # D2wxt_bis[iss + k*nxs2] = D2wxt[iss,k]
        D2pxt_bis[iss + k*nys] = D2pxt[iss,k]
        D2uxt_bis[iss + k*nys] = D2uxt[iss,k]
        D2vxt_bis[iss + k*nys] = D2vxt[iss,k]
        D2wxt_bis[iss + k*nys] = D2wxt[iss,k]


file1007 = open('../BDD_py/pp.ytp.dat', 'w')
file1008 = open('../BDD_py/ux.ytp.dat', 'w')
file1009 = open('../BDD_py/uy.ytp.dat', 'w')
file1010 = open('../BDD_py/uz.ytp.dat', 'w')
# for k in range(0,nt*nxs2):
for k in range(0,nt*nys):
    file1007.write(str(D2pxt_bis[k])+'\n')
    file1008.write(str(D2uxt_bis[k])+'\n')
    file1009.write(str(D2vxt_bis[k])+'\n')
    file1010.write(str(D2wxt_bis[k])+'\n')
file1007.close()
file1008.close()
file1009.close()
file1010.close()
print(' RIGHT  surface writing:  OK')

# LEFT ------------------------------------------------------
iss = -1
for jo in range(jsm,jsp+1):
    iss = iss + 1
    lo = int((nyo+1)/2) - 1  # Attention to indexes!
    # for io in range(0,nxo):
    r1 = (xo[ism] - XS)
    r2 = (yo[jo]  - YS)
    r3 = (0.     - ZS)
    xb = np.sqrt( r1**2. + b2*r2**2. + b2*r3**2. )
    if xb == 0.:
        phi      = 0.
        dphisdt  = 0.
        dphisdxo = 0.
        dphisdyo = 0.
        dphisdzo = 0.
    else:
        # phi      = fact * \
        #             np.exp( ii*omega*t - (ii*k_onde/b2)*(xb-Mach*r1))/xb
        # dphisdt  = ii*omega*phi
        # dphisdxo = fact*( np.exp( ii*omega*t)/xb**2)* \
        #                 (-xb*ii*k_onde/b2 * ((r1/xb)-Mach) - (r1/xb))* \
        #                 np.exp(-(ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdyo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r2/xb)  - b2 * (r2/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
        # dphisdzo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
        #                 (r3/xb)  - b2 * (r3/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))

        phi      = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb)
        dphisdt  = ii*omega*phi
        dphisdxo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb) *ii*Mach*k_onde/b2 \
                             + Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*\
                             k_onde/b2 * r1/xb
        dphisdyo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*k_onde * r2/xb
        dphiszo = 0.0


    Re = ( dphisdt + (Mach*co)*dphisdxo ).real
    D2pxt[iss,:] = -rho0 * Re + p0
    D2uxt[iss,:] = Mach*co + dphisdxo.real
    D2vxt[iss,:] = dphisdyo.real
    D2wxt[iss,:] = dphisdzo.real

# D2pxt_bis = np.zeros(nxs2*nt)
# D2uxt_bis = np.zeros(nxs2*nt)
# D2vxt_bis = np.zeros(nxs2*nt)
# D2wxt_bis = np.zeros(nxs2*nt)
D2pxt_bis = np.zeros(nys*nt)
D2uxt_bis = np.zeros(nys*nt)
D2vxt_bis = np.zeros(nys*nt)
D2wxt_bis = np.zeros(nys*nt)
for k in range(0,nt):
    for iss in range(0,nys):
    # for iss in range(0,nxs2):
        # D2pxt_bis[iss + k*nxs2] = D2pxt[iss,k]
        # D2uxt_bis[iss + k*nxs2] = D2uxt[iss,k]
        # D2vxt_bis[iss + k*nxs2] = D2vxt[iss,k]
        # D2wxt_bis[iss + k*nxs2] = D2wxt[iss,k]
        D2pxt_bis[iss + k*nys] = D2pxt[iss,k]
        D2uxt_bis[iss + k*nys] = D2uxt[iss,k]
        D2vxt_bis[iss + k*nys] = D2vxt[iss,k]
        D2wxt_bis[iss + k*nys] = D2wxt[iss,k]

file1007 = open('../BDD_py/pp.ytm.dat', 'w')
file1008 = open('../BDD_py/ux.ytm.dat', 'w')
file1009 = open('../BDD_py/uy.ytm.dat', 'w')
file1010 = open('../BDD_py/uz.ytm.dat', 'w')
# for k in range(0,nt*nxs2):
for k in range(0,nt*nys):
    file1007.write(str(D2pxt_bis[k])+'\n')
    file1008.write(str(D2uxt_bis[k])+'\n')
    file1009.write(str(D2vxt_bis[k])+'\n')
    file1010.write(str(D2wxt_bis[k])+'\n')
file1007.close()
file1008.close()
file1009.close()
file1010.close()
print(' LEFT   surface writing:  OK')


file3001 = open('../BDD_py/xs.dat','w')
file3002 = open('../BDD_py/ys.dat','w')
# for io in range(0,nxs2):
for io in range(ism,isp+1):
    file3001.write(str(xo[io])+'\n')
for jo in range(jsm,jsp+1):
    file3002.write(str(yo[jo])+'\n')
file3001.close()
file3002.close()
file3003 = open('../BDD_py/zs.dat','w')
for lo in range(lsm,lsp + 1):
    file3003.write(str(zo[lo])+'\n')
file3003.close()

# Pression analytique


D2pxt_an = np.zeros([nxo,nyo,nt])
D2uxt_an = np.zeros([nxo,nyo,nt])
D2vxt_an = np.zeros([nxo,nyo,nt])
D2wxt_an = np.zeros([nxo,nyo,nt])

for jo in range(0, nyo):
    print(' Fluid variables 2D: ' +str(jo + 1) + ' of ' + str(nyo))
    lo = int((nyo+1)/2) - 1  # Attention to indexes!
    for io in range(0,nxo):
        r1 = (xo[io] - XS)
        r2 = (yo[jo] - YS)
        r3 = (zo[lo] - ZS)
        xb = np.sqrt( r1**2 + b2*r2**2 + b2*r3**2 )
        if xb == 0.:
            phi      = 0.
            dphisdt  = 0.
            dphisdxo = 0.
            dphisdyo = 0.
            dphisdzo = 0.
        else:
            # phi      = fact * \
            #             np.exp( ii*omega*t - (ii*k_onde/b2)*(xb-Mach*r1))/xb
            # dphisdt  = ii*omega*phi
            # dphisdxo = fact*( np.exp( ii*omega*t)/xb**2)* \
            #                 (-xb*ii*k_onde/b2 * ((r1/xb)-Mach) - (r1/xb))* \
            #                 np.exp(-(ii*k_onde/b2)*(xb-Mach*r1))
            # dphisdyo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
            #                 (r2/xb)  - b2 * (r2/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))
            # dphisdzo = fact*(np.exp( + ii*omega*t)/xb**2.)*(-xb*ii*k_onde * \
            #                 (r3/xb)  - b2 * (r3/xb))*np.exp(- (ii*k_onde/b2)*(xb-Mach*r1))

            phi      = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb)
            dphisdt  = ii*omega*phi
            dphisdxo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * hankel2(0,k_onde/b2*xb) *ii*Mach*k_onde/b2 \
                                 + Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*\
                                 k_onde/b2 * r1/xb
            dphisdyo = Amplitude * ii / (4.*b2**0.5) * np.exp(ii*omega*t + ii*Mach*k_onde*r1/b2) * (- hankel2(1,k_onde/b2*xb))*k_onde * r2/xb
            dphiszo = 0.0

        Re = ( dphisdt + (Mach*co)*dphisdxo ).real
        D2pxt_an[io,jo,:] = -rho0 * Re + p0
        D2uxt_an[io,jo,:] = Mach*co + dphisdxo.real
        D2vxt_an[io,jo,:] = dphisdyo.real
        D2wxt_an[io,jo,:] = dphisdzo.real

for k in range(0,nt):
    file1007 = FortranFile('../ANALYTIQUE_py/pac.T' + str(numcar3(k+1)) + '.bin', 'w') # '+1' to be coherent with fortran code
    file1007.write_record(np.transpose(D2pxt_an[:,:,k]))
    file1007.close()



print(' FLUID VARIABLES 2D         : OK')



print('\n')
print('Nous sommes le ' + datetime.now().strftime('%d-%m-%Y') + '.')
print('Il est ' + datetime.now().strftime('%H:%M') + ' et '
                + datetime.now().strftime('%S.%f') + ' secondes.')
print('\n')
