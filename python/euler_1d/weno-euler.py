""" WENO scheme for the Euler 1D equation
The flux implemented is the ECUSP flux
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
import matplotlib.pyplot as plt

from utils import P_from_Ev, rhoE_from_Pv
from computeECUSPFlux import compute_ecusp_flux
from computeLFFlux import compute_lf_flux
from defineCase import defineCase

# Initial conditions
def f_0(U):
    b = (zvec < 0.5)
    U[:,0][b] = left[0]
    U[:,1][b] = left[1] * left[0]
    U[:,2][b] = rhoE_from_Pv(left[2],U[:,0][b],U[:,1][b]/U[:,0][b])
    
    b = (zvec>=0.5)
    U[:,0][b] = right[0]
    U[:,1][b] = right[1] * right[0]
    U[:,2][b] = rhoE_from_Pv(right[2],U[:,0][b],U[:,1][b]/U[:,0][b])


###############
#### SETUP ####
###############
# Grid
npt = 100
L = 1
dz = L/npt
zvec = np.linspace(dz/2,L-dz/2,npt)

#Fluid
GAM = 1.4

# Scheme
# Flux can be 'LF', 'LW', 'FORCE' ,'FLIC', ECUSP
order = 5
flux_type = 'LF'

# Data holders
# [rho,rho*u,E]
U = np.zeros([len(zvec),3])

# Case definition
caseNum = 7
left, right, cfl, tmax = defineCase(caseNum)
cfl = 0.25
f_0(U)
U0 = np.copy(U)

# Time
P0 = P_from_Ev(U[:,2]/U[:,0],U[:,0],U[:,1]/U[:,0])
a0 = np.sqrt(GAM * P0 / U[:,0])
v0 = U[:,1]/U[:,0]
lam0 = np.max(v0+a0)

dt = dz /lam0 * cfl
tc = 0

# Plots
f, axarr = plt.subplots(3, 2)
axarr[0,0].set_title("Density")
axarr[0,1].set_title("Pressure")
axarr[1,0].set_title("Temperature")
axarr[1,1].set_title("Velocity")
axarr[2,0].set_title("Mach number")

#######################
#### TIME MARCHING ####
#######################
idx = 0

if flux_type == 'ECUSP':
    compute_flux = lambda U,U0,dz,order: compute_ecusp_flux(U,U0,dz,order)
elif flux_type == 'LF':
    compute_flux = lambda U,U0,dz,order: compute_lf_flux(U,U0,dz,order)

while tc<tmax:
    UN = np.copy(U)  

    U1 = UN + dt * compute_flux(UN,U0,dz,order)  
    U1[0,:] = U0[0,:]
    U1[-1,:] = U0[-1,:]    
    
    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1,U0,dz,order) 
    U2[0,:] = U0[0,:]
    U2[-1,:] = U0[-1,:]   

    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2,U0,dz,order)
    UNP1[0,:] = U0[0,:]
    UNP1[-1,:] = U0[-1,:]   
       
    U = np.copy(UNP1)

    # Make sure we don't exceed the CFL condition by changing the time step
    rho = U[:,0]
    v = U[:,1] / U[:,0]
    E = U[:,2] / U[:,0]
    P = P_from_Ev(E,rho,v)
    asos = np.sqrt(GAM * P / rho)
    lam = np.max(np.abs(v)+asos)
        
    dt = cfl*dz/lam;
    
    if(tc + dt > tmax):
        dt = tmax - tc
    
    tc = tc+dt




axarr[0,0].plot(zvec,U[:,0],'o')
axarr[0,1].plot(zvec,P,'o')
axarr[1,0].plot(zvec,P/U[:,0],'o')
axarr[1,1].plot(zvec,U[:,1]/U[:,0],'o')
axarr[2,0].plot(zvec,U[:,1]/U[:,0]/asos,'o')

############### TRY WITH ECUSP
U = np.copy(U0)
tc = 0
flux_type = 'ECUSP'
dt = dz /lam0 * cfl
idx = 0

if flux_type == 'ECUSP':
    compute_flux = lambda U,U0,dz,order: compute_ecusp_flux(U,U0,dz,order)
elif flux_type == 'LF':
    compute_flux = lambda U,U0,dz,order: compute_lf_flux(U,U0,dz,order)

while tc<tmax:
    UN = np.copy(U)  
#    UNP1 = UN + dt * compute_flux(UN,U0,dz,order)
#    
    U1 = UN + dt * compute_flux(UN,U0,dz,order)   
    U1[0,:] = U0[0,:]
    U1[-1,:] = U0[-1,:]
    
    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1,U0,dz,order) 
    U2[0,:] = U0[0,:]
    U2[-1,:] = U0[-1,:]    
    
    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2,U0,dz,order)
    UNP1[0,:] = U0[0,:]
    UNP1[-1,:] = U0[-1,:] 
    
    U = np.copy(UNP1)

    # Make sure we don't exceed the CFL condition by changing the time step
    rho = U[:,0]
    v = U[:,1] / rho
    E = U[:,2] / rho
    P = P_from_Ev(E,rho,v)
    asos = np.sqrt(GAM * P / rho)
    lam = np.max(np.abs(v)+asos)
        
    dt = cfl*dz/lam;
    
    if(tc + dt > tmax):
        dt = tmax - tc
    
    tc = tc+dt

#axarr[0,0].plot(zvec,U[:,0],'x')
#axarr[0,1].plot(zvec,P,'x')
#axarr[1,0].plot(zvec,P/U[:,0],'x')
#axarr[1,1].plot(zvec,U[:,1]/U[:,0],'x')
#axarr[2,0].plot(zvec,U[:,1]/U[:,0]/asos,'x')

# Exact solution
caseStr = 'exact/' + str(caseNum) + '/case.csv'
data = np.genfromtxt(caseStr, delimiter =',',names = True)
#plt.plot(data['xe'],data['Me'],'k-')
#plt.plot(data['xe'],data['rhoe'],'k-')

axarr[0,0].plot(data['xe'],data['rhoe'],'k-')
axarr[0,1].plot(data['xe'],data['pe'],'k-')
axarr[1,0].plot(data['xe'],data['pe']/data['rhoe'],'k-')
axarr[1,1].plot(data['xe'],data['ue'],'k-')
axarr[2,0].plot(data['xe'],data['Me'],'k-')



# MATLAB solution
#qmatlab = np.genfromtxt('q.csv',delimiter=',')
#plt.plot(zvec,qmatlab[:,0],'^')




