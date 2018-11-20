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

# Initial conditions
def f_0(U):
    b = (zvec < 0.5)
    U[:,0][b] = 1
    U[:,1][b] = 0
    U[:,2][b] = rhoE_from_Pv(1,U[:,0][b],U[:,1][b]/U[:,0][b])
    
    b = (zvec>=0.5)
    U[:,0][b] = 0.125
    U[:,1][b] = 0
    U[:,2][b] = rhoE_from_Pv(0.1,U[:,0][b],U[:,1][b]/U[:,0][b])
    
    
f_0(U)
U0 = U

# Time
P0 = P_from_Ev(U[:,2]/U[:,0],U[:,0],U[:,1]/U[:,0])
a0 = np.sqrt(GAM * P0 / U[:,0])
v0 = U[:,1]/U[:,0]
lam = np.max(v0+a0)

cfl = 0.9
dt = dz /lam * cfl
tmax = 0.1
#tmax = 0.1
tc = 0


#######################
#### TIME MARCHING ####
#######################
idx = 0

if flux_type == 'ECUSP':
    compute_flux = lambda U,U0,dz,order: compute_ecusp_flux(U,U0,dz,order)
elif flux_type == 'LF':
    compute_flux = lambda U,U0,dz,order: compute_lf_flux(U,U0,dz,order)

while tc<tmax:
    UN = U  

    U1 = UN + dt * compute_flux(UN,U0,dz,order)  
    U1[0,:] = U0[0,:]
    U1[-1,:] = U0[-1,:]    
    
    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1,U0,dz,order) 
    U2[0,:] = U0[0,:]
    U2[-1,:] = U0[-1,:]   

    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2,U0,dz,order)
    UNP1[0,:] = U0[0,:]
    UNP1[-1,:] = U0[-1,:]   
       
    U = UNP1

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

plt.plot(zvec,U[:,0],'o')
#plt.plot(zvec,U0[:,0],'k')

############### TRY WITH ECUSP
UN = U0
tc = 0
flux_type = 'ECUSP'
idx = 0

if flux_type == 'ECUSP':
    compute_flux = lambda U,U0,dz,order: compute_ecusp_flux(U,U0,dz,order)
elif flux_type == 'LF':
    compute_flux = lambda U,U0,dz,order: compute_lf_flux(U,U0,dz,order)

while tc<tmax:
    UN = U  
#    UNP1 = UN + dt * compute_flux(UN,U0,dz,order)
    U1 = UN + dt * compute_flux(UN,U0,dz,order)   
    U1[0,:] = U0[0,:]
    U1[-1,:] = U0[-1,:]
    
    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1,U0,dz,order) 
    U2[0,:] = U0[0,:]
    U2[-1,:] = U0[-1,:]    
    
    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2,U0,dz,order)
    UNP1[0,:] = U0[0,:]
    UNP1[-1,:] = U0[-1,:] 
    
    U = UNP1

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

    
    tc = tc+dt


plt.plot(zvec,U[:,0],'x')

# Exact solution
rhoe = np.genfromtxt('rhoe.csv',delimiter=',')
xe = np.genfromtxt('xe.csv',delimiter=',')
plt.plot(xe,rhoe,'k-')

# MATLAB solution
qmatlab = np.genfromtxt('qsod.csv',delimiter=',')
plt.plot(zvec,qmatlab[:,0],'^')

