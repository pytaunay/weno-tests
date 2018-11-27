""" WENO scheme for the Euler 2D equation
The flux implemented is the ECUSP flux
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
import matplotlib.pyplot as plt

from utils import P_from_Ev, rhoE_from_Pv,calculate_dt
from computeLFCFlux import compute_lfc_flux
from defineCase import defineCase
from initialize import initialize

###############
#### SETUP ####
###############
# Fluid
GAM = 1.4

# Scheme
order = 5

# Case definition
case = 'double-mach'
options = defineCase(case)

U0,xvec,yvec,grid,dx,dy = initialize(options)
xv, yv = np.meshgrid(xvec,yvec)

### Time
tmax = options['tmax']
cfl = options['cfl']
dt = calculate_dt(U0,dx,dy,cfl)
tc = 0

# Resample the data on a mesh grid
# dataZ = griddata(grid, U0[:,1], (xv, yv), method='nearest')

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
compute_flux = lambda U,U0,dz,order: compute_lfc_flux(U,U0,dz,order)

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
axarr[2,1].plot(zvec,U[:,2]/U[:,0],'o')





