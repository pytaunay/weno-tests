""" WENO scheme for the Euler 2D equation
The flux implemented is the ECUSP flux
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
import matplotlib.pyplot as plt

from utils import P_from_Ev, rhoE_from_Pv, calculate_dt
from computeLFCFlux import compute_lfc_flux
from dmStencilCharacteristics import dmStencilCharacteristics
from dmBoundaryCondition import dmBoundaryCondition
from defineCase import defineCase
from initialize import initialize

from scipy.interpolate import griddata

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
Nx = len(xvec)
Ny = len(yvec)

### Time
tmax = options['tmax']
cfl = options['cfl']
dt = calculate_dt(U0,dx,dy,cfl)
tc = 0
tmax = dt

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
lambda_calc_char = lambda U,flx,U0,flx0,order,Lh,alpha,direction,tc: dmStencilCharacteristics(U,flx,U0,flx0,order,Lh,alpha,Nx,Ny,direction,options,tc)
compute_flux = lambda U,tc: compute_lfc_flux(U,U0,dx,dy,Nx,Ny,order,options,tc,lambda_calc_char)
lambda_bc = lambda U,tc: dmBoundaryCondition(U,U0,options,tc)


U = np.copy(U0)
UNP1 = np.zeros(U0.shape)

nelem = U.shape[0]
while tc<tmax:
    print(tc,dt)
    UN = np.copy(U)  

    UNP1 = UN + dt * compute_flux(UN,tc)  
#    U1 = UN + dt * compute_flux(UN,tc)  
#    lambda_bc(U1,tc)
##    
#    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1,tc)
#    lambda_bc(U2,tc)
#    
#    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2,tc)
#    lambda_bc(UNP1,tc)
            
    lambda_bc(UNP1,tc)
    U = np.copy(UNP1)

    # Make sure we don't exceed the CFL condition by changing the time step
    rho = U[:,0]
    u = U[:,1] / U[:,0]
    v = U[:,2] / U[:,0]
    E = U[:,3] / U[:,0]
    P = P_from_Ev(E,rho,u,v)
    asos = np.sqrt(GAM * P / rho)
    lam = np.max(np.abs(v)+asos)
        
    
    
#    dt = calculate_dt(U,dx,dy,cfl)
#    print(tc,dt)
    if(tc + dt > tmax):
        dt = tmax - tc
    
    tc = tc+dt


rhoZ = griddata(grid, rho, (xv, yv), method='nearest')
uZ = griddata(grid, u, (xv, yv), method='nearest')
vZ = griddata(grid, v, (xv, yv), method='nearest')
EZ = griddata(grid, E, (xv, yv), method='nearest')
PZ = griddata(grid, P, (xv, yv), method='nearest')
TZ = PZ/rhoZ
aZ = np.sqrt(GAM*PZ/rhoZ)
MZ = np.sqrt(uZ**2+vZ**2)/aZ

axarr[0,0].contour(xv,yv,uZ)
axarr[0,1].contour(xv,yv,PZ)
axarr[1,0].contour(xv,yv,TZ)
axarr[1,1].contour(xv,yv,uZ)
axarr[2,0].contour(xv,yv,MZ)
axarr[2,1].contour(xv,yv,EZ)





