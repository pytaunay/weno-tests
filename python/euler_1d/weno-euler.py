""" WENO scheme for the Euler 1D equation
The flux implemented is the ECUSP flux
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
import matplotlib.pyplot as plt

from utils import P_from_Ev, E_from_Pv

###############
#### SETUP ####
###############
# Grid
npt = 100
L = 1
dz = L/npt
zvec = np.linspace(dz/2,L-dz/2,npt)

# Time
dt = dz / 1 * 0.4
tmax = 0.14
tc = 0

#Fluid
GAM = 5/3

# Scheme
# Flux can be 'LF', 'LW', 'FORCE' ,'FLIC', ECUSP
order = 3
flux_type = 'ECUSP'

# Data holders
# [rho,rho*u,E]
U = np.zeros([len(zvec),3])

# Initial conditions
def f_0(U):
    b = (zvec < 0.5)
    U[:,0][b] = 1
    U[:,1][b] = 0
    U[:,2][b] = E_from_Pv(1,U[:,0][b],U[:,1][b])
    
    b = (zvec>=0.5)
    U[:,0][b] = 0.125
    U[:,1][b] = 0
    U[:,2][b] = E_from_Pv(0.1,U[:,0][b],U[:,1][b])
    
    
f_0(U)
U0 = U

#######################
#### TIME MARCHING ####
#######################
idx = 0

def compute_flux_lf(uL,uR):
    ### Left, right
    fL = flux(uL)
    fR = flux(uR)
    
    alpha = 1 # Derivative of flux
    
    return 1/2*(fL+fR-alpha*(uR-uL))    

def compute_flux_lw(uL,uR):
    alpha = 1
    u_lw = 1/2 * (uL+uR) - 1/2*alpha*(flux(uR)-flux(uL))
    
    return flux(u_lw)

def compute_flux_force(uL,uR):
    f_lf = compute_flux_lf(uL,uR)
    f_lw = compute_flux_lw(uL,uR)
    
    return 1/2*(f_lf + f_lw)
    
#
while tc<tmax:
    UN = U  
    U1 = UN + dt * compute_flux(UN)   
    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1) 
    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2)
    
    U = UNP1
    
    tc = tc+dt
#    
#
#plt.plot(zvec,u0,'-')
#plt.plot(zvec,uvec,'o')
