""" File: computeECUSPFlux.py
Description: functions to calculate the ECUSP flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np
from utils import P_from_Ev, E_from_Pv
from wenoCore import compute_lr

def compute_ah(aL,aR):
    return 1/2*(aL+aR)

def compute_flux(u):
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

    # Don't forget the BC for this problem
    up1[-1,:] = U0[-1,:]
    um1[0,:] = U0[0,:]
    
    ### Reconstruct the data on the stencil
    uL, uR = compute_lr(up1,u,um1)
            
    # Compute the RHS flux
    up1h = np.roll(uR,-1) # This will contain u_{i+1/2}^R
    um1h = np.roll(uL,1) # This will contain u_{i-1/2}^L
    
    # Don't forget the BC
    up1h[-1,:] = U0[-1,:]
    um1h[0,:] = U0[0,:] 
    
    ### Necessary intermediate steps
    rho = u[:,0]
    v = u[:,1]/u[:,0]
    E = u[:,2]/u[:,0]
    P = P_from_Ev(E,rho,v)
    
    # What's the local speed of sound and Mach number everywhere?
    a = np.sqrt(P/rho)
    aR = np.roll(a,-1)
    aL = np.roll(a,1)
    P0L = P_from_Ev(U0[0,-1]/U0[0,0],U0[0,0],U0[0,1]/U0[0,0])
    P0R = P_from_Ev(U0[-1,-1]/U0[-1,0],U0[-1,0],U0[-1,1]/U0[-1,0])
    
    aR[-1] = np.sqrt(/rho0) 
    aL[0] = np.sqrt(P0/rho0)
    
    ah = 1/2*(aR+aL)
    
    Ma = v/a
    MaR = np.roll(Ma,-1)
    MaL = np.roll(Ma,1)
    MaR[-1] = 0
    MaL[0] = 0
    
    
    
    
    fpR = 0
    fpL = 0
    
#    if flux_type == 'LF':
#        fpR = compute_flux_lf(uL,up1h)    
#        fpL = compute_flux_lf(um1h,uR)  
#    elif flux_type == 'ECUSP':
#        fpR = compute_flux_ecusp(uL,up1h)
#        fpL = compute_flux_ecusp(um1h,uR)

    return -1/dz * (fpR-fpL)

def compute_flux_ecusp():
    return 1