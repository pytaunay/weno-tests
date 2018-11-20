""" File: computeLFFlux.py
Description: functions to calculate the Lax Friedrichs flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np
from utils import P_from_Ev
from wenoCore import compute_lr

GAM = 5/3

def compute_euler_flux(U):
    rho = U[:,0]
    v = U[:,1] / rho
    E = U[:,2] / rho
    
    P = P_from_Ev(E,rho,v)
   
    flx = np.zeros(U.shape)
    
    flx[:,0] = rho*v
    flx[:,1] = rho*v**2 + P
    flx[:,2] = rho*E*v + P*v
    
    return flx
    

def compute_lf_flux(u,U0,dz,order):
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

    # Don't forget the BC for this problem
    up1[-1,:] = 1e20 * np.ones((1,3)) #U0[-1,:]
    um1[0,:] = 1e20 * np.ones((1,3)) #U0[0,:]
         
    ### Reconstruct the data on the stencil
    up1hL, um1hR = compute_lr(up1,u,um1,order)
                
    # Compute the RHS flux
    up1hR = np.roll(um1hR,-1,axis=0) # This will contain u_{i+1/2}^R
    um1hL = np.roll(up1hL,1,axis=0) # This will contain u_{i-1/2}^L
    
   
    # Don't forget the BC
    up1hR[-1,:] = U0[-1,:]
    um1hL[0,:] = U0[0,:] 
       
    ### i + 1/2   
    vL = up1hL[:,1]/up1hL[:,0]
    vR = up1hR[:,1]/up1hR[:,0]
    
    PL = P_from_Ev(up1hL[:,-1]/up1hL[:,0],up1hL[:,0],vL)
    PR = P_from_Ev(up1hR[:,-1]/up1hR[:,0],up1hR[:,0],vR)
    
    ap1hL = np.sqrt(GAM*PL/up1hL[:,0])
    ap1hR = np.sqrt(GAM*PR/up1hR[:,0])

    # We want the max eigenvalue of the Jacobian evaluated both at L and R     
    l1L = np.abs(vL-ap1hL)
    l1R = np.abs(vR-ap1hR)
    
    l2L = np.abs(vL)
    l2R = np.abs(vR)

    l3L = np.abs(vL+ap1hL)
    l3R = np.abs(vR+ap1hR)
    
    eig = np.vstack((l1L,l1R,l2L,l2R,l3L,l3R))
    eig = np.transpose(eig)
    alpha = np.max(eig,axis=1)
    
    fR = compute_euler_flux(up1hL) + compute_euler_flux(up1hR)
    for idx in np.arange(0,3,1):
        fR[:,idx] -= alpha * (up1hR[:,idx]-up1hL[:,idx])
    
    fR *= 1/2

    ### i - 1/2
    vL = um1hL[:,1]/um1hL[:,0]
    vR = um1hR[:,1]/um1hR[:,0]
    
    PL = P_from_Ev(um1hL[:,-1]/um1hL[:,0],um1hL[:,0],vL)
    PR = P_from_Ev(um1hR[:,-1]/um1hR[:,0],um1hR[:,0],vR)
    am1hL = np.sqrt(GAM*PL/um1hL[:,0])
    am1hR = np.sqrt(GAM*PR/um1hR[:,0])
 
    # We want the max eigenvalue of the Jacobian evaluated both at L and R     
    l1L = np.abs(vL-am1hL)
    l1R = np.abs(vR-am1hR)
    
    l2L = np.abs(vL)
    l2R = np.abs(vR)

    l3L = np.abs(vL+am1hL)
    l3R = np.abs(vR+am1hR)
    
    eig = np.vstack((l1L,l1R,l2L,l2R,l3L,l3R))
    eig = np.transpose(eig)
    alpha = np.max(eig,axis=1)
    
    fL = compute_euler_flux(um1hL) + compute_euler_flux(um1hR)
    for idx in np.arange(0,3,1):
        fL[:,idx] -= alpha * (um1hR[:,idx]-um1hL[:,idx])    
    
    fL *= 1/2
    return -1/dz * (fR-fL)
