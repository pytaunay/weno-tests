""" File: computeLFFlux.py
Description: functions to calculate the Lax Friedrichs flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np
from scipy.linalg import block_diag

from utils import P_from_Ev
from wenoCore import compute_lr

GAM = 1.4

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
    ### Calculate the max eigenvalue of all cells
    rho = u[:,0]
    v = u[:,1]/rho
    E = u[:,2]/rho
     
    P = P_from_Ev(E,rho,v)
    
    a = np.sqrt(GAM*P/rho)
    
    # We want the max eigenvalue  
    l1 = np.abs(v-a)   
    l2 = np.abs(v)
    l3 = np.abs(v+a)
    
    eig = np.vstack((l1,l2,l3))
    eig = np.transpose(eig)
    alpha = np.max(eig) # GLOBAL MAX

    ### LF Flux splitting: f = f^+ + f^- where
    ### f^{+-} = 1/2*(F[U] + alpha * U)
    ### We do reconstruction on the FLUX
    ### f^+
    # First calculate the flux
    flx = 1/2* (compute_euler_flux(u) + alpha*u)
    
    # f_{i+1}, f_{i-1}
    fp1 = np.roll(flx,-1,axis=0)
    fm1 = np.roll(flx,1,axis=0)     
    
    # WENO Reconstruction
    fp1hL, fm1hR = compute_lr(fp1,flx,fm1,order)
    fm1hL = np.roll(fp1hL,1,axis=0)
    
    # Flux
    fp1h = fp1hL
    fm1h = fm1hL
    
    ### f^-
    # First calculate the flux
    flx = 1/2* (compute_euler_flux(u) - alpha*u)
    
    # f_{i+1}, f_{i-1}
    fp1 = np.roll(flx,-1,axis=0)
    fm1 = np.roll(flx,1,axis=0)     
    
    # WENO Reconstruction
    fp1hL, fm1hR = compute_lr(fp1,flx,fm1,order)
    fp1hR = np.roll(fm1hR,-1,axis=0)   

    fp1h += fp1hR
    fm1h += fm1hR
    
    fR = fp1h
    fL = fm1h
    
    return -1/dz * (fR-fL)
