""" File: computeLFFlux.py
Description: functions to calculate the Lax Friedrichs flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np

from utils import P_from_Ev
from wenoCharacteristicsCore import compute_lr
from computeLFCFlux import to_characteristics

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

def compute_eigenvector(U):
    nelem = U.shape[0]
    Rjlist = []
    Rjinvlist = []
    
    for idx in range(nelem):
        Rj = np.eye(3)
        Rjinv = np.eye(3)

        Rjlist.append(Rj)
        Rjinvlist.append(Rjinv)
        
    Rj = Rjlist
    Rjinv = Rjinvlist

    return Rj,Rjinv 

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
    alpha = np.max(eig,axis = 1) # GLOBAL MAX

    ### LF Flux splitting: f = f^+ + f^- where
    ### f^{+-} = 1/2*(F[U] + alpha * U)
    ### We do reconstruction on the FLUX
    flx = compute_euler_flux(u)
    vlf = np.matmul(np.diag(alpha),u.T).T   

    ### Characteristics
    # Eigenvectors: Compute the block matrix at each right eigenvector
    up1h = u
    Rjp1h,Rjinvp1h = compute_eigenvector(up1h);
    
    ### Compute the flux
    flx = compute_euler_flux(u)
    flx0 = compute_euler_flux(U0)
    
    nelem = u.shape[0]
    nunk = u.shape[1]
    
    ######
    # I + 1/2
    ######   
    ### Transform into the characteristic domain
    v,g,vlf = to_characteristics(u,flx,U0,flx0,order,Rjinvp1h,alpha)
            
    ### Compute f+, f- at i+1/2 for all elements on the stencil
    FLXp = np.zeros((nelem,order,nunk))
    FLXm = np.zeros((nelem,order,nunk))
    
    FLXp = 1/2*(g[:,0:order,:] + vlf[:,0:order,:])
    FLXm = 1/2*(g[:,1:order+1,:] - vlf[:,1:order+1,:])
    

    ### Reconstruct WENO in the characteristics domain    
    # Reconstruct the data on the stencil
    fp1hL, fm1hR = compute_lr(FLXp,order)
    fp3hL, fp1hR = compute_lr(FLXm,order)
    
    # Compute the flux
    FLXp1h = fp1hL + fp1hR
    FLXm1h = np.roll(FLXp1h,1,axis=0)

    ### Go back into the normal domain
    for idx in range(nelem):
        FLXp1h[idx] = np.matmul(Rjp1h[idx],FLXp1h[idx])
        if idx > 0:
            FLXm1h[idx] = np.matmul(Rjp1h[idx-1],FLXm1h[idx])
        else:
            FLXm1h[idx] = np.matmul(Rjp1h[0],FLXm1h[idx])
    
    return -1/dz*(FLXp1h - FLXm1h)
