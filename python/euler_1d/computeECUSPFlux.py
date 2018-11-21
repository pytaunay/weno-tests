""" File: computeECUSPFlux.py
Description: functions to calculate the ECUSP flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np
from utils import P_from_Ev
from wenoCore import compute_lr
from computeECUSPConvectiveFlux import flux_convective
from computeECUSPPressureFlux import flux_pressure

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
    

def compute_ah(aL,aR):
    return 1/2*(aL+aR)

def compute_ecusp_flux(u,U0,dz,order):
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  
    
    ### Reconstruct the data on the stencil
    up1hL, um1hR = compute_lr(up1,u,um1,order)  
    up1hR = np.roll(um1hR,-1,axis=0)
                
    # Compute the RHS flux
    um1hL = np.roll(up1hL,1,axis=0) # This will contain u_{i-1/2}^L
       
    ### i + 1/2
    # Calculate the convective flux
    rhoL = up1hL[:,0]
    rhoR = up1hR[:,0]

    
    vL = up1hL[:,1]/rhoL
    vR = up1hR[:,1]/rhoR
    
    EL = up1hL[:,2]/rhoL
    ER = up1hR[:,2]/rhoR
    
    PL = P_from_Ev(EL,rhoL,vL)
    PR = P_from_Ev(ER,rhoR,vR)
    
    ap1hL = np.sqrt(GAM*PL/rhoL)
    ap1hR = np.sqrt(GAM*PR/rhoR)
    
    ah = 1/2*(ap1hL+ap1hR)
 
    MaL = vL/ah
    MaR = vR/ah

    
    fp1hc = flux_convective(up1hL,up1hR,ah,MaL,MaR)    
    fp1hp = flux_pressure(up1hL,up1hR,ah,MaL,MaR)
    
    ### i - 1/2
    rhoL = um1hL[:,0]
    rhoR = um1hR[:,0]
    
    vL = um1hL[:,1]/rhoL
    vR = um1hR[:,1]/rhoR
    
    EL = um1hL[:,2]/rhoL
    ER = um1hR[:,2]/rhoR
    
    PL = P_from_Ev(EL,rhoL,vL)
    PR = P_from_Ev(ER,rhoR,vR)
    am1hL = np.sqrt(GAM*PL/rhoL)
    am1hR = np.sqrt(GAM*PR/rhoR)
    
    ah = 1/2*(am1hL+am1hR)

    MaL = vL/ah
    MaR = vR/ah
    
    fm1hc = flux_convective(um1hL,um1hR,ah,MaL,MaR)
    fm1hp = flux_pressure(um1hL,um1hR,ah,MaL,MaR)


    ### Calculate all types of fluxes
    # Subsonic
    fRsub = (fp1hc + fp1hp)
    fLsub = (fm1hc + fm1hp)
    
    
    # Supersonic traveling right
    fRsupr = compute_euler_flux(up1hL)
    fLsupr = compute_euler_flux(um1hL)
    
    # Supersonic travelight left
    fRsupl = compute_euler_flux(up1hR)
    fLsupl = compute_euler_flux(um1hR)    
    
    
    ### Make sure we use the applicable formula, which depends on the local speed
    # What's the local speed of sound?
    vloc = u[:,1]/u[:,0]
    Eloc = u[:,2]/u[:,0]
    rholoc = u[:,0]
    
    Ploc = P_from_Ev(Eloc,rholoc,vloc)
    aloc = np.sqrt(GAM*Ploc/rholoc)

    b1 = (np.abs(vloc) <= aloc)
    b2 = vloc > aloc
    b3 = vloc < -aloc
    
    fR = np.zeros(u.shape)
    fL = np.zeros(u.shape)
        
    for idx in np.arange(0,3,1):
        fR[:,idx] = fRsub[:,idx] * b1 + fRsupr[:,idx] * b2 + fRsupl[:,idx] * b3
        fL[:,idx] = fLsub[:,idx] * b1 + fLsupr[:,idx] * b2 + fLsupl[:,idx] * b3

#
#    print(b1,b2,b3)
#    print("------------")

#    print(fR-fL)

    return -1/dz * (fR-fL)
#    return -1/dz * (fRsub-fLsub)
