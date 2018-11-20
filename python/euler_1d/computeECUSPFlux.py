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

def compute_ah(aL,aR):
    return 1/2*(aL+aR)

def compute_ecusp_flux(u,U0,dz,order):
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

    # Don't forget the BC for this problem
    up1[-1,:] = U0[-1,:]
    um1[0,:] = U0[0,:]
    
    
    
    ### Reconstruct the data on the stencil
    up1hL, um1hR = compute_lr(up1,u,um1,order)
                
    # Compute the RHS flux
    up1hR = np.roll(um1hR,-1,axis=0) # This will contain u_{i+1/2}^R
    um1hL = np.roll(up1hL,1,axis=0) # This will contain u_{i-1/2}^L
    
    # Don't forget the BC
    up1hR[-1,:] = U0[-1,:]
    um1hL[0,:] = U0[0,:] 
       
    ### i + 1/2
    # Calculate the convective flux
    vL = up1hL[:,1]/up1hL[:,0]
    vR = up1hR[:,1]/up1hR[:,0]
    
    PL = P_from_Ev(up1hL[:,-1]/up1hL[:,0],up1hL[:,0],vL)
    PR = P_from_Ev(up1hR[:,-1]/up1hR[:,0],up1hR[:,0],vR)
    
    ap1hL = np.sqrt(PL/up1hL[:,0])
    ap1hR = np.sqrt(PR/up1hR[:,0])
    
    ah = 1/2*(ap1hL+ap1hR)
    
    MaL = vL / ap1hL
    MaR = vR / ap1hR
    
    fp1hc = flux_convective(up1hL,up1hR,ah,MaL,MaR)
    
    # Calculate the pressure flux
    fp1hp = flux_pressure(up1hL,up1hR,ah,MaL,MaR)
    
    ### i - 1/2
    vL = um1hL[:,1]/um1hL[:,0]
    vR = um1hR[:,1]/um1hR[:,0]
    
    PL = P_from_Ev(um1hL[:,-1]/um1hL[:,0],um1hL[:,0],vL)
    PR = P_from_Ev(um1hR[:,-1]/um1hR[:,0],um1hR[:,0],vR)
    ap1hL = np.sqrt(PL/um1hL[:,0])
    ap1hR = np.sqrt(PR/um1hR[:,0])
    
    ah = 1/2*(ap1hL+ap1hR)
    
    MaL = vL / ap1hL
    MaR = vR / ap1hR
    
    fm1hc = flux_convective(um1hL,um1hR,ah,MaL,MaR)
    fm1hp = flux_pressure(um1hL,um1hR,ah,MaL,MaR)


    ### Sum the two
    fR = fp1hc + fp1hp
    fL = fm1hc + fm1hp
    # Add bool condition here  
    
    ### Make sure we use the applicable formula, which depends on the local speed

    # Supersonic traveling right
    
    # Supersonic traveling left
    


    return -1/dz * (fR-fL)
