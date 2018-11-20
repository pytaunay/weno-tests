""" File: computeECUSPConvectiveFlux.py
Description: calculates the convective flux that appears in the E-CUSP method
See https://arc.aiaa.org/doi/pdf/10.2514/1.3113
Zha and Hu, "Calculation of Transonic Internal Flows Using an Efficient
High-Resolution Upwind Scheme", AIAA Journal, 42, 2, 2004.

Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np

def compute_rhoh(rhoL,rhoR,uLp,uRp):
    return rhoL * uLp + rhoR * uRp

def compute_alphaL(UL,UR):
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    vL = UL[:,1]
    vR = UR[:,1]
    
    EL = UL[:,2]
    ER = UR[:,2]

    PL = P_from_Ev(EL/rhoL,rhoL,vL)
    PR = P_from_Ev(ER/rhoR,rhoR,vR)
    
    num = 2*(PL/rhoL)
    den = (PL/rhoL) + (PR/rhoR)
    
    return num/den

def compute_alphaR(UL,UR):
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    vL = UL[:,1]
    vR = UR[:,1]
    
    EL = UL[:,2]
    ER = UR[:,2]

    PL = P_from_Ev(EL/rhoL,rhoL,vL)
    PR = P_from_Ev(ER/rhoR,rhoR,vR)
    
    num = 2*(PR/rhoR)
    den = (PL/rhoL) + (PR/rhoR)
    
    return num/den

def compute_uLp(ah,MaL,UL,UR):
    aL = compute_alphaL(UL,UR)
    
    t1 = 1/2*(MaL + np.abs(MaL))
    t2 = aL * (1/4*(MaL+1)**2-1/2*(MaL+np.abs(MaL)))  
    
    return ah*(t1+t2)

def compute_uRp(ah,MaR,UL,UR):
    aR = compute_alphaR(UL,UR)
    
    t1 = 1/2*(MaR - np.abs(MaR))
    t2 = aR * (-1/4*(MaR-1)**2-1/2*(MaR-np.abs(MaR)))  
    
    return ah*(t1+t2)    

def flux_convective(UL,UR,ah,MaL,MaR):
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    # Compute uLplus, uRplus
    uLp = compute_uLp(ah,MaL,UL,UR)
    uRp = compute_uRp(ah,MaR,UL,UR)
    
    # Calculate rho_u1h
    rhoh = compute_rhoh(rhoL,rhoR,uLp,uRp)
    
    # Calculate qL, qR
    qL = UL
    qL[:,0] /= UL[:,0]
    qL[:,1] /= UL[:,1]
    qL[:,2] /= UL[:,2]
    
    qR = UR
    qR[:,0] /= UR[:,0]
    qR[:,1] /= UR[:,1]
    qR[:,2] /= UR[:,2]
    
    return 1/2 * (rhoh*(qR+qL) - np.abs(rhoh)*(qR-qL))