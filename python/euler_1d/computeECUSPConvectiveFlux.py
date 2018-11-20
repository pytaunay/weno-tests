""" File: computeECUSPConvectiveFlux.py
Description: calculates the convective flux that appears in the E-CUSP method
See https://arc.aiaa.org/doi/pdf/10.2514/1.3113
Zha and Hu, "Calculation of Transonic Internal Flows Using an Efficient
High-Resolution Upwind Scheme", AIAA Journal, 42, 2, 2004.

Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
from utils import P_from_Ev

#def compute_rhoh(rhoL,rhoR,uLp,uRm):
#    return rhoL * uLp + rhoR * uRm

def compute_alphaL(UL,UR):
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    vL = UL[:,1] / rhoL
    vR = UR[:,1] / rhoR
    
    EL = UL[:,2] / rhoL
    ER = UR[:,2] / rhoR

    PL = P_from_Ev(EL,rhoL,vL)
    PR = P_from_Ev(ER,rhoR,vR)
    
    num = 2*(PL/rhoL)
    den = (PL/rhoL) + (PR/rhoR)
    
    return num/den

def compute_alphaR(UL,UR):
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    vL = UL[:,1] / rhoL
    vR = UR[:,1] / rhoR
    
    EL = UL[:,2] / rhoL
    ER = UR[:,2] / rhoR

    PL = P_from_Ev(EL,rhoL,vL)
    PR = P_from_Ev(ER,rhoR,vR)
    
    num = 2*(PR/rhoR)
    den = (PL/rhoL) + (PR/rhoR)
    
    return num/den

def compute_uLp(ah,MaL,UL,UR):
    aL = compute_alphaL(UL,UR)
    
    t1 = 1/2*(MaL + np.abs(MaL))
    t2 = aL * (1/4*(MaL+1)**2-1/2*(MaL+np.abs(MaL)))  
    
    return ah*(t1+t2)

def compute_uRm(ah,MaR,UL,UR):
    aR = compute_alphaR(UL,UR)
    
    t1 = 1/2*(MaR - np.abs(MaR))
    t2 = aR * (-1/4*(MaR-1)**2-1/2*(MaR-np.abs(MaR)))  
    
    return ah*(t1+t2)    

def flux_convective(UL,UR,ah,MaL,MaR):
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    # Compute uLplus, uRplus
    uLp = compute_uLp(ah,MaL,UL,UR)
    uRm = compute_uRm(ah,MaR,UL,UR)
    
    # Calculate rho_u1h
    rhouh = rhoL * uLp + rhoR * uRm
    
    # Calculate qL, qR
    qL = UL
    qL[:,0] /= rhoL
    qL[:,1] /= rhoL
    qL[:,2] /= rhoL
    
    qR = UR
    qR[:,0] /= rhoR
    qR[:,1] /= rhoR
    qR[:,2] /= rhoR
    
    qsum = qR + qL
    qdif = qR - qL
    
    for idx in np.arange(0,3,1):
        qsum[:,idx] *= rhouh
        qdif[:,idx] *= np.abs(rhouh)
    
    return 1/2 * (qsum - qdif)