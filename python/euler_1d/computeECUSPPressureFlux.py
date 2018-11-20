""" File: computeECUSPPressureFlux.py
Description: calculates the pressure flux that appears in the E-CUSP method
See https://arc.aiaa.org/doi/pdf/10.2514/1.3113
Zha and Hu, "Calculation of Transonic Internal Flows Using an Efficient
High-Resolution Upwind Scheme", AIAA Journal, 42, 2, 2004.

Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
from utils import P_from_Ev

GAM = 1.4
"""
Description: calculates P^+ and P^- (Eq. 22)
P^+ is calculated on the left of the interface
P^- is calculated on the right
"""
def compute_Ppm(MaL,MaR):
    alpha = 3/16
    Pp = 1/4 * (MaL+1)**2*(2-MaL) + alpha * MaL * (MaL**2-1)**2
    
    Pm = 1/4 * (MaR-1)**2*(2+MaR) - alpha * MaR * (MaR**2-1)**2

    return Pp,Pm


def flux_pressure(UL,UR,ah,MaL,MaR):
#    rhoL = UL[:,0]
#    rhoR = UR[:,0]
#    
#    vL = UL[:,1]/rhoL
#    vR = UR[:,1]/rhoR
#    
#    EL = UL[:,2]/rhoL
#    ER = UR[:,2]/rhoR
#    
#    PL = P_from_Ev(EL,rhoL,vL)
#    PR = P_from_Ev(ER,rhoR,vR)
#    aL = np.sqrt(GAM*PL/rhoL)
#    aR = np.sqrt(GAM*PR/rhoR)  
#    
#    MaLloc = vL/aL
#    MaRloc = vR/aR
#    
#    qL = np.zeros(UL.shape)
#    qR = np.zeros(UR.shape)
#    
#    qL[:,1] = PL * MaLloc
#    qL[:,2] = PL * aL
#    
#    qR[:,1] = PR * MaRloc
#    qR[:,2] = PR * aR
#    
#    return 1/2 * (qR-qL)
    
    # Get the pressure from primitive variables
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    vL = UL[:,1]/rhoL
    vR = UR[:,1]/rhoR
    
    EL = UL[:,2] / rhoL
    ER = UR[:,2] / rhoR

    PL = P_from_Ev(EL,rhoL,vL)
    PR = P_from_Ev(ER,rhoR,vR)
    
    Pp, Pm = compute_Ppm(MaL,MaR)
    
    Fp = np.zeros(UL.shape)
    Fp[:,1] = Pp * PL + Pm * PR
    Fp[:,2] = 1/2 *( PL*(vL + ah) + PR*(vR-ah) )
       
    return Fp