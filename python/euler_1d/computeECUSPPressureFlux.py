""" File: computeECUSPPressureFlux.py
Description: calculates the pressure flux that appears in the E-CUSP method
See https://arc.aiaa.org/doi/pdf/10.2514/1.3113
Zha and Hu, "Calculation of Transonic Internal Flows Using an Efficient
High-Resolution Upwind Scheme", AIAA Journal, 42, 2, 2004.

Author: Pierre-Yves Taunay
Date: November 2018

"""


def flux_pressure(UL,UR):
    # Get the pressure from primitive variables
    rhoL = UL[:,0]
    rhoR = UR[:,0]
    
    vL = UL[:,1]
    vR = UR[:,1]
    
    EL = UL[:,2]
    ER = UR[:,2]

    PL = P_from_Ev(EL/rhoL,rhoL,vL)
    PR = P_from_Ev(ER/rhoR,rhoR,vR)
    
    return 1