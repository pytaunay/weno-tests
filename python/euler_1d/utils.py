""" File: utils.py
Description: usefule functions
Author: Pierre-Yves Taunay
Date: November 2018

"""

GAM = 5/3

# Useful functions
def P_from_Ev(E,rho,v):
    return (GAM-1)*(E-1/2*v**2)
    
def E_from_Pv(P,rho,v):
    return P/(GAM-1) + 1/2*v**2