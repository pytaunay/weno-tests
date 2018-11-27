""" File: utils.py
Description: usefule functions
Author: Pierre-Yves Taunay
Date: November 2018

"""

GAM = 1.4

# Useful functions
def P_from_Ev(E,rho,u,v):
    return (GAM-1)*(rho*E-1/2*rho*(u**2+v**2))
    
def rhoE_from_Pv(P,rho,u,v):
    return P/(GAM-1) + 1/2*rho*(u**2+v**2)