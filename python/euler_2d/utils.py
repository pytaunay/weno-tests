""" File: utils.py
Description: usefule functions
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np

GAM = 1.4

# Useful functions
def P_from_Ev(E,rho,u,v):
    return (GAM-1)*(rho*E-1/2*rho*(u**2+v**2))
    
def rhoE_from_Pv(P,rho,u,v):
    return P/(GAM-1) + 1/2*rho*(u**2+v**2)

def calculate_dt(U,dx,dy,cfl):
    ### Calculate the physical quantities
    rho = U[:,0]
    u = U[:,1]/rho
    v = U[:,2]/rho
    E = U[:,3]/rho
    
    P = P_from_Ev(E,rho,u,v)    
    a = np.sqrt(GAM*P/rho)
    
    ### Max eigenvalues
    # x-direction
    l1 = np.abs(u-a)   
    l2 = np.abs(u)
    l3 = np.abs(u+a)
    
    eig = np.vstack((l1,l2,l3))
    
    lx_max = np.max(eig)
    
    # y-direction
    l1 = np.abs(v-a)   
    l2 = np.abs(v)
    l3 = np.abs(v+a)
    
    eig = np.vstack((l1,l2,l3))
    
    ly_max = np.max(eig)
        
    ### CFL condition
    tx = 1/dx*lx_max
    ty = 1/dy*ly_max
#    
    dt = cfl/(tx+ty)
#    dt = 0.145*dx/np.max(a)

    return dt