""" File: eulerCharacteristics.py
Description: calculates the characteristics of the 2D Euler equation.
This includes the flux and the eigenvectors associated with it
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
from utils import P_from_Ev

GAM = 1.4

def compute_euler_flux(U,direction):
    rho = U[:,0]
    u = U[:,1] / rho
    v = U[:,2] / rho
    E = U[:,3] / rho
    
    P = P_from_Ev(E,rho,u,v)
   
    flx = np.zeros(U.shape)
    
    if direction == 'dx':
        flx[:,0] = rho*u
        flx[:,1] = rho*u**2 + P
        flx[:,2] = rho*u*v
        flx[:,3] = rho*E*u + P*u

    elif direction == 'dy':
        flx[:,0] = rho*v
        flx[:,1] = rho*u*v
        flx[:,2] = rho*v**2 + P
        flx[:,3] = rho*E*v + P*v
    
    return flx

def eigenvector_x(u,v,a,q,h,nunk):
    Rj = np.zeros((nunk,nunk))
    Lj = np.zeros((nunk,nunk))
    
    # Right eigenvector
    Rj[0,:] = np.ones((1,nunk))
    Rj[0,-1] = 0
            
    Rj[1,0] = u - a
    Rj[1,1] = u
    Rj[1,2] = u + a
    Rj[1,3] = 0
            
    Rj[2,0:3] = v
    Rj[2,3] = -1
     
    Rj[3,0] = h - a * u  
    Rj[3,1] = q
    Rj[3,2] = h + a * u
    Rj[3,3] = -v 

    # Left eigenvector
    Lj[0,0] = (GAM-1)*q + a*u
    Lj[0,1] = (1-GAM)*u - a
    Lj[0,2] = (1-GAM)*v
    Lj[0,3] = (GAM-1)
    Lj[0,:] /= (2*a**2)
    
    Lj[1,0] = a**2 - (GAM-1)*q
    Lj[1,1] = (GAM-1)*u
    Lj[1,2] = (GAM-1)*v
    Lj[1,3] = (1-GAM)
    Lj[1,:] /= a**2
    
    Lj[2,0] = (GAM-1)*q - a*u
    Lj[2,1] = (1-GAM)*u + a
    Lj[2,2] = (1-GAM)*v
    Lj[2,3] = (GAM-1)
    Lj[2,:] /= (2*a**2)
    
    Lj[3,0] = v
    Lj[3,2] = -1

    return Rj, Lj    

def eigenvector_y(u,v,a,q,h,nunk):
    Rj = np.zeros((nunk,nunk))
    Lj = np.zeros((nunk,nunk))
    
    # Right eigenvector
    Rj[0,:] = np.ones((1,nunk))
    Rj[0,-1] = 0
            
    Rj[1,0:3] = u 
    Rj[1,3] = 1
            
    Rj[2,0] = v - a
    Rj[2,1] = v
    Rj[2,2] = v + a
    Rj[2,3] = 0
     
    Rj[3,0] = h - a * v
    Rj[3,1] = q
    Rj[3,2] = h + a * v
    Rj[3,3] = u

    # Left eigenvector
    Lj[0,0] = (GAM-1)*q + a*v
    Lj[0,1] = (1-GAM)*u 
    Lj[0,2] = (1-GAM)*v - a
    Lj[0,3] = (GAM-1)
    Lj[0,:] /= (2*a**2)
    
    Lj[1,0] = a**2 - (GAM-1)*q
    Lj[1,1] = (GAM-1)*u
    Lj[1,2] = (GAM-1)*v
    Lj[1,3] = (1-GAM)
    Lj[1,:] /= a**2
    
    Lj[2,0] = (GAM-1)*q - a*v
    Lj[2,1] = (1-GAM)*u
    Lj[2,2] = (1-GAM)*v + a
    Lj[2,3] = (GAM-1)
    Lj[2,:] /= (2*a**2)
    
    Lj[3,0] = -u
    Lj[3,1] = 1

    return Rj, Lj    


def compute_eigenvector(U,U0,direction):
    rho = U[:,0]
    u = U[:,1] / rho
    v = U[:,2] / rho
    E = U[:,3] / rho
    
    P = P_from_Ev(E,rho,u,v)
    
    nunk = U.shape[1]
    nelem = U.shape[0]
    
    a = np.sqrt(GAM*P/rho)
    q = 1/2*(u**2 + v**2) # Dynamic pressure
    h = a**2/(GAM-1) + q # Enthalpy


    rho0 = U[:,0]
    u0 = U[:,1] / rho0
    v0 = U[:,2] / rho0
    E0 = U[:,3] / rho0
    
    P0 = P_from_Ev(E0,rho0,u0,v0)
    
    nunk = U0.shape[1]
    nelem = U0.shape[0]
    
    a0 = np.sqrt(GAM*P0/rho0)
    q0 = 1/2*(u0**2 + v0**2) # Dynamic pressure
    h0 = a0**2/(GAM-1) + q0 # Enthalpy

    
    Rjlist = []
    Ljlist = []
    
    if direction == 'dx':
        Rlhs0, Llhs0 = eigenvector_x(u0[0],v0[0],a0[0],q0[0],h0[0],nunk)
        for idx in range(nelem):
            Rj, Lj = eigenvector_x(u[idx],v[idx],a[idx],q[idx],h[idx],nunk)
            Rjlist.append(Rj)
            Ljlist.append(Lj)
        Rlhs0pre = None
        Llhs0pre = None
        
    
    elif direction == 'dy':
        # For the y-direction, the bottom boundary can either be pre or post-shock
        Rlhs0, Llhs0 = eigenvector_y(u0[0],v0[0],a0[0],q0[0],h0[0],nunk)
        Rlhs0pre, Llhs0pre = eigenvector_y(u0[-1],v0[-1],a0[-1],q0[-1],h0[-1],nunk)
        for idx in range(nelem):
            Rj, Lj = eigenvector_y(u[idx],v[idx],a[idx],q[idx],h[idx],nunk)
            Rjlist.append(Rj)
            Ljlist.append(Lj)
            
    Rj = Rjlist
    Lj = Ljlist

    return Rj,Lj,Rlhs0,Llhs0,Rlhs0pre,Llhs0pre
 

def to_characteristics(U,flx,U0,flx0,order,Lh,alpha,Nx,Ny,direction,options,tc,lambda_calc_char):
    nelem = U.shape[0]
    nunk = U.shape[1]    
    
    # Matrix holders
    V = np.zeros((nelem,order+1,nunk))
    VLF = np.zeros((nelem,order+1,nunk))
    H = np.zeros((nelem,order+1,nunk))
    
    # For all elements, evaluate R_{i+1/2}^-1 * [STENCIL]   
    # The conditional work for r = 2 (order 5)
    # We do the characteristics calculation for all elements for the whole stencil
    V,H,VLF = lambda_calc_char(U,flx,U0,flx0,order,Lh,alpha,direction,tc)
    


    return V,H,VLF