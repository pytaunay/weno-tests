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
    Lj[0,1] = (1-GAM)*u + a
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
    Lj[2,:] /= 2*a*2
    
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
    Lj[2,:] /= 2*a*2
    
    Lj[3,0] = -u
    Lj[3,2] = 1

    return Rj, Lj    


def compute_eigenvector(U,direction):
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
    
    Rjlist = []
    Ljlist = []
    
    if direction == 'dx':
        for idx in range(nelem):
            Rj, Lj = eigenvector_x(u[idx],v[idx],a[idx],q[idx],h[idx],nunk)
            Rjlist.append(Rj)
            Ljlist.append(Lj)
    
    elif direction == 'dy':
        for idx in range(nelem):
            Rj, Lj = eigenvector_y(u[idx],v[idx],a[idx],q[idx],h[idx],nunk)
            Rjlist.append(Rj)
            Ljlist.append(Lj)
            
    Rj = Rjlist
    Lj = Ljlist

    return Rj,Lj     
 

def to_characteristics(U,flx,U0,flx0,order,Lh,alpha,Nx,Ny,direction):
    nelem = U.shape[0]
    nunk = U.shape[1]    
    rweno = (int)((order-1)/2)    
    
    # Matrix holders
    V = np.zeros((nelem,order+1,nunk))
    VLF = np.zeros((nelem,order+1,nunk))
    H = np.zeros((nelem,order+1,nunk))
    
    # For all elements, evaluate R_{i+1/2}^-1 * [STENCIL]   
    # The conditional work for r = 1 and 2 (order 3 and 5). NOT HIGHER
    
    if direction == 'dx':
        for idx in range(nelem):
            Utmp = np.zeros((order+1,nunk))
            Ftmp = np.zeros((order+1,nunk))
                
            if order == 5:
                if idx % Nx == 0: # Left boundary, 2 ghost cells --- inlet
                    Utmp[0] = U0[idx]
                    Utmp[1] = U0[idx]
                    Ftmp[0] = flx0[idx]
                    Ftmp[1] = flx0[idx]
                    
                    Utmp[2:6] = U[idx:idx+rweno+2]
                    Ftmp[2:6] = flx[idx:idx+rweno+2]
                    
                elif idx % Nx == 1: # Left boundary, 1 ghost cell --- inlet
                    Utmp[0] = U0[idx]
                    Ftmp[0] = flx0[idx]
                    
                    Utmp[1:6] = U[idx-1:idx+rweno+2]
                    Ftmp[1:6] = flx[idx-1:idx+rweno+2]
                    
                elif idx % Nx == Nx-3: # Right boundary, 1 ghost cell --- outlet Neumann
                    Utmp[-1] = U[idx]
                    Ftmp[-1] = np.zeros((1,nunk))
                    
                    Utmp[0:5] = U[idx-rweno:idx+rweno+1]
                    Ftmp[0:5] = flx[idx-rweno:idx+rweno+1]
                    
                elif idx % Nx == Nx-2: # Right boundary, 2 ghost cells --- outlet Neumann
                    Utmp[-1] = U[idx]
                    Utmp[-2] = U[idx]
                    Ftmp[-1] = np.zeros((1,nunk))
                    Ftmp[-2] = np.zeros((1,nunk))
                    
                    Utmp[0:4] = U[idx-rweno:idx+rweno]
                    Ftmp[0:4] = flx[idx-rweno:idx+rweno]
                    
                elif idx % Nx == Nx - 1: # Right boundary, 3 ghost cells --- outlet Neumann
                    Utmp[-1] = U[idx]
                    Utmp[-2] = U[idx]
                    Utmp[-3] = U[idx]
                    Ftmp[-1] = np.zeros((1,nunk))
                    Ftmp[-2] = np.zeros((1,nunk))
                    Ftmp[-3] = np.zeros((1,nunk))
                    
                    Utmp[0:3] = U[idx-rweno:idx+rweno-1]
                    Ftmp[0:3] = flx[idx-rweno:idx+rweno-1]
                else: # Otherwise we are center of domain
                    Utmp[0:6] = U[idx-rweno:idx+rweno+2]
                    Ftmp[0:6] = flx[idx-rweno:idx+rweno+2]
                    
            # END IF ORDER == 5
            # TODO: ORDER 3
            # Perform the matrix multiplication
            V[idx,:,:] = np.matmul(Lh[idx],Utmp.T).T
            H[idx,:,:] = np.matmul(Lh[idx],Ftmp.T).T
            VLF[idx,:,:] = np.matmul(np.diag(alpha),V[idx,:,:].T).T  
        # END FOR NELEM LOOP
    
    elif direction == 'dy':
        for idx in range(nelem):
            Utmp = np.zeros((order+1,nunk))
            Ftmp = np.zeros((order+1,nunk))
            
            yidx = (int)(idx/Nx) # Increases by 1 for every row
            
            if order == 5:
                if yidx == 0: # Bottom boundary, 2 ghost cells --- reflective
                    Utmp[0] = U[idx]
                    Utmp[1] = U[idx]
                    Ftmp[0] = np.zeros((1,nunk))
                    Ftmp[1] = np.zeros((1,nunk))
                    
                    Utmp[2:6] = U[idx::Nx][0:4]
                    Ftmp[2:6] = flx[idx::Nx][0:4]
                    
                elif yidx == 1: # Bottom boundary, 1 ghost cell --- reflective
                    Utmp[0] = U[idx]
                    Ftmp[0] = np.zeros((1,nunk))
                    
                    Utmp[1:6] = U[idx::Nx][0:5]
                    Ftmp[1:6] = flx[idx::Nx][0:5]
                    
                elif yidx == Ny-3: # Top boundary, 1 ghost cell --- reflective
                    Utmp[-1] = U[idx + 2*Nx]
                    Ftmp[-1] = np.zeros((1,nunk))
                    
                    # Grab the bottom neighbors with idx::-Nx; limit to the last 5 with 0::5
                    # Since ::-Nx inverts the list, put it back in correct order with ::-1
                    # Example:
                    # top row index = 100, Nx = 10, current idx = 80
                    # idx + 2 *Nx::-10 = 100::-10 -> [100,90,80,70,60,50,40,30,20,10,0]
                    # 0:5 -> [100,90,80,70,60]
                    # ::-1 -> [60,70,80,90,100]
                    Utmp[0:5] = U[idx+2*Nx::-Nx][0:5][::-1]
                    Ftmp[0:5] = flx[idx+2*Nx::-Nx][0:5][::-1]
                    
                elif yidx == Ny-2: # Top boundary, 2 ghost cells --- reflective
                    Utmp[-1] = U[idx+Nx]
                    Utmp[-2] = U[idx+Nx]
                    Ftmp[-1] = np.zeros((1,nunk))
                    Ftmp[-2] = np.zeros((1,nunk))
                    
                    Utmp[0:4] = U[idx+Nx::-Nx][0:4][::-1]
                    Ftmp[0:4] = flx[idx+Nx::-Nx][0:4][::-1]
                    
                elif yidx == Ny - 1: # Top boundary, 3 ghost cells --- reflective
                    Utmp[-1] = U[idx]
                    Utmp[-2] = U[idx]
                    Utmp[-3] = U[idx]
                    Ftmp[-1] = np.zeros((1,nunk))
                    Ftmp[-2] = np.zeros((1,nunk))
                    Ftmp[-3] = np.zeros((1,nunk))
                    
                    Utmp[0:3] = U[idx::-Nx][0:3][::-1]
                    Ftmp[0:3] = flx[idx::-Nx][0:3][::-1]
                else:
                    Utmp[0:3] = U[idx::-Nx][0:3][::-1]
                    Utmp[3:6] = U[idx+Nx::Nx][0:3]
                    Ftmp[0:3] = flx[idx::-Nx][0:3][::-1]
                    Ftmp[3:6] = flx[idx+Nx::Nx][0:3]
                    
                
            # END IF ORDER == 5
            # TODO: ORDER 3
            
            # Perform the matrix multiplication
            V[idx,:,:] = np.matmul(Lh[idx],Utmp.T).T
            H[idx,:,:] = np.matmul(Lh[idx],Ftmp.T).T
            VLF[idx,:,:] = np.matmul(np.diag(alpha),V[idx,:,:].T).T    
        # END FOR NELEM LOOP
    # END DIRECTION CONDITIONAL


    return V,H,VLF