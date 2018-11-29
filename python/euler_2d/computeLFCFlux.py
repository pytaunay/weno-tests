""" File: computeLFFlux.py
Description: functions to calculate the Lax Friedrichs flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np

from utils import P_from_Ev
from wenoCharacteristicsCore import compute_lr
from eulerCharacteristics import compute_eigenvector, compute_euler_flux, to_characteristics

GAM = 1.4

def roe_average(u,U0):
    ### i + 1/2
    rhoL = u[:,0]
    vL = u[:,1] / rhoL
    EL = u[:,2] / rhoL
    
    rhoR = np.roll(rhoL,-1)
    vR = np.roll(vL,-1)
    ER = np.roll(EL,-1)
    
    rhoR[-1] = U0[-1,0]
    vR[-1] = U0[-1,1] / U0[-1,0]
    ER[-1] = U0[-1,2] / U0[-1,0]

    sL = np.sqrt(rhoL)
    sR = np.sqrt(rhoR)

    rhoave = np.sqrt(rhoL*rhoR)
    vave = (sL*vL + sR*vR)/(sL+sR)
    Eave = (sL*EL + sR*ER)/(sL+sR)

    up1h = np.copy(u)
    up1h[:,0] = rhoave
    up1h[:,1] = rhoave * vave
    up1h[:,2] = rhoave * Eave

    ### i - 1/2
    rhoR = u[:,0]
    vR = u[:,1] / rhoR
    ER = u[:,2] / rhoR
    
    rhoL = np.roll(rhoR,-1)
    vL = np.roll(vR,-1)
    EL = np.roll(ER,-1)

    rhoL[0] = U0[0,0]
    vL[0] = U0[0,1]/U0[0,0]
    EL[0] = U0[0,2]/U0[0,0]


    sL = np.sqrt(rhoL)
    sR = np.sqrt(rhoR)

    rhoave = np.sqrt(rhoL*rhoR)
    vave = (sL*vL + sR*vR)/(sL+sR)
    Eave = (sL*EL + sR*ER)/(sL+sR)

    um1h = np.copy(u)
    um1h[:,0] = rhoave
    um1h[:,1] = rhoave * vave
    um1h[:,2] = rhoave * Eave

    return up1h,um1h   
    
### Calculate the flux in one-direction
def compute_lfc_flux_1D(U,U0,dx,dy,Nx,Ny,order,direction,options,tc,lambda_calc_char):
    nelem = U.shape[0]
    nunk = U.shape[1]
    
    ### Compute the average at the border      
    Up1h = compute_average(U,U0,dx,dy,Nx,Ny,direction)

    ### Calculate the LF split flux in all cells
    rho = U[:,0]
    u = U[:,1]/rho
    v = U[:,2]/rho
    E = U[:,3]/rho
     
    P = P_from_Ev(E,rho,u,v)
    a = np.sqrt(GAM*P/rho)
    
    # We want the max eigenvalue 
    if direction == 'dx':
        l1 = np.abs(u-a)   
        l2 = np.abs(u)
        l3 = np.abs(u+a)
        l4 = np.abs(u)
        
    elif direction == 'dy':
        l1 = np.abs(v-a)   
        l2 = np.abs(v)
        l3 = np.abs(v+a)
        l4 = np.abs(v)
    
    eig = np.vstack((l1,l2,l3,l4))
    alpha = np.max(eig,axis=1) # GLOBAL MAX

    ### Characteristics
    # Eigenvectors: Compute the block matrix at each right eigenvector
    if direction == 'dx':
        Rh,Lh,Rlhs0,Llhs0,_,_ = compute_eigenvector(Up1h,U0,direction)
    elif direction == 'dy':
        Rh,Lh,Rlhs0,Llhs0,Rlhs0pre,Llhs0pre = compute_eigenvector(Up1h,U0,direction)
    
    ### Compute the flux
    flx = compute_euler_flux(U,direction)
    flx0 = compute_euler_flux(U0,direction)

    ### Transform into the characteristic domain
    # V = R^-1 U, H = R^-1 * F or R^-1 * G, VLF = alpha * V
#    V,H,VLF = to_characteristics(U,flx,U0,flx0,order,Lh,alpha,Nx,Ny,direction)
    V,H,VLF = to_characteristics(U,flx,U0,flx0,order,Lh,alpha,Nx,Ny,direction,options,tc,lambda_calc_char)      
    
    ### Compute f+, f- at i+1/2 for all elements on the stencil
    FLXp = np.zeros((nelem,order,nunk))
    FLXm = np.zeros((nelem,order,nunk))
    
    FLXp = 1/2*(H[:,0:order,:] + VLF[:,0:order,:])
    FLXm = 1/2*(H[:,1:order+1,:] - VLF[:,1:order+1,:])

    fp1hL, fm1hR = compute_lr(FLXp,order)
    fp3hL, fp1hR = compute_lr(FLXm,order)
    
    # Compute the flux
    FLXp1h = fp1hL + fp1hR
    
    if direction == 'dx':
        FLXm1h = np.roll(FLXp1h,1,axis=0)
    elif direction == 'dy':
        FLXm1h = np.roll(FLXp1h,Nx,axis=0)
    
    ### Go back into the normal domain    
    if direction == 'dx':
        for idx in range(nelem):
            FLXp1h[idx] = np.matmul(Rh[idx],FLXp1h[idx])
            if idx % Nx != 0:
                FLXm1h[idx] = np.matmul(Rh[idx-1],FLXm1h[idx])
            else:
                FLXm1h[idx] = np.matmul(Rlhs0,FLXm1h[idx])    
    elif direction == 'dy':
        for idx in range(nelem):            
            FLXp1h[idx] = np.matmul(Rh[idx],FLXp1h[idx])
            if idx > Nx-1: # Any cell with an index higher than Nx is NOT on the bottom boundary
                FLXm1h[idx] = np.matmul(Rh[idx-Nx],FLXm1h[idx])
            else:
                FLXm1h[idx] = np.matmul(Rlhs0,FLXm1h[idx]) # We don't care about what's at the bottom boundary: it will be overwritten
    
    if direction =='dx':
        dz = dx
    elif direction =='dy':
        dz = dy
    
    return -1/dz*(FLXp1h - FLXm1h)
    
### Calculate the average at i (or j) + 1/2
def compute_average(U,U0,dx,dy,Nx,Ny,direction):
    Up1 = np.copy(U)
    if direction == 'dx':
        # Roll the x-direction
        Up1 = np.roll(U,-1,axis=0)

        # The right-hand edge now has the LHS values... 
        # This fixes it by setting the i+1 value to the i value for the cells
        # on the right boundary
        Up1[Nx-1::Nx] = U[Nx-1::Nx]
        
    elif direction == 'dy':
        # Roll the y-direction
        Up1 = np.roll(U,-Nx,axis=0)
        
        # But be careful about the top edge
        for ii in range(Nx):
            Up1[ii + (Ny-1)*Nx] = U[ii + (Ny-1)*Nx]
        
    Up1h = (U+Up1) / 2
    
    return Up1h
    

def compute_lfc_flux(U,U0,dx,dy,Nx,Ny,order,options,tc,lambda_calc_char):    
    ### Flux in the x and y direction
    FLXX = compute_lfc_flux_1D(U,U0,dx,dy,Nx,Ny,order,'dx',options,tc,lambda_calc_char)
    
    FLXY = compute_lfc_flux_1D(U,U0,dx,dy,Nx,Ny,order,'dy',options,tc,lambda_calc_char)

    return FLXX+FLXY