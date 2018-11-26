""" File: computeLFFlux.py
Description: functions to calculate the Lax Friedrichs flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np
from scipy.linalg import block_diag

from utils import P_from_Ev
from wenoCore import compute_lr

GAM = 1.4

def compute_euler_flux(U):
    rho = U[:,0]
    v = U[:,1] / rho
    E = U[:,2] / rho
    
    P = P_from_Ev(E,rho,v)
   
    flx = np.zeros(U.shape)
    
    flx[:,0] = rho*v
    flx[:,1] = rho*v**2 + P
    flx[:,2] = rho*E*v + P*v
    
    return flx
    

def compute_right_eigenvector(U):
    rho = U[:,0]
    v = U[:,1] / rho
    E = U[:,2] / rho
    
    P = P_from_Ev(E,rho,v)
    
    a = np.sqrt(GAM*P/rho)
    Rjlist = []
    Rjinvlist = []
    
    for idx in np.arange(0,len(rho),1):
        Rj = np.zeros((3,3))
        Rj[0,:] = np.ones((1,3))
        
        Rj[1,0] = v[idx]
        Rj[1,1] = v[idx] - a[idx]
        Rj[1,2] = v[idx] + a[idx]
        
        h = v[idx]**2/2 + a[idx]**2/(GAM-1)        
        Rj[2,0] = v[idx]**2/2
        Rj[2,1] = h - v[idx]*a[idx]
        Rj[2,2] = h + v[idx]*a[idx]
        
        Rjlist.append(Rj)
        Rjinvlist.append(np.linalg.inv(Rj))
        
    Rj = Rjlist
    Rjinv = Rjinvlist

    return Rj,Rjinv     
 
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
    
def compute_lf_flux(u,U0,dz,order):
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

#    ### Roe average at u^{i+-1/2}
#    up1h, um1h = roe_average(u,U0)
#
#    ### LF Flux splitting: f = f^+ + f^- where
#    ### f^{+-} = 1/2*(F[U] + alpha * U)
#    ### We do reconstruction on the FLUX
#    
#    ### Calculate the LF split flux in all cells
#    rho = u[:,0]
#    v = u[:,1]/rho
#    E = u[:,2]/rho
#     
#    P = P_from_Ev(E,rho,v)
#    
#    a = np.sqrt(GAM*P/rho)
#    
#    # We want the max eigenvalue  
#    l1 = np.abs(v-a)   
#    l2 = np.abs(v)
#    l3 = np.abs(v+a)
#    
#    eig = np.vstack((l1,l2,l3))
#    eig = np.transpose(eig)
#    alpha = np.max(eig) # GLOBAL MAX
#
#    ### Characteristics
#    # Eigenvectors: Compute the block matrix at each right eigenvector
#    Rjp1h,Rjinvp1h = compute_right_eigenvector(up1h);
#    Rjm1h,Rjinvm1h = compute_right_eigenvector(um1h);
#
#    # Transform into the characteristic domain
#    v = np.copy(u)
#    g = np.copy(u) 
#    
#    flx = compute_euler_flux(u)
#    
#    ######
#    # I + 1/2
#    ######
#    nelem = u.shape[0]
#    for idx in np.arange(0,nelem,1):
#        v[idx,:] = np.matmul(Rjinvp1h[idx],u[idx,:])
#        g[idx,:] = np.matmul(Rjinvp1h[idx],flx[idx,:])
#
#    ### Reconstruct WENO in the characteristics domain
#    FLXp = np.zeros(u.shape)
#    FLXm = np.zeros(u.shape)
#    for idx in np.arange(0,3,1):
#        FLXp[:,idx] = 1/2*(g[:,idx] + alpha*v[:,idx])
#        FLXm[:,idx] = 1/2*(g[:,idx] - alpha*v[:,idx])
#    
#    ### i + 1/2 ^+
#    fp1 = np.roll(FLXp,-1,axis=0)
#    fm1 = np.roll(FLXp,1,axis=0)
#    
#    # Reconstruct the data on the stencil
#    fp1hL, fm1hR = compute_lr(fp1,FLXp,fm1,order)
#    fp1hR = np.roll(fm1hR,-1,axis=0)
#    #up1hL, um1hR = compute_lr(up1,u,um1,order)  
#    #up1hR = np.roll(um1hR,-1,axis=0)
#                
#    # Compute the RHS flux
#    fiph1_P = fp1hL
#    
#    ### i+1/2 ^-
#    fp1 = np.roll(FLXm,-1,axis=0)
#    fm1 = np.roll(FLXm,1,axis=0)
#    
#    # Reconstruct the data on the stencil
#    fp1hL, fm1hR = compute_lr(fp1,FLXm,fm1,order)
#    fp1hR = np.roll(fm1hR,-1,axis=0)
#    fiph1_M = fp1hR
#    
#    ### Go back to component domain
#    for idx in np.arange(0,nelem,1):
#        fiph1_P[idx,:] = np.matmul(Rjp1h[idx],fiph1_P[idx,:])
#        fiph1_M[idx,:] = np.matmul(Rjp1h[idx],fiph1_M[idx,:])
#
#    FLXp1h = fiph1_P + fiph1_M
#    
#    ######
#    # I - 1/2
#    ######
#    nelem = u.shape[0]
#    for idx in np.arange(0,nelem,1):
#        v[idx,:] = np.matmul(Rjinvm1h[idx],u[idx,:])
#        g[idx,:] = np.matmul(Rjinvm1h[idx],flx[idx,:])
#
#    ### Reconstruct WENO in the characteristics domain
#    FLXp = np.zeros(u.shape)
#    FLXm = np.zeros(u.shape)
#    for idx in np.arange(0,3,1):
#        FLXp[:,idx] = 1/2*(g[:,idx] + alpha*v[:,idx])
#        FLXm[:,idx] = 1/2*(g[:,idx] - alpha*v[:,idx])
#    
#    ### i - 1/2 ^+
#    fp1 = np.roll(FLXp,-1,axis=0)
#    fm1 = np.roll(FLXp,1,axis=0)
#    
#    # Reconstruct the data on the stencil
#    fp1hL, fm1hR = compute_lr(fp1,FLXp,fm1,order)
#    fp1hR = np.roll(fm1hR,-1,axis=0)
#    #up1hL, um1hR = compute_lr(up1,u,um1,order)  
#    #up1hR = np.roll(um1hR,-1,axis=0)
#                
#    # Compute the RHS flux
#    fiph1_M = fm1hR
#    
#    ### i-1/2 ^-
#    fp1 = np.roll(FLXm,-1,axis=0)
#    fm1 = np.roll(FLXm,1,axis=0)
#    
#    # Reconstruct the data on the stencil
#    fp1hL, fm1hR = compute_lr(fp1,FLXm,fm1,order)
##    fp1hR = np.roll(fm1hR,-1,axis=0)
#    fm1hL = np.roll(fp1hL,1,axis=0)
#    fiph1_P = fm1hL
#    
#    ### Go back to component domain
#    for idx in np.arange(0,nelem,1):
#        fiph1_P[idx,:] = np.matmul(Rjm1h[idx],fiph1_P[idx,:])
#        fiph1_M[idx,:] = np.matmul(Rjm1h[idx],fiph1_M[idx,:])    
#        
#    FLXm1h = fiph1_P + fiph1_M
#
#    return -1/dz*(FLXp1h - FLXm1h)

    ### Calculate the max eigenvalue of all cells
    rho = u[:,0]
    v = u[:,1]/rho
    E = u[:,2]/rho
     
    P = P_from_Ev(E,rho,v)
    
    a = np.sqrt(GAM*P/rho)
    
    # We want the max eigenvalue  
    l1 = np.abs(v-a)   
    l2 = np.abs(v)
    l3 = np.abs(v+a)
    
    eig = np.vstack((l1,l2,l3))
    eig = np.transpose(eig)
    alpha = np.max(eig) # GLOBAL MAX

    ### LF Flux splitting: f = f^+ + f^- where
    ### f^{+-} = 1/2*(F[U] + alpha * U)
    ### We do reconstruction on the FLUX
    ### f^+
    # First calculate the flux
    flx = 1/2* (compute_euler_flux(u) + alpha*u)
    
    # f_{i+1}, f_{i-1}
    fp1 = np.roll(flx,-1,axis=0)
    fm1 = np.roll(flx,1,axis=0)     
    
    # WENO Reconstruction
    fp1hL, fm1hR = compute_lr(fp1,flx,fm1,order)
    fm1hL = np.roll(fp1hL,1,axis=0)
    
    # Flux
    fp1h = fp1hL
    fm1h = fm1hL
    
    ### f^-
    # First calculate the flux
    flx = 1/2* (compute_euler_flux(u) - alpha*u)
    
    # f_{i+1}, f_{i-1}
    fp1 = np.roll(flx,-1,axis=0)
    fm1 = np.roll(flx,1,axis=0)     
    
    # WENO Reconstruction
    fp1hL, fm1hR = compute_lr(fp1,flx,fm1,order)
    fp1hR = np.roll(fm1hR,-1,axis=0)   

    fp1h += fp1hR
    fm1h += fm1hR
    
    fR = fp1h
    fL = fm1h
    
    return -1/dz * (fR-fL)
