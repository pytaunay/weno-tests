""" File: computeLFFlux.py
Description: functions to calculate the Lax Friedrichs flux.

Author: Pierre-Yves Taunay
Date: November 2018
"""

import numpy as np
from scipy.linalg import block_diag

from utils import P_from_Ev
from wenoCharacteristicsCore import compute_lr

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
    
def to_characteristics(u,flx,U0,flx0,order,Rjinvp1h,alpha):
    nelem = u.shape[0]
    nunk = u.shape[1]    
    rweno = (int)((order-1)/2)    
    
    # Matrix holders
    v = np.zeros((nelem,order+1,nunk))
    vlf = np.zeros((nelem,order+1,nunk))
    g = np.zeros((nelem,order+1,nunk))
    
    # For all elements, evaluate R_{i+1/2}^-1 * [STENCIL]   
    # The conditional work for r = 1 and 2 (order 3 and 5). NOT HIGHER
    for idx in range(nelem):
        utmp = np.zeros((order+1,nunk))
        ftmp = np.zeros((order+1,nunk))
        
        # Grab the correct stencil data for the center domain
        if idx >= rweno and idx < nelem-rweno-1:
            utmp = u[idx-rweno:idx+rweno+2,:]
            ftmp = flx[idx-rweno:idx+rweno+2,:]
        
        # Beware of edges: treat the orders 3 and 5 separately
        # TODO: Parameterize!
        if order == 5:
            if idx == 0:
                utmp[0] = U0[0]
                utmp[1] = U0[0]
                ftmp[0] = flx0[0]
                ftmp[1] = flx0[0]
                utmp[2:6] = u[0:4]
                ftmp[2:6] = flx[0:4]
            elif idx == 1:
                utmp[0] = U0[0]
                ftmp[0] = flx0[0]
                utmp[1:6] = u[0:5]
                ftmp[1:6] = flx[0:5]  
            elif idx == nelem - 3:
                utmp[0:5] = u[nelem-5:nelem]
                ftmp[0:5] = flx[nelem-5:nelem]
                utmp[5] = U0[-1]
                ftmp[5] = flx0[-1]     
            elif idx == nelem - 2:
                utmp[0:4] = u[nelem-4:nelem]
                ftmp[0:4] = flx[nelem-4:nelem]
                utmp[4] = U0[-1]
                utmp[5] = U0[-1]
                
                ftmp[4] = flx0[-1]
                ftmp[5] = flx0[-1]     
            elif idx == nelem - 1:
                utmp[0:3] = u[nelem-3:nelem]
                ftmp[0:3] = flx[nelem-3:nelem]
                utmp[3] = U0[-1]
                utmp[4] = U0[-1]
                utmp[5] = U0[-1]
                ftmp[3] = flx0[-1]
                ftmp[4] = flx0[-1]
                ftmp[5] = flx0[-1]
        
        # TODO: Adjust order 3
        elif order == 3:
            if idx == 0:
                utmp[0] = U0[0]
                ftmp[0] = flx0[0]
                utmp[1:4] = u[0:3]
                ftmp[1:4] = flx[0:3]
            elif idx == nelem - 2:
                utmp[0:3] = u[nelem-3:nelem]
                ftmp[0:3] = flx[nelem-3:nelem]
                utmp[3] = U0[-1]
                ftmp[3] = flx0[-1]
                
            elif idx == nelem - 1:
                utmp[0:2] = u[nelem-2:nelem]
                ftmp[0:2] = flx[nelem-2:nelem]
                utmp[2] = U0[-1]
                utmp[3] = U0[-1]
                ftmp[2] = flx0[-1]
                ftmp[3] = flx0[-1]
        
        # Perform the matrix multiplication
        v[idx,:,:] = np.matmul(Rjinvp1h[idx],utmp.T).T
        g[idx,:,:] = np.matmul(Rjinvp1h[idx],ftmp.T).T
        vlf[idx,:,:] = np.matmul(np.diag(alpha),v[idx,:,:].T).T

    return v,g,vlf

def compute_eigenvector(U):
    rho = U[:,0]
    v = U[:,1] / rho
    E = U[:,2] / rho
    
    P = P_from_Ev(E,rho,v)
    
    a = np.sqrt(GAM*P/rho)
    Rjlist = []
    Rjinvlist = []
    
    for idx in np.arange(0,len(rho),1):
        Rj = np.zeros((3,3))
        Rjinv = np.zeros((3,3))
        # Right eigenvector
        Rj[0,:] = np.ones((1,3))
        
        Rj[1,0] = v[idx] - a[idx] 
        Rj[1,1] = v[idx]
        Rj[1,2] = v[idx] + a[idx]
        
        h = v[idx]**2/2 + a[idx]**2/(GAM-1)        
        Rj[2,0] = h - v[idx]*a[idx]  
        Rj[2,1] = v[idx]**2/2
        Rj[2,2] = h + v[idx]*a[idx]
        
        # Left eigenvector
        voa = v[idx]/a[idx] # v over a
        a2 = a[idx]**2
        Rjinv[0,0] = 1/2*(1/2*(GAM-1)*voa**2 + voa)
        Rjinv[0,1] = -1/(2*a[idx]) * ((GAM-1)*voa+1)
        Rjinv[0,2] = (GAM-1)/(2*a2)
        
        Rjinv[1,0] = 1-1/2*(GAM-1)*voa**2
        Rjinv[1,1] = (GAM-1)*v[idx]/a2
        Rjinv[1,2] = -(GAM-1)/a2
        
        Rjinv[2,0] = 1/2*(1/2*(GAM-1)*voa**2 - voa)
        Rjinv[2,1] = -1/(2*a[idx]) * ((GAM-1)*voa - 1)
        Rjinv[2,2] = (GAM-1)/(2*a2)
    

        Rjlist.append(Rj)
        Rjinvlist.append(Rjinv)
        
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
    
def compute_lfc_flux(u,U0,dz,order):
    ### Indices
    nelem = u.shape[0] # Number of elements
    nunk = u.shape[1] # Number of unknowns
    
    # What's r that corresponds to our WENO order?
    rweno = int((order-1)/2)
    
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

    ### Roe average at u^{i+-1/2}
#    up1h, um1h = roe_average(u,U0)
    up1h = (u+up1)*1/2
    

    ### LF Flux splitting: f = f^+ + f^- where
    ### f^{+-} = 1/2*(F[U] + alpha * U)
    ### We do reconstruction on the FLUX
    
    ### Calculate the LF split flux in all cells
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
    alpha = np.max(eig,axis=1) # GLOBAL MAX

    ### Characteristics
    # Eigenvectors: Compute the block matrix at each right eigenvector
    Rjp1h,Rjinvp1h = compute_eigenvector(up1h);
    
    ### Compute the flux
    flx = compute_euler_flux(u)
    flx0 = compute_euler_flux(U0)
    
    ######
    # I + 1/2
    ######   
    ### Transform into the characteristic domain
    v,g,vlf = to_characteristics(u,flx,U0,flx0,order,Rjinvp1h,alpha)
            
    ### Compute f+, f- at i+1/2 for all elements on the stencil
    FLXp = np.zeros((nelem,order,nunk))
    FLXm = np.zeros((nelem,order,nunk))
    
    FLXp = 1/2*(g[:,0:order,:] + vlf[:,0:order,:])
    FLXm = 1/2*(g[:,1:order+1,:] - vlf[:,1:order+1,:])
    
#    if order == 5:
#        FLXp = 1/2*(g[:,0:5,:] + vlf[:,0:5,:])
#        FLXm = 1/2*(g[:,1:6,:] - vlf[:,1:6,:])
#    elif order == 3:
#        FLXp = 1/2*(g[:,0:3,:] + vlf[:,0:3,:])
#        FLXm = 1/2*(g[:,1:6,:] - vlf[:,1:6,:])
#    FLXm = 1/2*(g - vlf)

    ### Reconstruct WENO in the characteristics domain    
    # Reconstruct the data on the stencil
    fp1hL, fm1hR = compute_lr(FLXp,order)
    fp3hL, fp1hR = compute_lr(FLXm,order)
    
    # Compute the flux
    FLXp1h = fp1hL + fp1hR
    FLXm1h = np.roll(FLXp1h,1,axis=0)

#    t1 = np.copy(FLXp1h)
#    t2 = np.copy(FLXm1h)

    ### Go back into the normal domain
    for idx in range(nelem):
        FLXp1h[idx] = np.matmul(Rjp1h[idx],FLXp1h[idx])
        if idx > 0:
            FLXm1h[idx] = np.matmul(Rjp1h[idx-1],FLXm1h[idx])
        else:
            FLXm1h[idx] = np.matmul(Rjp1h[0],FLXm1h[idx])
    
    return -1/dz*(FLXp1h - FLXm1h)
