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
    # What's r that corresponds to our WENO order?
    rweno = int((order-1)/2)
    stencil_size = order
    
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

    ### Roe average at u^{i+-1/2}
#    up1h, um1h = roe_average(u,U0)
    up1h = (u+up1)*1/2
    um1h = (u+um1)*1/2
    

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
    Rjm1h,Rjinvm1h = compute_eigenvector(um1h);

    # Transform into the characteristic domain
    nelem = u.shape[0] # Number of elements in domain
    nunk = u.shape[1] # Number of unknowns
    
    v = np.zeros((nelem,stencil_size,nunk))
    vlf = np.zeros((nelem,stencil_size,nunk))
    g = np.zeros((nelem,stencil_size,nunk))
    
    flx = compute_euler_flux(u)
    flx0 = compute_euler_flux(U0)
    
    ######
    # I + 1/2
    ######   
    # For all elements, evaluate R_{i+1/2}^-1 * [STENCIL]   
    # The conditional work for r = 1 and 2 (order 3 and 5). NOT HIGHER
    for idx in range(nelem):
        utmp = np.zeros((stencil_size,nunk))
        ftmp = np.zeros((stencil_size,nunk))
        
        # Grab the correct stencil data
        if idx >= rweno and idx < nelem-rweno:
            utmp = u[idx-rweno:idx+rweno+1,:]
            flxtmp = flx[idx-rweno:idx+rweno+1,:]
        
        # Beware of edges
        elif idx == rweno-1 
            if rweno-1>0: # Order 5
                utmp[0,:] = U0[0,:]
                ftmp[0,:] = flx0[0,:]
                utmp[rweno-1:2*rweno+1,:] = u[idx-rweno+1:idx+rweno+1,:]
                ftmp[rweno-1:2*rweno+1,:] = flx[idx-rweno+1:idx+rweno+1,:]  
            elif rweno-1 == 0: # Order 3
                utmp[0,:] = U0[0,:]
                ftmp[0,:] = flx0[0,:]
                utmp[rweno-1+1:idx+rweno+2,:] = u[idx-rweno+1:idx+rweno+2,:]
                ftmp[rweno-1+1:idx+rweno+2,:] = flx[idx-rweno+1:idx+rweno+2,:] 
                
        elif idx == 0 and rweno > 1: # Order 5: we treated the case 0 above
            utmp[0,:] = U0[0,:]
            ftmp[0,:] = flx0[0,:]
            utmp[1,:] = U0[0,:]
            ftmp[1,:] = flx0[0,:]            
            
            utmp[rweno:2*rweno+1,:] = u[idx:idx+rweno,:]
            ftmp[rweno:2*rweno+1,:] = flx[idx:idx+rweno,:]
            
        elif idx == nelem - rweno - 1 and nelem - rweno -1 not nelem - 1:
            
        elif idx == nelem - 1:
          
            
        v[idx,:,:] = np.matmul(Rjinvp1h[idx],utmp.T).T
        g[idx,:,:] = np.matmul(Rjinvp1h[idx],flxtmp.T).T
        vlf[idx,:,:] = np.matmul(np.diag(alpha),v[idx,:,:].T).T
            
            
    ### Reconstruct WENO in the characteristics domain
    # f^+, f^- at i+1/2
    FLXp = np.zeros(u.shape)
    FLXm = np.zeros(u.shape)
    for idx in np.arange(0,3,1):
        FLXp[:,idx] = 1/2*(g[:,idx] + vlf[:,idx])
        FLXm[:,idx] = 1/2*(g[:,idx] - vlf[:,idx])
    
    # i + 1/2 ^+
    fp1 = np.roll(FLXp,-1,axis=0)
    fm1 = np.roll(FLXp,1,axis=0)
    
    # Reconstruct the data on the stencil
    fp1hL, fm1hR = compute_lr(fp1,FLXp,fm1,order)
    fp1hR = np.roll(fm1hR,-1,axis=0)
                
    # Compute the RHS flux
    fiph1_P = fp1hL
    
    ### i+1/2 ^-
    fp1 = np.roll(FLXm,-1,axis=0)
    fm1 = np.roll(FLXm,1,axis=0)
    
    # Reconstruct the data on the stencil
    fp1hL, fm1hR = compute_lr(fp1,FLXm,fm1,order)
    fp1hR = np.roll(fm1hR,-1,axis=0)
    fiph1_M = fp1hR
    
    ### Go back to component domain
    for idx in np.arange(0,nelem,1):
        fiph1_P[idx,:] = np.matmul(Rjp1h[idx],fiph1_P[idx,:])
        fiph1_M[idx,:] = np.matmul(Rjp1h[idx],fiph1_M[idx,:])

    FLXp1h = fiph1_P + fiph1_M
#    FLXm1h = np.roll(FLXp1h,1,axis=0)
    
    ######
    # I - 1/2
    ######
    nelem = u.shape[0]
    for idx in np.arange(0,nelem,1):
        v[idx,:] = np.matmul(Rjinvm1h[idx],u[idx,:])
        g[idx,:] = np.matmul(Rjinvm1h[idx],flx[idx,:])
        vlf[idx,:] = np.matmul(Rjinvm1h[idx],u[idx,:])
        vlf[idx,:] = np.matmul(np.diag(alpha),vlf[idx,:])

    ### Reconstruct WENO in the characteristics domain
    # f^+, f^- at i+1/2
    FLXp = np.zeros(u.shape)
    FLXm = np.zeros(u.shape)
    for idx in np.arange(0,3,1):
        FLXp[:,idx] = 1/2*(g[:,idx] + vlf[:,idx])
        FLXm[:,idx] = 1/2*(g[:,idx] - vlf[:,idx])
    
    # i + 1/2 ^+
    fp1 = np.roll(FLXp,-1,axis=0)
    fm1 = np.roll(FLXp,1,axis=0)
    
    # Reconstruct the data on the stencil
    fp1hL, fm1hR = compute_lr(fp1,FLXp,fm1,order)
    fp1hR = np.roll(fm1hR,-1,axis=0)
                
    # Compute the RHS flux
    fiph1_M = fm1hR
    
    ### i+1/2 ^-
    fp1 = np.roll(FLXm,-1,axis=0)
    fm1 = np.roll(FLXm,1,axis=0)
    
    # Reconstruct the data on the stencil
    fp1hL, fm1hR = compute_lr(fp1,FLXm,fm1,order)
    fm1hL = np.roll(fp1hL,1,axis=0)
    fiph1_P = fm1hL
    
    ### Go back to component domain
    for idx in np.arange(0,nelem,1):
        fiph1_P[idx,:] = np.matmul(Rjm1h[idx],fiph1_P[idx,:])
        fiph1_M[idx,:] = np.matmul(Rjm1h[idx],fiph1_M[idx,:])

    FLXm1h = fiph1_P + fiph1_M

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

    return -1/dz*(FLXp1h - FLXm1h)
