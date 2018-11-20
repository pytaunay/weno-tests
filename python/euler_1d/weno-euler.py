""" WENO scheme for the Euler 1D equation
The flux implemented is the ECUSP flux
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np
import matplotlib.pyplot as plt

###############
#### SETUP ####
###############
# Grid
npt = 100
L = 1
dz = L/npt
zvec = np.linspace(dz/2,L-dz/2,npt)
EPS = 1e-16

# Time
dt = dz / 1 * 0.4
tmax = 0.14
tc = 0

#Fluid
GAM = 5/3

# Scheme
# Flux can be 'LF', 'LW', 'FORCE' ,'FLIC', ECUSP
order = 3
flux_type = 'ECUSP'

# Data holders
# [rho,rho*u,E]
U = np.zeros([len(zvec),3])

# Useful functions
def P_from_Ev(E,rho,v):
    return (GAM-1)*(E-1/2*v**2)
    
def E_from_Pv(P,rho,v):
    return P/(GAM-1) + 1/2*v**2


# Initial conditions
def f_0(U):
    b = (zvec < 0.5)
    U[:,0][b] = 1
    U[:,1][b] = 0
    U[:,2][b] = E_from_Pv(1,U[:,0][b],U[:,1][b])
    
    b = (zvec>=0.5)
    U[:,0][b] = 0.125
    U[:,1][b] = 0
    U[:,2][b] = E_from_Pv(0.1,U[:,0][b],U[:,1][b])
    
    
f_0(U)
U0 = U

#######################
#### TIME MARCHING ####
#######################
idx = 0
### WENO 3
# vi+1/2[0]^L : 1/2 v_i     + 1/2 v_{i+1}
# vi+1/2[1]^L : -1/2 v_{i-1} + 3/2 v_i
# vi-1/2[0]^R : 3/2 v_i     - 1/2 v_{i+1}
# vi-1/2[1]^R : 1/2 v_{i-1} + 1/2 v_i
def compute_weights(up1,u,um1):
    if order == 3:
        d0 = 2/3
        d1 = 1/3
        
        beta0 = (up1-u)**2
        beta1 = (u-um1)**2
        
        alpha0 = d0 / (EPS+beta0)**2
        alpha1 = d1 / (EPS+beta1)**2
        
        alphat0 = d1 / (EPS+beta0)**2
        alphat1 = d0 / (EPS+beta1)**2
        
        alphasum = alpha0+alpha1
        alphatsum = alphat0 + alphat1
        
        w0 = alpha0 / alphasum
        w1 = alpha1 / alphasum
        
        wt0 = alphat0 / alphatsum
        wt1 = alphat1 / alphatsum
    
        return w0,w1,wt0,wt1
    elif order == 5:
        up2 = np.roll(u,-2)
        um2 = np.roll(u,2)
        
        d0 = 3/10
        d1 = 3/5
        d2 = 1/10
        
        beta0 = 13/12*(u-2*up1+up2)**2 + 1/4*(3*u-4*up1+up2)**2
        beta1 = 13/12*(um1-2*u+up1)**2 + 1/4*(um1-up1)**2
        beta2 = 13/12*(um2-2*um1+u)**2 + 1/4*(um2-4*um1+3*u)**2
        
        alpha0 = d0/(EPS+beta0)**2
        alpha1 = d1/(EPS+beta1)**2
        alpha2 = d2/(EPS+beta2)**2
        
        alphat0 = d2/(EPS+beta0)**2
        alphat1 = d1/(EPS+beta1)**2
        alphat2 = d0/(EPS+beta2)**2
        
        alphasum = alpha0 + alpha1 + alpha2
        alphatsum = alphat0 + alphat1 + alphat2
        
        w0 = alpha0/alphasum
        w1 = alpha1/alphasum
        w2 = alpha2/alphasum
        
        wt0 = alphat0/alphatsum
        wt1 = alphat1/alphatsum
        wt2 = alphat2/alphatsum
        
        return w0,w1,w2,wt0,wt1,wt2

def compute_lr(up1,u,um1):
    if order == 3:
        u0p = 1/2*u + 1/2*up1
        u1p = -1/2*um1 + 3/2*u
        
        u0m = 3/2*u - 1/2*up1
        u1m = 1/2*um1 + 1/2*u
        
        w0,w1,wt0,wt1 = compute_weights(up1,u,um1)
        
        uL = w0*u0p + w1*u1p
        uR = wt0*u0m + wt1*u1m
    elif order == 5:
        up2 = np.roll(up1,-1)
        um2 = np.roll(um1,1)
        
        u0m = 11/6*u - 7/6*up1 + 1/3*up2
        u1m = 1/3*um1 + 5/6*u - 1/6*up1
        u2m = -1/6*um2 + 5/6*um1 + 1/3*u
        
        u0p = 1/3*u + 5/6*up1 - 1/6*up2
        u1p = -1/6*um1 + 5/6*u + 1/3*up1
        u2p = 1/3*um2 -7/6*um1 + 11/6*u
        
        w0,w1,w2,wt0,wt1,wt2 = compute_weights(up1,u,um1)
        
        uL = w0*u0p + w1*u1p + w2*u2p
        uR = wt0*u0m + wt1*u1m + wt2*u2m
    
    return uL,uR

def compute_rhoh(rhoL,rhoR,uLp,uRp):
    return rhoL * uLp + rhoR * uRp

def compute_ah(aL,aR):
    return 1/2*(aL+aR)

def compute_flux(u):
    # u_{i+1}, u_{i-1}
    up1 = np.roll(u,-1,axis=0)
    um1 = np.roll(u,1,axis=0)  

    # Don't forget the BC for this problem
    up1[-1,:] = U0[-1,:]
    um1[0,:] = U0[0,:]
    
    ### Reconstruct the data on the stencil
    uL, uR = compute_lr(up1,u,um1)
            
    # Compute the RHS flux
    up1h = np.roll(uR,-1) # This will contain u_{i+1/2}^R
    um1h = np.roll(uL,1) # This will contain u_{i-1/2}^L
    
    # Don't forget the BC
    up1h[-1,:] = U0[-1,:]
    um1h[0,:] = U0[0,:] 
    
    ### Necessary intermediate steps
    rho = u[:,0]
    v = u[:,1]/u[:,0]
    E = u[:,2]/u[:,0]
    P = P_from_Ev(E,rho,v)
    
    # What's the local speed of sound and Mach number everywhere?
    a = np.sqrt(P/rho)
    aR = np.roll(a,-1)
    aL = np.roll(a,1)
    P0L = P_from_Ev(U0[0,-1]/U0[0,0],U0[0,0],U0[0,1]/U0[0,0])
    P0R = P_from_Ev(U0[-1,-1]/U0[-1,0],U0[-1,0],U0[-1,1]/U0[-1,0])
    
    aR[-1] = np.sqrt(/rho0) 
    aL[0] = np.sqrt(P0/rho0)
    
    ah = 1/2*(aR+aL)
    
    Ma = v/a
    MaR = np.roll(Ma,-1)
    MaL = np.roll(Ma,1)
    MaR[-1] = 0
    MaL[0] = 0
    
    
    
    
    fpR = 0
    fpL = 0
    
#    if flux_type == 'LF':
#        fpR = compute_flux_lf(uL,up1h)    
#        fpL = compute_flux_lf(um1h,uR)  
#    elif flux_type == 'ECUSP':
#        fpR = compute_flux_ecusp(uL,up1h)
#        fpL = compute_flux_ecusp(um1h,uR)

    return -1/dz * (fpR-fpL)

def compute_flux_ecusp():
    return 1

def compute_flux_lf(uL,uR):
    ### Left, right
    fL = flux(uL)
    fR = flux(uR)
    
    alpha = 1 # Derivative of flux
    
    return 1/2*(fL+fR-alpha*(uR-uL))    

def compute_flux_lw(uL,uR):
    alpha = 1
    u_lw = 1/2 * (uL+uR) - 1/2*alpha*(flux(uR)-flux(uL))
    
    return flux(u_lw)

def compute_flux_force(uL,uR):
    f_lf = compute_flux_lf(uL,uR)
    f_lw = compute_flux_lw(uL,uR)
    
    return 1/2*(f_lf + f_lw)
    
#
while tc<tmax:
    UN = U  
    U1 = UN + dt * compute_flux(UN)   
    U2 = 3/4*UN + 1/4*U1 + 1/4* dt * compute_flux(U1) 
    UNP1 = 1/3*UN + 2/3*U2 + 2/3 * dt * compute_flux(U2)
    
    U = UNP1
    
    tc = tc+dt
#    
#
#plt.plot(zvec,u0,'-')
#plt.plot(zvec,uvec,'o')
