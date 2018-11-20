""" File: wenoCore.py
Description: computes the WENO reconstruction with either order 3 or 5

Author: Pierre-Yves Taunay
Date: November 2018

"""
import numpy as np

EPS = 1e-16


### WENO 3
# vi+1/2[0]^L : 1/2 v_i     + 1/2 v_{i+1}
# vi+1/2[1]^L : -1/2 v_{i-1} + 3/2 v_i
# vi-1/2[0]^R : 3/2 v_i     - 1/2 v_{i+1}
# vi-1/2[1]^R : 1/2 v_{i-1} + 1/2 v_i
def compute_weights(up1,u,um1,order):
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
        up2 = np.roll(up1,-1)
        um2 = np.roll(um1,1)
        
        up2[-1,:] = 1e20 * np.ones((1,3))
        um2[0,:] = 1e20 * np.ones((1,3))
        
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

def compute_lr(up1,u,um1,order):
    if order == 3:
        u0p = 1/2*u + 1/2*up1
        u1p = -1/2*um1 + 3/2*u
        
        u0m = 3/2*u - 1/2*up1
        u1m = 1/2*um1 + 1/2*u
        
        w0,w1,wt0,wt1 = compute_weights(up1,u,um1,order)
        
        uL = w0*u0p + w1*u1p
        uR = wt0*u0m + wt1*u1m
    elif order == 5:
        up2 = np.roll(up1,-1)
        um2 = np.roll(um1,1)
        
        up2[-1,:] = 1e20 * np.ones((1,3))
        um2[0,:] = 1e20 * np.ones((1,3))
        
        u0m = 11/6*u - 7/6*up1 + 1/3*up2
        u1m = 1/3*um1 + 5/6*u - 1/6*up1
        u2m = -1/6*um2 + 5/6*um1 + 1/3*u
        
        u0p = 1/3*u + 5/6*up1 - 1/6*up2
        u1p = -1/6*um1 + 5/6*u + 1/3*up1
        u2p = 1/3*um2 -7/6*um1 + 11/6*u
        
        w0,w1,w2,wt0,wt1,wt2 = compute_weights(up1,u,um1,order)
        
        uL = w0*u0p + w1*u1p + w2*u2p
        uR = wt0*u0m + wt1*u1m + wt2*u2m
    
    return uL,uR
