""" File: initialize.py
Description: builds the grid for different cases and populates variables
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np

def initialize(options):
    if options['case'] == 'double-mach':
        U,xvec,yvec,grid,dx,dy = initialize_dm(options)
   
    return U,xvec,yvec,grid,dx,dy   
        
def initialize_dm(options):
    Nx = options['grid'][0]
    Ny = options['grid'][1]
    
    Nxtot = Nx + 2
    Nytot = Ny + 2
    
    # Total number of points WITH the boundaries included
    # The boundaries are at x = 0, y = 0, x = L, y = L
    nelem = Nxtot * Nytot
    nunk = options['Upost'].shape[0] # Number of unknowns
    
    Lx = options['xlim'][1]-options['xlim'][0]
    Ly = options['ylim'][1]-options['ylim'][0]
    
    dx = Lx/(Nxtot-1)
    dy = Ly/(Nytot-1)
    
    #xvec = np.linspace(dx/2,Lx-dx/2,Nx)
    #yvec = np.linspace(dy/2,Ly-dy/2,Ny)
    xvec = np.linspace(0,Lx,Nxtot)
    yvec = np.linspace(0,Ly,Nytot)

    U = np.zeros((nelem,nunk))
    grid = np.zeros((nelem,2))
    
    for jj in range(Nytot):
        yc = yvec[jj]
        for ii in range(Nxtot):
            xc = xvec[ii]
            
            slope = np.tan(options['angleshock'])
            ylim = slope*(xc - options['xshock'])
            
            if xc <= options['xshock']:
                U[ii+Nxtot*jj] = options['Upost']
            elif xc > options['xshock'] and yc >= ylim:
                U[ii+Nxtot*jj] = options['Upost']
            else:
                U[ii+Nxtot*jj] = options['Upre']
             
            grid[ii+Nxtot*jj] = np.array([xc,yc])    
                
    return U,xvec,yvec,grid,dx,dy
