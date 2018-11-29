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
    ### Grid
    # x direction
    Nx = options['grid'][0]    
    Nxtot = Nx + 2
    Lx = options['xlim'][1]-options['xlim'][0]
    dx = Lx/(Nxtot-1)
    xvec = np.linspace(0,Lx,Nxtot)
    
    # y direction
    # override the Ny parameter to have dx = dy
    Ly = options['ylim'][1]-options['ylim'][0]
    dy = dx    
    yvec = np.arange(0,Ly+dx,dx)
    Nytot = len(yvec)
    
    # Total number of points WITH the boundaries included
    # The boundaries are at x = 0, y = 0, x = L, y = L
    nelem = Nxtot * Nytot
    nunk = options['Upost'].shape[0] # Number of unknowns


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
