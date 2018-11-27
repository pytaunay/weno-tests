""" File: initialize.py
Description: builds the grid for different cases and populates variables
Author: Pierre-Yves Taunay
Date: November 2018

"""

import numpy as np

def initialize(options):
    if options['case'] == 'double-mach':
        U,xvec,yvec,grid = initialize_dm(options)
   
    return U,xvec,yvec,grid     
        
def initialize_dm(options):
    Nx = options['grid'][0]
    Ny = options['grid'][1]
    
    nelem = Nx*Ny # Total number of points
    nunk = options['Upost'].shape[0] # Number of unknowns
    
    Lx = options['xlim'][1]-options['xlim'][0]
    Ly = options['ylim'][1]-options['ylim'][0]
    
    dx = Lx/Nx
    dy = Ly/Ny
    
    xvec = np.linspace(dx/2,Lx-dx/2,Nx)
    yvec = np.linspace(dy/2,Ly-dy/2,Ny)

    U = np.zeros((nelem,nunk))
    grid = np.zeros((nelem,2))
    
    for jj in range(Ny):
        yc = yvec[jj]
        for ii in range(Nx):
            xc = xvec[ii]
            
            slope = np.tan(options['angleshock'])
            ylim = slope*(xc - options['xshock'])
            
            if xc <= options['xshock']:
                U[ii+Nx*jj] = options['Upost']
            elif xc > options['xshock'] and yc >= ylim:
                U[ii+Nx*jj] = options['Upost']
            else:
                U[ii+Nx*jj] = options['Upre']
             
            grid[ii+Nx*jj] = np.array([xc,yc])    
                
    return U,xvec,yvec,grid