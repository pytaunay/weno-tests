""" dmBoundaryCondition.py
Imposes the boundary conditions for the double Mach problem

"""
import numpy as np

def dmBoundaryCondition(U,U0,options,tc):
    Nx = options['grid'][0] + 2
    Ny = options['grid'][1] + 2
    
    dx = options['xlim'][1] / Nx
    xs0 = options['xshock'] + options['ylim'][1] / np.sqrt(3) # Shock location at t=0
    vs = 20 / np.sqrt(3) 
    xs = xs0 + vs * tc
    
    U[::Nx] = options['Upost'] # Left boundary
    U[Nx-1::Nx] = options['Upre'] # Right boundary
    
    for idx in range(1,Nx-1):
        xc = idx * dx
        
        # Bottom boundary
        if xc < options['xshock']: # Past the shock (upstream)
            U[idx] = options['Upost']
        else: # After the shock (downstream): reflective wall
            U[idx,0] = U[idx+Nx,0]    
            U[idx,1] = U[idx+Nx,1] 
            U[idx,2] = 0
            U[idx,3] = U[idx+Nx,3] 

        # Top boundary               
        if xc < xs:
            U[Nx*(Ny-1)+idx] = options['Upost']
        else:
            U[Nx*(Ny-1)+idx] = options['Upre']
                
