""" dmBoundaryCondition.py
Imposes the boundary conditions for the double Mach problem

"""
import numpy as np

def dmBoundaryCondition(U,U0,options,tc):
    Nx = options['grid'][0] + 2
    
    dx = options['xlim'][1] / options['grid'][0] 
    xs = options['xshock'] + (1+20*tc)/np.sqrt(3) # Shock location
    
    U[::Nx] = options['Upost'] # Left boundary
    U[Nx-1::Nx] = options['Upre'] # Right boundary
    
    
    for idx in range(1,Nx-1):
        xc = (idx % Nx) * dx
        # Bottom boundary
        if xc <= options['xshock']: # Past the shock (upstream)
            U[idx] = options['Upost']
        else: # After the shock (downstream): reflective wall
            U[idx,0] = U[idx+Nx,0]    
            U[idx,1] = U[idx+Nx,1] 
            U[idx,3] = U[idx+Nx,3] 

        # Top boundary               
        if xc < xs:
            U[idx] = options['Upost']
        else:
            U[idx] = options['Upre']
                
