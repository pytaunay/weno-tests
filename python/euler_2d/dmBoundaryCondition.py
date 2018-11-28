""" dmBoundaryCondition.py
Imposes the boundary conditions for the double Mach problem

"""
import numpy as np

def dmBoundaryCondition(U,U0,options,tc):
    nelem = U.shape[0]
    Nx = options['grid'][0]
    Ny = options['grid'][1]
    
    dx = options['xlim'][1] / options['grid'][0] 
    xs = options['xshock'] + (1+20*tc)/np.sqrt(3) # Shock location
    
    for idx in range(nelem):
        yidx = (int)(idx/Nx) # Increases by 1 for every row
        xc = dx/2 + (idx%Nx) * dx
        
        if idx % Nx == 0: # Left boundary
            U[idx] = options['Upost']
        elif yidx == 0:
            if xc < options['xshock']:
                U[idx] = options['Upost']
        elif yidx == Ny-1:
            if xc < xs:
                U[idx] = options['Upost']
            else:
                U[idx] = options['Upre']