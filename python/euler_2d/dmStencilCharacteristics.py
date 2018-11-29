""" dnmtencilCharacteristics.py
Calculates the stencil characteristics decomposition for all elements

"""

import numpy as np

from eulerCharacteristics import compute_euler_flux

def dmStencilCharacteristics(U,flx,U0,flx0,order,Lh,alpha,Nx,Ny,direction,options,tc):
    nelem = U.shape[0]
    nunk = U.shape[1]    
    rweno = (int)((order-1)/2)    
    
    # Matrix holders
    V = np.zeros((nelem,order+1,nunk))
    VLF = np.zeros((nelem,order+1,nunk))
    H = np.zeros((nelem,order+1,nunk))    
    
    if direction == 'dx':
        for idx in range(nelem):
            Utmp = np.zeros((order+1,nunk))
            Ftmp = np.zeros((order+1,nunk))
                
            if order == 5:
                if idx % Nx == 0: # Left boundary, 2 ghost cells --- inlet
                    Utmp[0] = U0[idx]
                    Utmp[1] = U0[idx]
                    Ftmp[0] = flx0[idx]
                    Ftmp[1] = flx0[idx]
                    
                    Utmp[2:6] = U[idx:idx+rweno+2]
                    Ftmp[2:6] = flx[idx:idx+rweno+2]
                    
                elif idx % Nx == 1: # Left boundary, 1 ghost cell --- inlet
                    Utmp[0] = U0[idx]
                    Ftmp[0] = flx0[idx]
                    
                    Utmp[1:6] = U[idx-1:idx+rweno+2]
                    Ftmp[1:6] = flx[idx-1:idx+rweno+2]
                    
                elif idx % Nx == Nx-3: # Right boundary, 1 ghost cell --- outlet Neumann
                    Utmp[-1] = U[idx+2]
#                    Ftmp[-1] = np.zeros((1,nunk))
                    Ftmp[-1] = flx[idx+2]
                    
                    Utmp[0:5] = U[idx-rweno:idx+rweno+1]
                    Ftmp[0:5] = flx[idx-rweno:idx+rweno+1]
                    
                elif idx % Nx == Nx-2: # Right boundary, 2 ghost cells --- outlet Neumann
                    Utmp[-1] = U[idx+1]
                    Utmp[-2] = U[idx+1]
#                    Ftmp[-1] = np.zeros((1,nunk))
#                    Ftmp[-2] = np.zeros((1,nunk))
                    Ftmp[-1] = flx[idx+1] 
                    Ftmp[-2] = flx[idx+1]
                    
                    Utmp[0:4] = U[idx-rweno:idx+rweno]
                    Ftmp[0:4] = flx[idx-rweno:idx+rweno]
                    
                elif idx % Nx == Nx - 1: # Right boundary, 3 ghost cells --- outlet Neumann
                    Utmp[-1] = U[idx]
                    Utmp[-2] = U[idx]
                    Utmp[-3] = U[idx]
#                    Ftmp[-1] = np.zeros((1,nunk))
#                    Ftmp[-2] = np.zeros((1,nunk))
#                    Ftmp[-3] = np.zeros((1,nunk))
                    Ftmp[-1] = flx[idx]
                    Ftmp[-2] = flx[idx]
                    Ftmp[-3] = flx[idx]
                    
                    Utmp[0:3] = U[idx-rweno:idx+rweno-1]
                    Ftmp[0:3] = flx[idx-rweno:idx+rweno-1]
                else: # Otherwise we are center of domain
                    Utmp[0:6] = U[idx-rweno:idx+rweno+2]
                    Ftmp[0:6] = flx[idx-rweno:idx+rweno+2]
                    
            # END IF ORDER == 5
            # TODO: ORDER 3
            # Perform the matrix multiplication
            V[idx,:,:] = np.matmul(Lh[idx],Utmp.T).T
            H[idx,:,:] = np.matmul(Lh[idx],Ftmp.T).T
            VLF[idx,:,:] = np.matmul(np.diag(alpha),V[idx,:,:].T).T  
        # END FOR NELEM LOOP
    
    elif direction == 'dy':
        for idx in range(nelem):
            Utmp = np.zeros((order+1,nunk))
            Ftmp = np.zeros((order+1,nunk))
            
            yidx = (int)(idx/Nx) # Increases by 1 for every row
            dx = options['xlim'][1] / (Nx-1)
#            dy = options['ylim'][1] / options['grid'][1]
            xc = (idx%Nx) * dx # Location of cell along x axis
#            yc = dy/2 + yidx*dy
            
            # Current shock location
            xs = options['xshock'] + (1+20*tc)/np.sqrt(3)

            if idx == 524-Nx:
                print("532")
            
            if order == 5:
                if yidx == 0: # Bottom boundary, 2 ghost cells
                    if xc < options['xshock']: # If we are before the shock, enforce the upstream B.C.
                        Utmp[0] = U0[0]
                        Utmp[1] = U0[0]
                        Ftmp[0] = flx0[0]
                        Ftmp[1] = flx0[0]
                    else:
                        tmp = U[idx]
                        tmp[2] = 0 # Zero y-velocity
                        flxtmp = compute_euler_flux(np.array([tmp]),direction)
                        flxtmp = flxtmp[0]
                        
                        Utmp[0] = tmp
                        Utmp[1] = tmp
                        
                        Ftmp[0] = flxtmp
                        Ftmp[1] = flxtmp
                    
                    Utmp[2:6] = U[idx::Nx][0:4]
                    Ftmp[2:6] = flx[idx::Nx][0:4]
                    
                elif yidx == 1: # Bottom boundary, 1 ghost cell
                    if xc < options['xshock']: # If we are before the shock, enforce the upstream B.C.
                        Utmp[0] = U0[0]
                        Ftmp[0] = flx0[0]
                    else:
                        tmp = U[idx-Nx]
                        tmp[2] = 0 # Zero y-velocity
                        flxtmp = compute_euler_flux(np.array([tmp]),direction)
                        flxtmp = flxtmp[0]
                        
                        Utmp[0] = tmp
                        Ftmp[0] = flxtmp

                    Utmp[1:6] = U[idx-Nx::Nx][0:5]
                    Ftmp[1:6] = flx[idx-Nx::Nx][0:5]
                    
                elif yidx == Ny-3: # Top boundary, 1 ghost cell
                    if xc < xs:
                        Utmp[-1] = U0[0]
                        Ftmp[-1] = flx0[0]
                    else:
                        Utmp[-1] = U0[-1] # Top right value. LIMITATION: CAN'T RUN PAST THE END OF DOMAIN
                        Ftmp[-1] = flx0[-1]
                    
                    # Grab the bottom neighbors with idx::-Nx; limit to the last 5 with 0::5
                    # Since ::-Nx inverts the list, put it back in correct order with ::-1
                    # Example:
                    # top row index = 100, Nx = 10, current idx = 80
                    # idx + 2 *Nx::-10 = 100::-10 -> [100,90,80,70,60,50,40,30,20,10,0]
                    # 0:5 -> [100,90,80,70,60]
                    # ::-1 -> [60,70,80,90,100]
                    Utmp[0:5] = U[idx+2*Nx::-Nx][0:5][::-1]
                    Ftmp[0:5] = flx[idx+2*Nx::-Nx][0:5][::-1]
                    
                elif yidx == Ny-2: # Top boundary, 2 ghost cells
                    if xc < xs:
                        Utmp[-1] = U0[0]
                        Utmp[-2] = U0[0]
                        Ftmp[-1] = flx0[0]
                        Ftmp[-2] = flx0[0]
                    else:
                        Utmp[-1] = U0[-1]
                        Utmp[-2] = U0[-1]
                        Ftmp[-1] = flx0[-1]
                        Ftmp[-2] = flx0[-1]                        
                    
                    Utmp[0:4] = U[idx+Nx::-Nx][0:4][::-1]
                    Ftmp[0:4] = flx[idx+Nx::-Nx][0:4][::-1]
                    
                elif yidx == Ny - 1: # Top boundary, 3 ghost cells
                    if xc < xs:
                        Utmp[-1] = U0[0]
                        Utmp[-2] = U0[0]
                        Utmp[-3] = U0[0]
                        Ftmp[-1] = flx0[0]
                        Ftmp[-2] = flx0[0]
                        Ftmp[-3] = flx0[0]
                    else:
                        Utmp[-1] = U0[-1]
                        Utmp[-2] = U0[-1]
                        Utmp[-3] = U0[-1]
                        Ftmp[-1] = flx0[-1]
                        Ftmp[-2] = flx0[-1]
                        Ftmp[-3] = flx0[-1]
                    
                    Utmp[0:3] = U[idx::-Nx][0:3][::-1]
                    Ftmp[0:3] = flx[idx::-Nx][0:3][::-1]
                else:
                    Utmp[0:3] = U[idx::-Nx][0:3][::-1]
                    Utmp[3:6] = U[idx+Nx::Nx][0:3]
                    Ftmp[0:3] = flx[idx::-Nx][0:3][::-1]
                    Ftmp[3:6] = flx[idx+Nx::Nx][0:3]
                    
                
            # END IF ORDER == 5
            # TODO: ORDER 3
            
            # Perform the matrix multiplication
            V[idx,:,:] = np.matmul(Lh[idx],Utmp.T).T
            H[idx,:,:] = np.matmul(Lh[idx],Ftmp.T).T
            VLF[idx,:,:] = np.matmul(np.diag(alpha),V[idx,:,:].T).T    
        # END FOR NELEM LOOP
    # END DIRECTION CONDITIONAL
    return V,H,VLF
