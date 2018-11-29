""" File: defineCase.py
Description: defines a test case

"""

import numpy as np
from utils import rhoE_from_Pv

""" Function: defineCase
Returns an array that defines the case in the form:
    [rhoL uL PL]
    [rhoR uR PR]
"""
def defineCase(case):
    left = np.zeros(3)
    right = np.zeros(3)

    # Contains xmin, xmax, ymin, ymax
    xlim = np.zeros(2)     
    ylim = np.zeros(2)
    grid = np.zeros(2)
    
    options = {}

    # Double mach reflection 
    if case == 'double-mach':
        xlim = np.array([0.0,4.0])
        ylim = np.array([0.0,1.0])
        grid = np.array([256,64])
        tmax = 0.2
        cfl = 0.6

        ### Initial conditions
        # Post-shock
        q1 = 8
        u = 4.125*np.sqrt(3)
        v = -4.125
        P = 116.5
        
        q2 = q1*u
        q3 = q1*v
        q4 = rhoE_from_Pv(P,q1,u,v)
        options['Upost'] = np.array([q1,q2,q3,q4]) # U AFTER shock        
        
        # Pre-shock
        q1 = 1.4
        u = 0
        v = 0
        P = 1
        q2 = q1*u
        q3 = q1*v
        q4 = rhoE_from_Pv(P,q1,u,v)  
        options['Upre'] = np.array([q1,q2,q3,q4]) # U BEFORE shock

        options['xshock'] = 1/6
        options['angleshock'] = 60*np.pi/180.
        
    # TODO: Flow over forward facing step
    elif case == 'facing-step':
        xlim = np.array([0.0,3.0])
        ylim = np.array([0.0,1.0])
        grid = np.array([242,79])

        options['xystep'] = np.array([0.6,0.2]) # Location of step 

            
    options['case'] = case
    options['cfl'] = cfl
    options['tmax'] = tmax
    options['xlim'] = xlim
    options['ylim'] = ylim
    options['grid'] = grid

    return options 
