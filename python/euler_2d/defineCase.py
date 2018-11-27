""" File: defineCase.py
Description: defines a test case

"""

import numpy as np

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

    # Double mach reflection 
    if case == 'double-mach':
        xlim = np.array([0.0,4.0])
        ylim = np.array([0.0,1.0])
        grid = np.array([260,80])
        tmax = 0.2
        cfl = 0.6

        options['xshock'] = 1/6
        options['angleshock'] = 60*np.pi/180.
        options['Upre'] = np.array([1.4,0,0,1]) # U BEFORE shock
        options['Upost'] = np.array([8,4.125*np.sqrt(3),-4.125,116.5]) # U AFTER shock

    # Flow over forward facing step
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
