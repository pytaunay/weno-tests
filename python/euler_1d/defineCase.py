""" File: defineCase.py
Description: defines a test case

"""

import numpy as np

""" Function: defineCase
Returns an array that defines the case in the form:
    [rhoL uL PL]
    [rhoR uR PR]
"""
def defineCase(number):
    left = np.zeros(3)
    right = np.zeros(3)
    
    # Sod's problem
    if number == 1:
        left = np.array([1,0,1])
        right = np.array([0.125,0,0.1])
        cfl = 0.9
        tmax = 0.1
        
    # Sod's problem 2 JCP 27:1 1978
    elif number == 7:
        left = np.array([1,0.75,1])
        right = np.array([0.125,0,0.1])    
        cfl = 0.9
        tmax = 0.17
    
    # Lax problem 
    elif number == 8:
        left = np.array([0.445,0.698,3.528])
        right = np.array([0.5,0,0.571])    
        cfl = 0.9
        tmax = 0.15
        
    # Mach 3 test case
    elif number == 9:
        left = np.array([3.857,0.92,10.333])
        right = np.array([1,3.55,1])    
        cfl = 0.9
        tmax = 0.09        
            
    return left,right,cfl,tmax