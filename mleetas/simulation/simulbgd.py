"""
Functions to generate background catalogs

Luc Moutote - ITES - lmoutote@unistra.fr
"""
import numpy as np
import random as rd

def poisson_stat(T,mu,mag):
    """ 
    Generated un synthetic catalog where the number of event over a time-window
    is distributed along a Poisson law. (i.e. waiting time between succesives
    events is distributed along an exponential law. Magnitude are randomly 
    sampled in a list. 
    
    Parameters:
    -----------
    T : float
        Duration of the catalog to generate 
    mu : float
        Background rate of event. Number of events per unit of time
    mag : numpy array
        A array containing magnitudes distribution that will be randomly 
        sampled during the generation process. For usual case, magnitudes 
        are distributed along Gutenberg-Richter law. Use "gutricht()" to 
        generate such distribution. 

    Returns:
    --------
    tday : numpy array
        Origin times of the simulated backgound catalog. (Time = ]0,T]).
    mag : numpy array
        Magnitudes of the simulated background catalog.
    """
     
    t    = [0]
    m    = [0]
        
    # While the last event (t,m) is < T   
    while True: 
        # Waiting time for the new event: Inverse transform exponential law
        dt      =-1/mu*np.log(rd.random())
        # Add the new waiting time to the last event
        newt    = t[-1] + dt

        # Add the new event to the list
        if newt < T:
            t.append(newt)
            m.append( rd.choice(mag) )

        # Case if the first dt is > T ( to have at least one event in the cat.)
        elif len(t)==1:
            continue

        # If last event (t,m) > T
        else:
            break
    
    return np.array(t[1:]),np.array(m[1:])


