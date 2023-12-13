"""
Function to generate magnitude distributions

Author:
    Luc Moutote -- luc.moutote@gmail.com -- 30/11/23
"""
import numpy as np
import random as rd

def gutricht(b,mc,N):
    """
    Sample magnitudes from the Gutemberg Richter law. 

    Magnitude are sample using inverse transform approach:
    ------------------------------------------------------
        f(m) = bln(10)exp(-bln(10)(m-mc))   [Gutemberg-Richter]
        =>  P(M<=m) = exp(-bln(10)(m-mc))
        =>  m       =-1/beta*ln(P(M<=m))
    
    Parameters:
    -----------
    b : float
        b-value of the Gutemberg Richter law
    mc : float
        Magnitude of completness
    N : int
        Number of magnitudes to sample (i.e. length of the output array)

    Returns:
    --------
    m : numpy array
        An array containing N magnitudes drawn from the Gutenberg-Richter law
    """
    beta    = b*np.log(10)
    m       = [-1/beta*np.log(rd.random()) + mc for _ in range(N)]
    
    return np.array(m)

def bvalue_aki(mag,mc):
    """
    Maximum likely-hood estimator of the b-value of the Gutenberg Richter law

    See a paper from Keitii Aki in the 60s

    Parameters:
    -----------
    mag: numpy array
        array of magnitudes of the catalogue (Only magnitudes 
        above the magnitude of completnes)

    Returns:
    --------
    b : float
        The maximum likely hood estimate of the G-R b-value for the given 
        magnitude distribution
    mc : float
        Magnitude of completness
    """
    return np.log10(np.e)/np.mean(mag-mc)


