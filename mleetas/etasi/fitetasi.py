""" Temporal ETASI model (With short-term incompletness) parameters inversion 
procedure. 

See Hainzl, S., 2021. ETAS-Approach Accounting for Short-Term Incompleteness of
Earthquake Catalogs. Bulletin of the Seismological Society of America. 
https://doi.org/10.1785/0120210146

Parameter space:
===============
    theta   = ( A, c, p, alpha, mu, b, Tb )

Author:
    Luc Moutote -- luc.moutote@gmail.com
"""

import numpy as np
import os
from scipy.optimize import minimize
from scipy.optimize import fmin_l_bfgs_b
from multiprocessing import Pool
from scipy import integrate
import pdb
import time as time
 
def nratio(theta,T,magmax):
    """
    Compute the branching ratio n for a catalog of length T and for a 
    troncated Gutenberg-Richter distribution. 

    The branching ratio here is define for the classic temporal ETAS model. It's 
    just to understand if current ETASI paramaters will allow to compute non-
    explosives synthtetics ETAS catalogs. 
    Reminder: To generate synthetics ETASI catalogs we first generate a classic 
    ETAS catalog, then we remove hidden event along Tb. 

    Parameters:
    ----------
    theta : tuple or list
        The ETASI parameters.
    T : float 
        Duration of the catalog in unit of time (Usually days).
    magmax : float
        Maximum magnitude of the catalog.

    Returns:
    --------
    n : float
        The branching ratio n
    """
    A, c, p, al, mu, b, Tb = theta                                                         
    beta    = b*np.log(10) 
    Naft    = A/(1-p)/(al-beta)*( np.exp( (al-beta)*magmax ) -1) *( (c+T)**(1-p) - c**(1-p))                  
    N       = -(1/beta) * ( np.exp(-beta*magmax) - 1)
    n       = Naft/N
    return n
 
def ll3multi(todo):
    """
    Compute LL3 for a set of samples (see Hainzl 2021). 

    LL3 = third member of the ETASI Log-Likelyhood function. The full ETASI 
    condtional intensity R(t) need to be numerically integrated to obtain LL3. 

    As this is long to compute over the full length of the catalog in one time,
    this funcion is called by mutliproceessing to integrate LL3 over several 
    batch of samples at the same time.

    Parameters:
    -----------
    todo : tuple
        A tuple containing (smartsample,tevent,mevent,theta)   
    smartsample : numpy array
        The time samples where to integrate LL3 (Need to be sorted)      
    tevent : numpy array
        Events origine time of the catalog (Need to be sorted)
    mevent : numpy array
        Events magnitude of the catalog 
    theta : tuple or list
        The ETASI parameters.

    Returns:
    --------
    LL3 : float
        LL3 numerically integrated over the time sample contained in 
        smartsample.
    """
    smartsample,tevent,mevent,theta = todo
    A, c, p, al, mu, b, Tb   = theta
    R       = smartsample*0
    for i,t in enumerate(smartsample):
        ii      = np.where(tevent < t)[0]
        R0      = mu + np.sum( A * np.exp(al*mevent[ii]) * (c + t -tevent[ii])**-p )
        R[i]    = (1/Tb) * (1 - np.exp(-Tb*R0))

    LL3 = integrate.trapz(R,smartsample)
    return LL3
 
def cost(theta,tevent,mevent,LB,UB,smartsample):
    """
    Compute ETASI Log-Likelyhood function LL (See Hainzl 2021). 
    LL  = LL12 (First and second member) + LL3 (Third member) 

    LL3 need to be numerically integrated and can be computationally expensive 
    for long catalogs.
    LL3 integration use multiprocessing (Pool).
    
    Parameters:
    -----------
    theta : tuple or list
        The ETASI parameters.
    tevent : numpy array
        Events origine time of the catalog (Need to be sorted)
    mevent : numpy array
        Events magnitude of the catalog 
    UB : tuple or list
        The upper bounds of the ETASI parameters
    LB : tuple or list
        The lower bounds of the ETASI parameters
    smartsample : numpy array
        The time samples where to integrate LL3 (Need to be sorted)      
    
    Returns:
    --------
    LL : float
        The Log-Likelyhood for the current ETASI parameter.
    """

    """ Convert parameters back to their true value for LLH computation """
    A, c, p, al, mu, b, Tb = theta
    # Log distributed between bounds
    A       = 10**(UB[0]*A  + LB[0]*(1-A))
    c       = 10**(UB[1]*c  + LB[1]*(1-c))
    mu      = 10**(UB[4]*mu + LB[4]*(1-mu))
    Tb      = 10**(UB[6]*Tb + LB[6]*(1-Tb))
    # Linearly distributed between bounds
    b       = UB[5]*b   + LB[5]*(1-b)
    p       = UB[2]*p   + LB[2]*(1-p)
    al      = UB[3]*al  + LB[3]*(1-al)

    theta   = A, c, p, al, mu, b, Tb

    """ Compute LLH 1st and 2nd member """
    beta    = b*np.log(10)
    LL1     = 0
    LL2     = 0
    for i,t in enumerate(tevent):
        R0      = mu + A*np.sum( np.exp(al*mevent[:i])*(c+t-tevent[:i])**(-p) )  
        LL1    += np.log( beta  * R0 * Tb * (np.exp(-beta*mevent[i]) * np.exp(- R0 * Tb * np.exp(-beta*mevent[i]))) / (1-np.exp(-R0*Tb))   )
        LL2    += np.log( (1/Tb) * (1 - np.exp(-Tb*R0)) )

    """ Compute LLH  3rd member"""
    # Multiproc intructions to feed the multiproc function
    ncores  = os.cpu_count()
    sublen  = int(np.ceil(len(smartsample)/ncores))
    todo    = []
    for n in range(ncores):
        istart  = n * sublen

        if (n+1) * sublen +1 <= len(smartsample):
            iend    = (n+1) * sublen + 1
        else:
            iend   = len(smartsample)

        todo.append((smartsample[istart:iend],tevent,mevent,theta))

    # Lauch multiproc
    pool    = Pool(ncores) #select the nbr of core to work on
    outs    = pool.map(ll3multi,todo)
    pool.close()
    pool.join()

    # Finally the integral 
    LL3     = np.sum(outs)

    """ Sum LLH members """
    # L-BFGS-B had some issure when LL12 is -np.inf 
    # from np.log(of something to 0)
    if LL1 == -np.inf or LL2 == -np.inf:
        LL = LL3 - (-1e10)
    else:    
        LL = LL3 - LL1 - LL2

    return LL



def fitETASI(tevent,mevent,theta):
    """
    Estimate ETASI parameters from a temporal catalog of seismicity. Following 
    Hainzl 2021 a Log-Likelyhood function (Eq.11) is maximized. We use the 
    L-BFGS-B algorithm to minimize '- Log-Likelyhood'.
    
    ETASI conditional intensity:
    ----------------------------
    R(t) = (1/Tb) * (1 - np.exp(-Tb*R0))
    With:
        R0 = mu + np.sum( A * np.exp(al*(m_i<t)) * (c + t - (t_i<t) )**-p )

    See Hainzl, S., 2021. ETAS-Approach Accounting for Short-Term i
    Incompleteness of Earthquake Catalogs. Bulletin of the Seismological 
    Society of America. https://doi.org/10.1785/0120210146

    Parameters:
    -----------
    tevent : numpy array
        Events origine time relative to 0 (tevent[0]=0)
    mevent : numpy array
        Magnitude of events. (WARNING) Mc must be = 0. Remember to remove your 
        Mc from your magnitudes before calling the function.
    theta : tuple or list
        Initial guess for the ETASI parameters (A,c,p,al,mu,b,Tb).
        A : float
            Average aftershock productivity
        c : float
            From Omori-Utsu law
        p : float
            From Omori-Utsu law
        al : float
            Scale the aftershock productivity of magnitude 
        mu : float
            Background rate
        b : float
            b-value of the Gutenberg-Richter law
        Tb : float
            Shadow time before one can observe an event t2 with a 
            magnitude m2<m1 ater an event t1 of magnitude m1.  

    Returns:
    --------
    theta : tuple
        The inverted ETASI parameter that minimize '- Log-likelyhood'.
    LLH : float
        The final value of '- Log-Likelyhood'.
    disp : dict
        Convergence messages from the L-BFGS-B solver.

    Notes:
    ------ 
    The third member of the cost function (Log likelyhood Eq.11) is 
    paralelized with Pool 
    """

    t1          = time.perf_counter()

    A, c, p, al, mu, b, Tb = theta
    tevent      = tevent - tevent[0]
    T           = tevent[-1]           
    magmax      = max(mevent) 
    n           = nratio(theta,T,magmax) 

    print('Inital guess:\n'
          '............. A  = {:.10f}\n'
          '............. c  = {:.10f}\n'
          '............. p  = {:.10f}\n'
          '............. al = {:.10f}\n'
          '............. mu = {:.10f}\n'
          '............. b  = {:.10f}\n'
          '............. Tb = {:.10f}\n'
          '............. n  = {:.10f}\n'.format(A,c,p,al,mu,b,Tb,n))

    # CONSTRAINTS ON PARAMETER SPACE
    min_A   ,max_A  = np.log10(1e-5)  , np.log10(1e1)     
    min_c   ,max_c  = np.log10(1e-5)  , np.log10(1e-1)    
    min_p   ,max_p  = 0.2             , 2.5              
    min_al  ,max_al = 0.2             , 3.5               
    min_mu  ,max_mu = np.log10(1e-4)  , np.log10(100)     
    min_b   ,max_b  = 0.2             , 2                 
    min_Tb  ,max_Tb = np.log10(1e-5)  , np.log10(1/24)    

    LB          = np.array([min_A,min_c,min_p,min_al,min_mu,min_b,min_Tb])
    UB          = np.array([max_A,max_c,max_p,max_al,max_mu,max_b,max_Tb])

    # Convertion of Values for the solver. Betwen [0 and 1]. Linearly distributed or Log distributed
    # Allow to sample each parameters between its bounds with a similar step
    A           = ( np.log10(A) - min_A )   / (max_A   -min_A )
    c           = ( np.log10(c) - min_c )   / (max_c   -min_c )
    p           = ( p           - min_p )   / (max_p   -min_p )
    al          = ( al          - min_al)   / (max_al  -min_al)
    mu          = ( np.log10(mu)- min_mu)   / (max_mu  -min_mu)
    b           = ( b           - min_b )   / (max_b   -min_b )
    Tb          = ( np.log10(Tb)- min_Tb)   / (max_Tb  -min_Tb)

    theta       = A, c, p, al, mu, b, Tb
    
    # Solve using L-BFGS-B
    # Define a smart sampling of R(t) for numerical integration of 3rd member of Eq.11
    smartsample = []  
    log         = np.array([0,1e-7,1e-3,3e-3,7e-3,1.5e-2,3e-2,5e-2,1.5e-1,.45,1.35]) 
    dt          = tevent[1:] - tevent[:-1]
    for i,t in enumerate(tevent[:-1]):
        ii      = np.where(log < dt[i])[0]
        tlog    = t + log[ii]
        smartsample.append(tlog)
    smartsample = np.concatenate(smartsample)
    
    bounds      = [(0,1)]*len(theta)    
    theta, LLH, disp = fmin_l_bfgs_b(cost                                           ,
                                     theta                                          ,
                                     fprime  = None                                 ,
                                     args    = (tevent,mevent,LB,UB,smartsample)    ,
                                     approx_grad = True                             ,
                                     bounds  = bounds                               ,
                                     m       = 10                                   ,
                                     factr   = 1e7                                  ,
                                     pgtol   = 1e-5                                 ,
                                     epsilon = 1e-4                                 ,
                                     iprint  = 1                                 ,
                                     maxfun  = 15000                                ,       
                                     maxiter = 15000                                ,  
                                     disp    = None                                 ,
                                     callback= None                                 ,
                                     maxls   = 10                       )

    A, c, p, al, mu, b, Tb = theta

    # Back to True Value 
    A       = 10**(UB[0]*A  + LB[0]*(1-A))
    c       = 10**(UB[1]*c  + LB[1]*(1-c))
    p       = UB[2]*p       + LB[2]*(1-p)
    al      = UB[3]*al      + LB[3]*(1-al)
    mu      = 10**(UB[4]*mu + LB[4]*(1-mu))
    b       = UB[5]*b       + LB[5]*(1-b)
    Tb      = 10**(UB[6]*Tb + LB[6]*(1-Tb))

    t2          = time.perf_counter()
    tCPUh       = (t2-t1)/60/60
    tCPUm       = (t2-t1)/60
    tCPUs       = (t2-t1)

    # Final solution
    theta       = A, c, p, al, mu, b, Tb

    # Branching ratio ETAS classic
    n           = nratio(theta,T,magmax) 

    print("Theta final:\n"
          "............ A  = {:.10f}\n"
          "............ c  = {:.10f}\n"
          "............ p  = {:.10f}\n"
          "............ al = {:.10f}\n"
          "............ mu = {:.10f}\n"
          "............ b  = {:.10f}\n"
          "............ Tb = {:.10f}\n"
          "............ n  = {:.10f}\n"
          "LLH Final:\n"
          "............ LLH = {:.5f}\n"
          "Full CPU time:\n"
          "............ Hours   = {:.2f}\n" 
          "............ Minutes = {:.2f}\n" 
          "............ Seconds = {:.2f}\n".format(A,c,p,al,mu,b,Tb,n,LLH,tCPUh,tCPUm,tCPUs))
    
    return theta, LLH, disp 
