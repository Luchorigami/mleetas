'''
Temporal ETAS model parameters inversion procedure.

Parameter space:
    theta   = ( A, c, p, alpha, mu, b )

Author:
    Luc Moutote -- luc.moutote@gmail.com -- 30/11/23
'''

import numpy as np
import os
from scipy.optimize import minimize
from scipy.optimize import fmin_l_bfgs_b
from multiprocessing import Pool
from scipy import integrate
import pdb
import time as time
            
def cost(theta,tevent,mevent,LB,UB):
    """
    Compute ETASI Log Likelyhood function
    """

    """ Convert parameters back to their true value for LLH computation """
    A, c, p, al, mu = theta

    # Log distributed between bounds
    A       = 10**(UB[0]*A  + LB[0]*(1-A))
    c       = 10**(UB[1]*c  + LB[1]*(1-c))
    mu      = 10**(UB[4]*mu + LB[4]*(1-mu))
    # Linearly distributed between bounds
    p       = UB[2]*p   + LB[2]*(1-p)
    al      = UB[3]*al  + LB[3]*(1-al)

    theta   = A, c, p, al, mu

    """ Compute LLH 1st and 2nd member """
    LL1    = 0
    for i,t in enumerate(tevent):
        R0      = mu + A*np.sum( np.exp(al*mevent[:i])*(c+t-tevent[:i])**(-p) )  
        LL1    += np.log( R0 )

    """ Compute LLH  3rd member"""
    # Finally the integral 
    LL2     = mu*tevent[-1] + np.sum(A/(1-p)*np.exp(al*mevent)*((tevent[-1]+c-tevent)**(1-p)-c**(1-p)))

    """ Compute LLGR """
    b       = np.log10(np.e)/np.mean(mevent)
    LLGR    = len(tevent)*np.log(np.log(10)*b)-np.log(10)*b*np.sum(mevent)

    """ Sum LLH members """
    # L-BFGS-B had some issure when LL12 is -np.inf 
    # from np.log(of something to 0)
    if LL1 == -np.inf:
        LL = LL3 - (-1e10) - LLGR
    else:    
        LL = LL2 - LL1 - LLGR

    return LL



def fitETAS(tevent,mevent,theta):
    """
    Estimate ETAS parameter from a temporal catalog. A Log-Likelyhood function
    is maximized with the L-BFGS-B algorithm.
    The conditional intensity function of the ETAS model follows:

    R(t) = mu + np.sum( A * np.exp(al*(m_i<t)) * (c + t - (t_i<t) )**-p )
    t is the time.


    tevent  : List of events origine time relative to 0 (tevent[0]=0)
    mevent  : Related magnitude of events. Warning Mc must be = 0. Remember to 
              remove Mc from your magnitude before calling the function
    theta   : ETASI parameter (A,c,p,al,mu,b,Tb)
              A     : Average aftershock productivity
              c     : From Omori-Utsu law
              p     : From Omori-Utsu law
              al    : Scale the aftershock productivity of magnitude 
              mu    : Background rate
              b     : From Gutemberg Richter law
                      magnitude m2<m1 ater an event t1 of magnitude m1.  

    """

    t1          = time.perf_counter()

    A, c, p, al, mu = theta
    tevent      = tevent-tevent[0]

    print('Inital guess:\n'
          '............. A  = {:.10f}\n'
          '............. c  = {:.10f}\n'
          '............. p  = {:.10f}\n'
          '............. al = {:.10f}\n'
          '............. mu = {:.10f}\n'.format(A,c,p,al,mu))

    # CONSTRAINTS ON PARAMETER SPACE
    min_A   ,max_A  = np.log10(1e-5)  , np.log10(1e1)     
    min_c   ,max_c  = np.log10(1e-5)  , np.log10(1e1)    
    min_p   ,max_p  = 0.2             , 2.5              
    min_mu  ,max_mu = np.log10(1e-4)  , np.log10(100)     
    min_al  ,max_al = 0.2       , 3.5               

    LB          = np.array([min_A,min_c,min_p,min_al,min_mu])
    UB          = np.array([max_A,max_c,max_p,max_al,max_mu])

    # Convertion of Values for the solver. Betwen [0 and 1]. Linearly distributed or Log distributed
    # Allow to sample each parameters between its bounds with a similar step
    A           = ( np.log10(A) - min_A )   / (max_A   -min_A )
    c           = ( np.log10(c) - min_c )   / (max_c   -min_c )
    p           = ( p           - min_p )   / (max_p   -min_p )
    al          = ( al          - min_al)   / (max_al  -min_al)
    mu          = ( np.log10(mu)- min_mu)   / (max_mu  -min_mu)

    theta       = A, c, p, al, mu
    
    # Solve using L-BFGS-B
    # Define a smart sampling of R(t) for numerical integration of 3rd member of Eq.11

    bounds      = [(0,1)]*len(theta)    

    theta, LLH, disp = fmin_l_bfgs_b(cost                                           ,
                                     theta                                          ,
                                     fprime  = None                                 ,
                                     args    = (tevent,mevent,LB,UB)                ,
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

    A, c, p, al, mu = theta

    # Back to True Value 
    A       = 10**(UB[0]*A  + LB[0]*(1-A))
    c       = 10**(UB[1]*c  + LB[1]*(1-c))
    p       = UB[2]*p       + LB[2]*(1-p)
    al      = UB[3]*al      + LB[3]*(1-al)
    mu      = 10**(UB[4]*mu + LB[4]*(1-mu))

    t2      = time.perf_counter()
    tCPUh   = (t2-t1)/60/60
    tCPUm   = (t2-t1)/60
    tCPUs   = (t2-t1)

    theta   = A, c, p, al, mu

    print("Theta final:\n"
          "............ A  = {:.10f}\n"
          "............ c  = {:.10f}\n"
          "............ p  = {:.10f}\n"
          "............ al = {:.10f}\n"
          "............ mu = {:.10f}\n"
          "LLH Final:\n"
          "............ LLH = {:.5f}\n"
          "Full CPU time:\n"
          "............ Hours   = {:.2f}\n" 
          "............ Minutes = {:.2f}\n" 
          "............ Seconds = {:.2f}\n".format(A,c,p,al,mu,LLH,tCPUh,tCPUm,tCPUs))

    
    
    return theta, LLH, disp 
