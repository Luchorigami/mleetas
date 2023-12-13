""" Stochastic simulation of the temporal Epidemic Type Aftershock Sequence 
(ETAS) model and the ETASI model (ETAS but accounting fot short-term 
Incompleteness). Basically build synthetics earthquake catalog (time and 
magnitudes)

References ETASI:
=================
    Hainzl, S., 2021. ETAS-Approach Accounting for Short-Term Incompleteness of
    Earthquake Catalogs. Bulletin of the Seismological Society of America. 
    https://doi.org/10.1785/0120210146

References ETAS:
=================
    Zhuang, J., Harte, D., Werner, M.J., Hainzl, S., Zhou, S., 2012. Basic 
    models of seismicity: Temporal models. Community Online Resource for 
    Statistical Seismicity Analysis Theme V.

Author:
    Luc Moutote -- luc.moutote@gmail.com
"""

import numpy as np
import random as rd
import time

from mleetas.simulation.simulbgd import poisson_stat


def etasi(A,c,p,al,mu,Tb,T,mc,mag,lim=None,custom_bgd=None):
    """
    Generate a synthetic temporal ETASI catalog. First generate a classic ETAS
    catalog, then remove event hidden by short-term incompletness. A event 
    (ti,mi) is hidden if it is preceded by a larger magnitude in [t-Tb,t[ 
    (see Hainzl 2021).
    
    ETASI conditional intensity:
    ---------------------------
    R(t) = (1/Tb) * (1 - np.exp(-Tb*R0(t)))
    With:
        R0(t) = mu + np.sum( A * np.exp(al*(m_i<t)) * (c + t - (t_i<t) )**-p )

    Parameters:
    -----------
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
    Tb : float
        Shadow time before one can observe an event t2 with a 
        magnitude m2<m1 ater an event t1 of magnitude m1.  
    T : float
        Duration of the catalog to generate
    mc : float
        Magnitude of completness of the catalog
    mag : numpy array
        A array containing magnitudes distribution that will be randomly 
        sampled during the generation process. For usual case, magnitudes 
        are distributed along Gutenberg-Richter law. Use "gutricht()" to 
        generate such distribution. 
    lim : int 
        If the generated catalog is > lim, the function return [],[] Avoid 
        non-stationary cases or too productive catalog. Deactivated (None) by 
        default.
    custom_bgd : tuple
         A tuple ([t0,..,ti],[m0,...,mi]) containing background event. ETAS 
         aftershock cascades will be generated over this list. If None, the 
         background events are generated with "bgdpoisson()" following a 
         poisson distribution hypothesis.

    Returns:
    --------
    tetasI : numpy array
        Origin times of the simulated ETASI catalog. (Time = ]0,T]).
    metasI : numpy array
        Magnitude of the ETASI catalog.
    tetas : numpy array
        Origin times of the simulated ETAS catalog before the short-term 
        incompleteness step.
    metas : numpy array
        Magnitude of the simulated ETAS catalog before the short-term 
        incompleteness step.
    """
    print('ETASI simulation:')
    print(A,c,p,al,mu,Tb,mc)
    # Generate a classic temporal ETAS catalog
    tetas, metas    = etas(A,c,p,al,mu,T,mc,mag,lim=lim,custom_bgd=custom_bgd)

    if len(tetas):
        # Remove event hidden by short-term incompletness along Tb
        print('Building short term incompletness...')
        ihide   = []
        for i,t in enumerate(tetas):
            ii      = np.where( (tetas >= t-Tb) * (tetas< t) * (metas > metas[i]) )[0]
            if len(ii):
                ihide   .append(i)

        if len(ihide):
            ihide   = np.array(ihide)
            tetasI  = np.delete(tetas,ihide)
            metasI  = np.delete(metas,ihide)
            print(f'ETASI simulation terminated: Nb event = {len(tetasI)} -- Mag max = {np.max(metasI)}\n')     

        else:
            tetasI  = tetas
            metasI  = metas
            print(f'ETASI simulation found zero hidden event in the ETAS simulation')     
            print(f'ETASI simulation terminated: Nb event = {len(tetasI)} -- Mag max = {np.max(metasI)}\n')     
    
    # Case tetas is empty : return tetasI as empty 
    else:
        tetasI  = tetas
        metasI  = metas

    return tetasI, metasI, tetas, metas


def etas(A,c,p,al,mu,T,mc,mag,lim=None,custom_bgd=None):
    """
    Generate a synthetic temporal ETAS catalog. (see Zhuang et al, 2012)
    
    ETAS conditional intensity:
    ---------------------------
    R(t) = mu + np.sum( A * np.exp(al*(m_i<t)) * (c + t - (t_i<t) )**-p )
    
    Parameters:
    -----------
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
    T : float
        Duration of the catalog to generate
    mc : float
        Magnitude of completness of the catalog
    mag : numpy array
        A array containing magnitudes distribution that will be randomly 
        sampled during the generation process. For usual case, magnitudes 
        are distributed along Gutenberg-Richter law. Use "gutricht()" to 
        generate such distribution. 
    lim : int 
        If the generated catalog is > lim, the function return [],[] Avoid 
        non-stationary cases or too productive catalog. Deactivated (None) by 
        default.
    custom_bgd : tuple
         A tuple ([t0,..,ti],[m0,...,mi]) containing background event. ETAS 
         aftershock cascades will be generated over this list. If None, the 
         background events are generated with "bgdpoisson()" following a 
         poisson distribution hypothesis.

    Returns:
    --------
    tetas : numpy array
        Origin times of the simulated ETAS catalog.
    metas : numpy array
        Magnitude of the simulated ETAS catalog.
    """
    print('Simulating ETASI...')

    if p == 1:
        print('p=1 not implemented yet')
        return 0        
    else:
        None
           
    ''' Background = catalogue poiss '''
    if custom_bgd:
        t, m   = custom_bgd 
    else:
        t, m   = poisson_stat(T, mu, mag)
    

    ''' Cascade ETAS'''
    # Initialization of list with the backroung event
    tetas, metas    = t.tolist(), m.tolist()
    
    # Index of the list event for the while True
    count           = 0
    # Draw aftershocks for each event in the list 
    while True:
        tmain   = tetas[count]
        mmain   = metas[count]
        # Compute the theorical number of expected aftershocks in ]tmain:T[
        kmag   = A*np.exp(al*(mmain-mc)) / (1-p) * ( (T-tmain + c)**(1-p) - c**(1-p) )
        # Draw a stocastic number of aftershock in a Poisson law of mean kmag
        Naft    = np.random.poisson(kmag)
          
        # For each aftershock, sample a origin in the omori-utsu law
        # over ]tmain:T[
        coef    = (T-tmain + c)**(1-p) - c**(1-p)
        newt    = [tmain + ( coef*rd.random() + c**(1-p) )**(1/(1-p)) - c for _ in range(Naft)]
        # Sample a magnitude for each aftershock
        newm    = [rd.choice(mag) for _ in range(Naft)] 
     
        # Add new events in the list
        if len(newt):
            [tetas.append(nt) for nt in newt]
            [metas.append(nm) for nm in newm]
            
            # Sort events after tmain
            tetasnosort, metasnosort = tetas[:count], metas[:count]
            tetastosort, metastosort = tetas[count:], metas[count:]
            ordre       = np.argsort(tetastosort)
            tetastosort = np.array(tetastosort)[ordre]
            metastosort = np.array(metastosort)[ordre]            
            tetas       = np.concatenate( [tetasnosort,tetastosort])
            metas       = np.concatenate( [metasnosort,metastosort])
       
            tetas   = tetas.tolist()
            metas   = metas.tolist()

        else:
            None
            
        # Loading text, allow to understand if simulation converge or diverge
        print('ETAS... : Curent event: {:d}, Current duration: {:.2f}, Total length: {:d}, Magmax: {:.1f}, Difference: {:d}     \r'.format(count, tmain, len(tetas),max(metas),len(tetas)-count),end="")

        # Update count to the next event in list
        count = count + 1
        
        # Stop rule : If we are at the last event in the list
        if count >= len(tetas):
            tetas   = np.array(tetas) 
            metas   = np.array(metas) 
            print(f'\nETAS simulation (T) terminated: Nb event = {len(tetas)} -- Mag max = {np.max(metas)}\n')     
            break
        
        # Stop rule : If we the length of the catalog is above lim 
        elif lim != None and len(tetas)>lim:
            print('\nETAS simulation stopped : len(tetas) > lim\n')
            tetas, metas = [], []
            break
        
        # Otherwise we continue the afterhock generation for the next event in 
        # the list
        else:
            continue

    return tetas, metas 
    

 
