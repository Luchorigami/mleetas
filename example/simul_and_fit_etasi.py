"""
How to simulate and fit the temporal ETASI model of Hainzl (2021)
The ETASI model can take into account the short term incpletness when the 
seismicity rate is high.

See Hainzl, S., 2021. ETAS-Approach Accounting for Short-Term Incompleteness of
Earthquake Catalogs. Bulletin of the Seismological Society of America. 
https://doi.org/10.1785/0120210146

Luc Moutote - ITES - luc.moutote@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

from mleetas.simulation import simuletas
from mleetas.simulation import simulmag 
from mleetas.etasi import fitetasi


"""
Parameters
"""
# Duration of catalog (arbitrary unit; lets say days)
T   = 2000

# Magnitude of completness
mc  = 2

# ETAS parameter
A   = 0.01
c   = 0.05
p   = 1.1
al  = 2
mu  = 1

# Blind time paramter
Tb  = 0.01

# Gutemberg Richter b-value; for magnitude distribution
b   = 1
mag = simulmag.gutricht(b,mc,10000)



"""
Generate a synthetic ETASI catalog
- (tetas  ,metas ) is a classic ETAS simulation before the incompletness step
- (tetasi ,metasi) is the corresponding ETASI simulation after removing hidden 
events from (tetas,metas).
"""
tetasi,metasi,tetas, metas = simuletas.etasi(A,c,p,al,mu,Tb,T,mc,mag,lim=20000)

# Plot
plt.figure()
plt.plot(tetas,metas,'r.',label='Events missed by incompletness')
plt.plot(tetasi,metasi,'k.')
plt.xlabel('Time')
plt.ylabel('Magnitude')
plt.show()

print('Number of event: ', len(tetas))



"""
Fit the ETASI model on the generated catalog

We use the array theta to gather all ETASI parameters (7)
Compared to ETAS we add 2 new parameter b, and Tb. 
b is the bvalue of the Gutenberg-Richter law
Tb is the blind time that model incompletness (see hainzl 2021)
We also add the b-value because in ETASI the magnitude distribution
and the other ETASI parameters are not independant. 

theta = np.array([A,c,p,al,mu,b,Tb])
"""
# Try to give a good starting point for the parameter search
theta0 = np.array([0.1,0.005,1.6,1.2,3,1.2,.1])

# Run Etasi fit, don't forget to remove mc from magnitudes
thetafit, LLH, disp = fitetasi.fitETASI(tetasi,metasi-mc,theta0)








