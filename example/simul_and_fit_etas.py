"""
How to simulate and fit the temporal ETAS model 

References ETAS:
=================
    Zhuang, J., Harte, D., Werner, M.J., Hainzl, S., Zhou, S., 2012. Basic 
    models of seismicity: Temporal models. Community Online Resource for 
    Statistical Seismicity Analysis Theme V.

Luc Moutote - ITES - luc.moutote@gmail.com - 2023
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

from mleetas.simulation import simuletas
from mleetas.simulation import simulmag 
from mleetas.etas import fitetas


"""
Parameters
"""
# Duration of catalog (arbitrary unit; lets say days)
T   = 2000

# Magnitude of completness
mc  = 2

# ETAS parameter
A   = 0.01
c   = 0.01
p   = 1.1
al  = 2
mu  = 1

# Gutemberg Richter b-value; for magnitude distribution
b   = 1
mag = simulmag.gutricht(b,mc,10000)



"""
Generate a synthetic ETAS catalog
"""
tetas, metas = simuletas.etas(A,c,p,al,mu,T,mc,mag,lim=20000)
print('Number of event: ', len(tetas))

# Plot
plt.figure()
plt.plot(tetas,metas,'k.')
plt.xlabel('Time')
plt.ylabel('Magnitude')
plt.show()


"""
Fit the ETAS model on the generated catalog
We use the array theta to gather all ETAS parameters
theta = np.array([A,c,p,al,mu])
"""

# In the classic ETAS model the b-value is independant from ETAS parameters
# It can be estimated aside from ETAS with the aki estimator
b_aki       = simulmag.bvalue_aki(mag,mc)
print('MLE b-value: ',b_aki)

# Try to give a good starting point for the parameter search
theta0 = np.array([0.1,0.005,1.6,1.2,3])

# Run Etas fit, don't forget to remove mc from magnitudes
thetafit, LLH, disp = fitetas.fitETAS(tetas,metas-mc,theta0)








