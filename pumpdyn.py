#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *
from scipy.integrate import ode

def solveODE(f, tarray, z0, p):
    print np.size(tarray)
    print np.size(z0)
    zarray = np.zeros((np.size(tarray), np.size(z0)))
    eqn = ode(f).set_integrator('vode')
    eqn.set_initial_value(z0, tarray[0]).set_f_params(p)
    zarray[0,:] = z0;
    idx=0
    for time in np.nditer(tarray[1:]):
        idx = idx+1
        eqn.integrate(time)
        zarray[idx,:] = eqn.y
    return zarray

def rhs(t, z, p):
    omega = z[0]
    KE = z[1]
    P = (p.Tin(omega)-p.Tout(t))*omega
    #print p.Tin(omega)
    #print t, omega, p.Tin(omega)
    if (KE <= 0) and (P<=0):
        zdot = np.array([-omega,-KE])
    elif (KE <=0):
        zdot = np.array([-omega,P])
    else:
        zdot = np.array([sqrt(2*KE/p.Idrive)-omega,P])
    return zdot

