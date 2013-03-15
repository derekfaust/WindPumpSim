# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 20:55:42 2013

@author: dfaust
"""

#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *
from scipy.integrate import ode

class Turbine:
    Idrive = 1
    TorqueCurve = np.array([[0,0],[20,5],[40,4],[60,3],[80,2],[100,1]])
    Aout = 2
    OutPeriod = 2
    def Tin(self, omega):
        return np.interp(omega, self.TorqueCurve[:,0], self.TorqueCurve[:,1])
    def Tout(self, t):
        return self.Aout*(1+sin(t/self.OutPeriod))

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

myTurbine = Turbine()
time = np.linspace(0,200,1000)
initSpeed = 50
initCond = np.array([initSpeed,.5*myTurbine.Idrive*initSpeed**2])
zarray = solveODE(rhs, time, initCond, myTurbine)

#Plot
fig = plot.figure(1)
plot.plot(time, zarray[:,0], 'r')
#plot.axis('equal')
plot.xlabel('X Position (m)')
plot.ylabel('Time')
plot.title('KE and Omega With Time')
#plot.savefig('12-1-10d.png')
plot.show()