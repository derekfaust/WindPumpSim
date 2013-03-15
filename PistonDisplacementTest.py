# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 21:36:48 2013

@author: dfaust
"""

#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *

class piston:
    r_c = .01
    r_r = .5

def xpos(theta, omega, omegadot, piston):
    comp1 = piston.r_c*cos(theta)
    comp2 = piston.r_r*sqrt(1-(piston.r_c/piston.r_r*sin(theta))**2)
    xpos = comp1+comp2
    return xpos

def xdot(theta, omega, omegadot, piston):
    comp1 = piston.r_c*sin(theta)
    denom = piston.r_r*sqrt(1-(piston.r_c/piston.r_r*sin(theta))**2)
    comp2 = (piston.r_c**2*sin(theta)*cos(theta))/denom
    xdot = omega*(comp1-comp2)
    return xdot
    
def xddot(theta, omega, omegadot, piston):
    denom = piston.r_r*sqrt(1-(piston.r_c/piston.r_r*sin(theta))**2)
    comp1 = piston.r_c*sin(theta)
    comp2 = piston.r_c**2*sin(theta)*cos(theta)
    omegadot_term = omegadot*(comp1-comp2/denom)
    num = (piston.r_c**2*cos(theta)*sin(theta))**2
    omega_term = omega**2*(piston.r_c/denom+num/denom**3+piston.r_c*cos(theta))
    xddot = omegadot_term + omega_term    
    return xddot

tFinal = 100
dt = .1
steps = int(floor(tFinal/dt))       #Find Number of Times
time = np.linspace(0,tFinal,steps)  #Create Time array
motion = np.zeros((1000, 3))
mypiston = piston()
theta = 0
omega = 1
for index in xrange(1,steps):
    t = time[index]
    omegadot = sin(t/2)
    theta = theta + omega*dt
    omega = omega + omegadot*dt
    motion[index,:] = np.array([xpos(theta,omega,omegadot,mypiston),xdot(theta,omega,omegadot,mypiston),xddot(theta,omega,omegadot,mypiston)])

#Plot
fig = plot.figure(1)
plot.plot(time, motion[:,0], 'r', time, motion[:,1], 'g')
#plot.axis('equal')
plot.xlabel('Time')
plot.ylabel('Position')
plot.title('Position With Time')
#plot.savefig('12-1-10d.png')
plot.show()