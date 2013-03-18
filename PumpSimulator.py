#! /usr/bin/env python2

#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *
#Import custom modules
import statepredictor as sp
from pump import *

#Simulation Parameters
timespan = [0,100]
timesteps = 1000
initial_theta= 0
initial_omega= 20

#Turbine Parameters
turbine_I = 1
omega_torque_curve = np.array([[0,  0],
                               [20, 5],
                               [40, 4],
                               [60, 3],
                               [80, 2],
                               [100,1],
                               [120,0]])
#Piston Parameters
crank_radius = .03
rod_length = .3
piston_radius = .03
piston_mass = .1
initial_piston_angle = 0

#Initial Conditions
times = np.linspace(timespan[0],timespan[1],timesteps)
initial_state = np.array([initial_theta,initial_omega])

#Create Objects to Simulate
the_turbine = Turbine(turbine_I, omega_torque_curve)
one_piston = Piston(crank_radius, rod_length, piston_radius, piston_mass,
                    initial_piston_angle)
the_pump = Pump(the_turbine, one_piston)

states = sp.predict(the_pump.statedot, times, initial_state, [])

#Plot the velocity vs time
fig = plot.figure(1)
plot.plot(times, states[:,1], 'k')
plot.xlabel('Time')
plot.ylabel('Omega')
plot.title('Pump')
plot.show()
