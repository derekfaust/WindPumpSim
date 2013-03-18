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
initial_omega= 10

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
initial_state = np.array([initial_theta,initial_omega, 0])

#Create Objects to Simulate
the_turbine = Turbine(turbine_I, omega_torque_curve)
one_piston = Piston(crank_radius, rod_length, piston_radius, piston_mass,
                    initial_piston_angle)
the_pump = Pump(the_turbine, one_piston)

states = sp.predict(the_pump.statedot, times, initial_state, [])

#Plot the angular velocity and water pumped vs time
fig = plot.figure(1)
plot.subplot('211')
plot.plot(times, states[:,1], 'k')
plot.xlabel('Time (s)')
plot.ylabel('Omega (rad/s)')
plot.title('Pump')
plot.subplot('212')
plot.plot(times, states[:,2]*1000, 'k')
plot.xlabel('Time (s)')
plot.ylabel('Volume Pumped (L)')
plot.title('Pump')
fig.tight_layout()
plot.show()
