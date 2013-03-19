#! /usr/bin/env python2

#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *
from sys import float_info
#Import custom modules
import statepredictor as sp
from pump import *

def waterpumped(r_crank, r_rod, r_piston):
    #Simulation Parameters
    timespan = [0,10]
    timesteps = 1000
    initial_theta= 0
    initial_omega= 8.1
    
    #Turbine Parameters
    turbine_I = .1
    in_sprocket_teeth = 9
    out_sprocket_teeth = 70
    gear_ratio = out_sprocket_teeth/in_sprocket_teeth
    omega_torque_curve = np.array([[0,  0],
                                   [42-float_info.epsilon, 0],
                                   [42, 1.28],
                                   [63, .64],
                                   [83, .32],
                                   [104, .16],
                                   [125, .08],
                                   [146, 0]])
    omega_torque_curve[:,0] = omega_torque_curve[:,0]/gear_ratio
    omega_torque_curve[:,1] = omega_torque_curve[:,1]*gear_ratio
    
    #Piston Parameters
    crank_radius = r_crank
    rod_length = r_rod
    piston_radius = r_piston
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
    
    return states[-1,2]

print waterpumped(.03, .3, .08)