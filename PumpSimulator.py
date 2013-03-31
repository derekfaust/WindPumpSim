#! /usr/bin/env python2

#Import required modules
from __future__ import division #Force float division
import numpy as np              #Array library
import pylab as plot            #Plotting functions
from math import *              #For trignometric functions
from sys import float_info      #For stall cutoff
#Import custom modules for the simulation
import statepredictor as sp     #Numerical Solution Interface
from pump import *              #Pump components

#Function ready for optimization
def waterpumped(r_crank, phase_offset):
    #Simulation Parameters
    timespan = [0,10]   #Beginning and ending times
    timesteps = 1000    #Number of timesteps
    initial_theta= 0    #Initial angle of turbine and pump
    initial_omega= 10  #Set the initial speed. Zero yeilds no torque
    
    #Turbine Parameters
    turbine_I = .1          #Turbine moment of inertia
    in_sprocket_teeth = 9   #Number of teeth on input sprocket
    out_sprocket_teeth = 70 #Number of teeth on output sprocket
    #Ratio of output torque to input torque
    gear_ratio = out_sprocket_teeth/in_sprocket_teeth
    #Measured torque curve [rad/s, Nm] from turbine
    omega_torque_curve = np.array([[0,  0],
                                   [56.55-float_info.epsilon, 0],
                                   [56.55, 1.01],
                                   [63.36, .94],
                                   [75.40, .79],
                                   [81.68, .72],
                                   [86.39, .65],
                                   [90.06, .58],
                                   [92.89, .51], 
                                   [99.48, .43],
                                   [103.67, .36],
                                   [107.23, .29],
                                   [110.48, .22],
                                   [114.25, .14],
                                   [122.00, .07],
                                   [122.00+float_info.epsilon, 0]])
    #Modify torque curve to convert to torque and omega at pump input
    omega_torque_curve[:,0] = omega_torque_curve[:,0]/gear_ratio
    omega_torque_curve[:,1] = omega_torque_curve[:,1]*gear_ratio
    
    #Piston Parameters
    crank_radius = r_crank
    rod_length = .254           #10 inches, Should be as long as will fit.
    piston_radius = .022        #.875in to meters.
    piston_mass = .1            #Estimated from volume and PVC density.
    initial_piston_angle = 0
    
    #Initial Conditions
    times = np.linspace(timespan[0],timespan[1],timesteps)
    initial_state = np.array([initial_theta,initial_omega, 0])
    
    #Create Objects (from pump module) to Simulate
    the_turbine = Turbine(turbine_I, omega_torque_curve)
    one_piston = Piston(crank_radius, rod_length, piston_radius, piston_mass,
                        initial_piston_angle)
    two_piston = Piston(crank_radius, rod_length, piston_radius, piston_mass,
                        initial_piston_angle+phase_offset)
    the_pump = Pump(the_turbine, [one_piston, two_piston])
    
    #Run the state prediction to get an array of states at requested times.
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
    
    #Return the amount of water pumped
    return states[-1,2]

#If this file is executed, do this:
if __name__ == '__main__':
    #Find the amount of water pumped at a given time.
    print waterpumped(.03, 0)
