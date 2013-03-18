#! /usr/bin/env python2

#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *
#Import custom modules
import statepredictor as sp
from pump import *

#Initial Conditions
times = np.linspace(0,100,1000)
initial_state = np.array([0,50])

#Create Objects to Simulate
the_turbine = Turbine(1, np.array([[0,0],[20,5],[40,4],[60,3],[80,2],[100,1]]))
one_piston = Piston(.03, .3, .03, .25, 0)
the_pump = Pump(the_turbine, one_piston)

states = sp.predict(the_pump.statedot, times, initial_state, [])

#Plot the velocity vs time
fig = plot.figure(1)
plot.plot(times, states[:,1], 'k')
plot.xlabel('Time')
plot.ylabel('Omega')
plot.title('Pump')
plot.show()
