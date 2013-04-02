#Import required modules
from __future__ import division
import numpy as np
from math import *
from scipy.integrate import ode

def predict(statedot, times, initial_state, parameters):
    #Predict the state at a series of times given initial state
    #and a differential function.
    
    #Initialize the state array to save to
    state_array = np.zeros((np.size(times), np.size(initial_state)))
    #Create ODE object
    eqn = ode(statedot).set_integrator('dopri5')
    #Set up the ODE for integration
    eqn.set_initial_value(initial_state, times[0]).set_f_params(parameters)
    #Initialize the first value in the state array
    state_array[0,:] = initial_state;
    
    #Iterate for each time in times and record the state to the state array
    for index in range(1,np.size(times)):
        time = times[index]             #Get the current time
        eqn.integrate(time)             #Iterate at this time
        state_array[index,:] = eqn.y    #Get the state
    return state_array
