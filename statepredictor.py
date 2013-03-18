#Import required modules
from __future__ import division
import numpy as np
from math import *
from scipy.integrate import ode

def predict(statedot, times, initial_state, parameters):
    print np.size(times)
    print np.size(initial_state)
    print initial_state
    state_array = np.zeros((np.size(times), np.size(initial_state)))
    eqn = ode(statedot).set_integrator('vode')
    eqn.set_initial_value(initial_state, times[0]).set_f_params(parameters)
    state_array[0,:] = initial_state;
    for index in range(1,np.size(times)):
        time = times[index]
        eqn.integrate(time)
        state_array[index,:] = eqn.y
    return state_array
