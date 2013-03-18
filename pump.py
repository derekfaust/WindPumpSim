#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *

class DummyTurbine:
    Idrive = 1
    TorqueCurve = np.array([[0,0],[20,5],[40,4],[60,3],[80,2],[100,1]])
    Aout = 2
    OutPeriod = 2
    def Tout(self, t):
        return self.Aout*(1+sin(t/self.OutPeriod))


class Turbine:
    def __init__(self, I, performance):
        self.I = I
        self.omega = performance[:,0]
        self.torque = performance[:,1]
    
    def Tin(self, omega):
        return np.interp(omega, self.TorqueCurve[:,0], self.TorqueCurve[:,1])

    def set_state(self, state):
        [theta, omega, omegadot] = state
        return true
    
    def get_state(self):
        return [theta, omega, omegadot]
    
    def power(self, state):
        return state[1]*self.Tin(state[1])
    
    def de_omegadot_coef(self, state):
        return self.I*state[1]

class Pump:
    g = 9.8
    tube_radius=.01
    tube_height=1
    rho_w=1000
    vol_pumped=0
    max_pressure=rho_w*g*tube_height
    tube_area = tube_radius**2*pi
    state = [0,0,0]
    
    def __init__(self, drivesystem, mechanism):
        self.drive = drivesystem
        self.mechanism = mechanism    
    
    def backpressure(self, Q):
        if Q > 0:
            return min(self.max_pressure, self.rho_w*self.g*(self.vol_pumped/self.tube_area))
        else:
            return 0

    def get_state(self):
        return state

    def get_torque(self, state):
        torque = mechanism.get_torque(state)
        return torque

    def net_power(self, state):
        bp = self.backpressure(mechanism.Q(state))
        mech_power = mechanism.power(state, bp)
        return drive.power(state) + mech_power

    def omegadot(self, state):
        pnet = self.netpower(state)
        c_omegadot_t = self.drive.de_omegadot_coeff(state)
        c_omegadot_m = self.mechanism.de_omegadot_coeff(state)
        b_de_m = self.mechanism.de_offset(state)
        omegadot = (pnet - b_de_m)/(c_omegadot_t+c_omegadot_m)
        return omegadot

    def statedot(t, state):
        thetadot = state[1]
        omegadot = self.omegadot(state)
        statedot = [thetadot, omegadot]
        return statedot

class Piston:
    state = [0, 0]
    mu = .5
    
    def __init__(self, r_crank, r_rod, r_piston, m_piston, theta):
        self.r_c = r_crank
        self.r_r = r_rod
        self.m_p = m_piston
        self.area = r_piston**2*pi
        self.theta0 = theta
    
    def alpha(self, theta):
        return asin(self.r_c/self.r_r*sin(theta))

    def xstate(self, state):
        return [xpos(state), xdot(state)]

    def xpos(self, state):
        theta = state[0]
        comp1 = self.r_c*cos(theta)
        comp2 = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        xpos = comp1+comp2
        return xpos

    def xdot(self, state):
        theta = state[0]
        omega = state[1]
        xdot = omega*self.xdotcoef(theta)
        return xdot
    
    def xdotcoef(self, theta):
        comp1 = self.r_c*sin(theta)
        denom = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        comp2 = (self.r_c**2*sin(theta)*cos(theta))/denom
        return comp1-comp2

    def xddotcoef(self, theta):
        num = (self.r_c**2*cos(theta)*sin(theta))**2
        denom = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        coef = (self.r_c**2/denom+num/(denom**3)+self.r_c*cos(theta))
        return coef
    
    def power(self, state, backpressure):
        pumpingpower = backpressure*self.area
        return pumpingpower

    def de_omegadot_coeff(self, state):
        coef = self.m*state[1]*xdotcoef(state[0])**2
        return coef

    def de_offset(self, state):
        offset = self.m*state[1]**3*xdotcoef(state[0])*xddotcoef(state[0])
        return offset
