#Import required modules
from __future__ import division
import numpy as np
from math import *

class Turbine:
    def __init__(self, I, performance):
        self.I = I
        self.omega = performance[:,0]
        self.torque = performance[:,1]
    
    def Tin(self, omega):
        return np.interp(omega, self.omega, self.torque)

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
    tube_height=1.5
    rho_w=1000
    max_pressure=rho_w*g*tube_height
    tube_area = tube_radius**2*pi
    state = [0,0]
    
    def __init__(self, drivesystem, mechanism):
        self.drive = drivesystem
        self.mechanism = mechanism    
    
    def backpressure(self, state, Q):
        vol_pumped = state[2]
        if Q > 0:
            backpressure = min(self.max_pressure, self.rho_w*self.g*(vol_pumped/self.tube_area))           
            return backpressure
        else:
            return 0

    def get_state(self):
        return self.state

    def get_torque(self, state):
        torque = self.mechanism.get_torque(state)
        return torque

    def net_power(self, state):
        bp = self.backpressure(state, self.mechanism.Q(state))
        mech_power = self.mechanism.power(state, bp, self.drive.Tin(state[1]))
        return self.drive.power(state) - mech_power

    def omegadot(self, state):
        pnet = self.net_power(state)
        c_omegadot_t = self.drive.de_omegadot_coef(state)
        c_omegadot_m = self.mechanism.de_omegadot_coef(state)
        b_de_m = self.mechanism.de_offset(state)
        omegadot = (pnet - b_de_m)/(c_omegadot_t+c_omegadot_m)
        return omegadot

    def statedot(self, t, state, p):
        #state = [theta, omega, vol]
        thetadot = state[1]
        omegadot = self.omegadot(state)
        voldot = self.mechanism.Q(state)
        statedot = [thetadot, omegadot, voldot]
        return statedot

class Piston:
    state = [0, 0]
    mu = .5
    
    def __init__(self, r_crank, r_rod, r_piston, m_piston, theta):
        self.r_c = r_crank
        self.r_r = r_rod
        self.mass = m_piston
        self.area = r_piston**2*pi
        self.theta0 = theta
    
    def alpha(self, theta):
        return asin(self.r_c/self.r_r*sin(theta))

    def xstate(self, state):
        return [xpos(state), xdot(state)]

    def xpos(self, state):
        theta = state[0]+self.theta0
        comp1 = self.r_c*cos(theta)
        comp2 = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        xpos = -comp1+comp2
        return xpos

    def xdot(self, state):
        theta = state[0]+self.theta0
        omega = state[1]
        xdot = -omega*self.xdotcoef(theta)
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

    def Q(self, state):
        xd = self.xdot(state)
        if xd > 0:
            q = self.xdot(state)*self.area
        else:
            q = 0
        return q
    
    def power(self, state, backpressure, torque):
        theta = state[0]+self.theta0
        pumpingpower = backpressure*self.area*self.xdot(state)
        T_rod = torque/(self.r_c*sin(pi-theta-self.alpha(theta)))
        power_lost_friction = abs(self.mu*T_rod*sin(theta)*self.xdot(state))
        power_lost_friction = 0
        power_out = pumpingpower + power_lost_friction
        return power_out

    def de_omegadot_coef(self, state):
        theta = state[0]+self.theta0
        omega = state[1]
        coef = self.mass*omega*self.xdotcoef(theta)**2
        return coef

    def de_offset(self, state):
        theta = state[0]+self.theta0
        omega = state[1]
        offset = self.mass*omega**3*self.xdotcoef(theta)*self.xddotcoef(theta)
        return offset
