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
    theta = 0
    omega = 0
    omegadot = 0
    def __init__(self, I, performance, output):
        self.I = I
        self.omega = performance[:,1]
        self.torque = performance[:,2]
        self.output = output
    def Tin(self, omega):
        return np.interp(omega, self.TorqueCurve[:,0], self.TorqueCurve[:,1])
    #Pump is in control, no need to get pump torque.
    def Tout(self, time):
        return output.get_torque(time)
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
    tube_area = tube_radius**2*math.pi
    state = [0,0,0]
    
    def __init__(self, drivesystem):
        self.drive = drivesystem
        mechanism = Piston(self, .03, .3, .03, .25, 0)
    
    def backpressure(self, Q):
        if Q > 0:
            return min(self.max_pressure, self.rho_w*self.g*(self.vol_pumped/self.tube_area)
        else
            return 0

    def get_state(self):
        return state
    
    #DO NOT USE
    #State-based, not time-based, for now.
    def make_state(self, time):
        if lastupdate == time:
            return false
        else
            [self.theta, self.omega, self.omegadot] = self.drive.get_state()
            self.theta += self.theta0
            self.lastupdate = time
            return true

    def get_torque(self, state):
        torque = mechanism.get_torque(state)
        return torque

    def net_power(self, state):
        return drive.power(state) + mechanism.power(state)

    def omegadot(self, state):
        pnet = self.netpower(state)
        c_omegadot_t = self.drive.de_omegadot_coeff(state)
        c_omegadot_m = self.mechanism.de_omegadot_coeff(state)
        b_de_m = self.mechanism.de_offset(state)
        omegadot = (pnet - b_de_m)/(c_omegadot_t+c_omegadot_m)
        return omegadot

    def delta_state(t, state):
        thetadot = state[1]
        omegadot = pump.omegadot(state)
        statedot = [thetadot, omegadot]
        return statedot

class Piston:
    state [0, 0]
    lastupdate = 0
    mu = .5
    
    def __init__(self, pump, r_crank, r_rod, r_piston, m_piston, theta)
        masterpump = pump
        self.r_c = r_crank
        self.r_r = r_rod
        self.m_p = m_piston
        self.area = r_piston**2*pi
        self.theta0 = theta
        #make_state()
    
    #DO NOT USE
    #State-based, not time-based, for now.
    def make_state(self, time)
        if lastupdate == time:
            return false
        else
            [self.theta, self.omega, self.omegadot] = self.masterpump.get_state()
            self.theta += self.theta0
            self.alpha = asin(self.r_c/self.r_r*sin(theta))
            self.xpos = xpos(self, theta, omega, omegadot)
            self.xdot = xdot(self, theta, omega, omegadot)
            #Incalculable
            #self.xddot = xddot(self, theta, omega, omegadot)
            self.lastupdate = time
            return true

    def alpha(self, theta):
        return asin(self.r_c/self.r_r*sin(theta))

    def xstate(self, state):
        return [xpos(state), xdot(state)]

    def xpos(self, state):
        theta = state[0]
        omega = state[1]
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
    
    #BROKEN
    def get_torque(self, time):
        make_state(time)
        Q = self.area*self.xdot
        P = masterpump.get_pressure(Q)
        #xddot is incalculable
        F_rod = (self.m_p*self.xddot + P*self.area)/(
                cos(self.alpha)+self.mu*sin(self.alpha)*(self.xdot/abs(self.xdot)))
        torque = -F_rod*self.r_c*sin(pi-self.alpha-self.theta)
        return torque

    def power(self, state):
        x = self.xpos(self, state)
        v = self.xdot(self, state)
        pumpingpower = masterpump.backpressure(self.area*v)
        return pumpingpower

    def de_omegadot_coeff(self, state):
        coef = self.m*state[1]*xdotcoef(state[0])**2
        return coef

    def de_offset(self, state):
        offset = self.m*state[1]**3*xdotcoef(state[0])*xddotcoef(state[0])
        return offset
