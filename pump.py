#Import required modules
from __future__ import division
import numpy as np
import pylab as plot
from math import *

class Turbine:
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
        # return self?
    def Tin(self, omega):
        return np.interp(omega, self.TorqueCurve[:,0], self.TorqueCurve[:,1])
    def Tout(self, time):
        return output.get_torque(time)
    def set_state(self, time, state):
        [theta, omega, omegadot] = state
        return true
    def get_state(self,time)
        return [theta, omega, omegadot]

class Pump:
    g = 9.8
    tube_radius=.01
    tube_height=1
    rho_w=1000
    vol_pumped=0
    max_pressure=rho_w*g*tube_height
    tube_area = tube_radius**2*math.pi
    theta = 0
    omega = 0
    omegadot = 0
    def __init__(self, drivesystem):
        self.drive = drivesystem
        mechanism = Piston(self, .03, .3, .03, .25, 0)
    def get_pressure(self, Q, Qdot):
        if Q > 0:
            return min(self.max_pressure, self.rho_w*self.g*(self.vol_pumped/self.tube_area)
        else
            return 0
    def get_state(self):
        return [self.theta, self.omega, self.omegadot]
    def make_state(self, time):
        if lastupdate == time:
            return false
        else
            [self.theta, self.omega, self.omegadot] = self.drive.get_state()
            self.theta += self.theta0
            self.lastupdate = time
            return true
    def get_torque(self, time):
        make_state(time)
        torque = mechanism.get_torque(time)
        return torque

class Piston:
    theta = 0
    omega = 0
    omegadot = 0
    alpha = 0
    xpos = 0
    xdot = 0
    xddot = 0
    lastupdate = 0
    mu = .5
    def __init__(self, pump, r_crank, r_rod, r_piston, m_piston, theta)
        masterpump = pump
        self.r_c = r_crank
        self.r_r = r_rod
        self.m_p = m_piston
        self.area = r_piston**2*pi
        self.theta0 = theta
        make_state()

    def make_state(self, time)
        if lastupdate == time:
            return false
        else
            [self.theta, self.omega, self.omegadot] = self.masterpump.get_state()
            self.theta += self.theta0
            self.alpha = asin(self.r_c/self.r_r*sin(theta))
            self.xpos = xpos(self, theta, omega, omegadot)
            self.xdot = xdot(self, theta, omega, omegadot)
            self.xddot = xddot(self, theta, omega, omegadot)
            self.lastupdate = time
            return true

    def xpos(self, theta, omega, omegadot):
        comp1 = self.r_c*cos(theta)
        comp2 = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        xpos = comp1+comp2
        return xpos

    def xdot(self, theta, omega, omegadot):
        comp1 = self.r_c*sin(theta)
        denom = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        comp2 = (self.r_c**2*sin(theta)*cos(theta))/denom
        xdot = omega*(comp1-comp2)
        return xdot
    
    def xddot(self, theta, omega, omegadot):
        denom = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        comp1 = self.r_c*sin(theta)
        comp2 = self.r_c**2*sin(theta)*cos(theta)
        omegadot_term = omegadot*(comp1-comp2/denom)
        num = (self.r_c**2*cos(theta)*sin(theta))**2
        omega_term = omega**2*(self.r_c/denom+num/denom**3+self.r_c*cos(theta))
        xddot = omegadot_term + omega_term    
        return xddot

    def get_torque(self, time):
        make_state(time)
        Q = self.area*self.xdot
        Qdot = self.area*self.xddot
        P = masterpump.get_pressure(Q,Qdot)
        F_rod = (self.m_p*self.xddot + P*self.area)/(
                cos(self.alpha)+self.mu*sin(self.alpha)*(self.xdot/abs(self.xdot)))
        torque = -F_rod*self.r_c*sin(pi-self.alpha-self.theta)
        return torque
