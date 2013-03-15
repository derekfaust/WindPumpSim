#!/usr/bin/env python2

import numpy as np
import math

class Turbine:
    def __init__(self, output):
        self.omega = output[:,1]
        self.torque = output[:,2]
    def torque(self, omega):
        return np.interp(omega, self.omega, self.torque, left=0, right=0)
    def omega(self, torque):
        return np.interp(torque, self.torque, self.omega, left=0, right=0)

class Pump:
    g = 9.8
    TubeRadius=.01
    TubeLength=1
    rho_w=1000
    vol_pumped=0
    max_pressure=rho_w*g*tube_length
    tube_area = self.TubeRadius**2*math.pi
    def __init__(self):
    def back_pressure(self):
        return min(self.max_pressure, self.rho_w*self.g*(self.vol_pumped/tube_area)

class SingleActionPistonPump(Pump):
    arm_length=.12
    piston_radius=.05
    crank_radius=.04
    piston_area=piston_radius**2*math.pi
    initial_angle=0                         #Zero Corresponds to Fully Extended
    angle=initial_angle
    friction_force=5
    zero_pos=crankradius+arm_length
    def _init__(self):
    def piston_position(self):
        crank_comp = self.crank_radius*math.cos(self.angle)
        arm_comp = (self.arm_length**2 - (self.crank_radius*math.sin(self.angle))**2)**.5
        return zero_pos-arm_comp-crank_comp
    def piston_force(self):
    def required_torque(self):

        
