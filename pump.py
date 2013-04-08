#Import required modules
from __future__ import division
import numpy as np
from math import *

class Turbine:
    #Turbine object for pump simulation
    def __init__(self, I, performance):
        #Initialize Turbine object with provided specifications
        self.I = I
        self.omega = performance[:,0]
        self.torque = performance[:,1]
    
    def Tin(self, omega):
        #Determine the torque provided by the wind
        Tin = np.interp(omega, self.omega, self.torque)
        return Tin
    
    def power(self, state):
        #Determine the power
        omega = state[1]
        return omega*self.Tin(omega)
    
    def de_omegadot_coef(self, state):
        #Determine angular acceleration coefficient for dE/dt expression
        omega = state[1]
        return self.I*omega

class Pump:
    #Pump object for pump simulation
    
    #Universal pump parameters
    g = 9.8                 #Gravity acceleration (m/s^2)
    tube_radius=.01         #Radius of the tube (m)
    tube_height=1.5         #Height of water column to overcome (m)
    tube_length=2           #Tube length (m)
    rho_w=1000              #Density of water [kg/m^3]
    friction_factor = .019  #Darcy friction factor
    tube_area = tube_radius**2*pi       #Cross-sectional area of tube
    max_pressure = rho_w*g*tube_height  #Maximum pressure due to gravity

    def __init__(self, drivesystem, mechanisms):
        #Initialize pump object with known drive and mechanisms
        self.drive = drivesystem        #Drive system, such as turbine
        self.mechanisms = mechanisms    #Mechanism for converting rotation to flow
    
    def backpressure(self, state, Q):
        #Calculate the backpressure that must be provided by the mechanism
        vol_pumped = state[2]   #Get volume pumped from state
        #Backpressure is the maximum, or due to column height
        head_loss = self.friction_factor*(self.tube_length/(2*self.tube_radius))*(self.rho_w*(Q/self.tube_area)**2/2)
        backpressure = min(self.rho_w*self.g*(vol_pumped/self.tube_area), self.max_pressure)
        pressure = head_loss+backpressure
        return pressure

    def net_power(self, state):
        #Calculate the net power into the system
        mech_power = 0
        for mechanism in self.mechanisms:
            #Find the backpressure
            bp = self.backpressure(state, mechanism.Q(state))
            #Find the power leaving the system at the mechanisms
            mech_power += mechanism.power(state, bp,
                                               self.drive.Tin(state[1]))
        #Net power is power entering the system minus power leaving        
        pnet = self.drive.power(state) - mech_power      
        return pnet

    def omegadot(self, state):
        #Calculate the change in angular velocity (which drives the state
        #change of the entire system)
        #Determine the net power
        pnet = self.net_power(state)
        #Determine the omegadot coefficients in the
        #change in kinetic energy expression
        c_omegadot_t = self.drive.de_omegadot_coef(state)
        c_omegadot_m = 0
        for mechanism in self.mechanisms:
            c_omegadot_m += mechanism.de_omegadot_coef(state)
        #Determine the component of the change in kinetic
        #energy of the mechanism independent of omegadot
        b_de_m = 0        
        for mechanism in self.mechanisms:
            b_de_m += mechanism.de_offset(state)
        #Calculate omegadot given power, offsets, and coefficients
        #(See equation write-up for details on theory)
        omegadot = (pnet - b_de_m)/(c_omegadot_t+c_omegadot_m)
        #print pnet - b_de_m, state[1]
        return omegadot

    def statedot(self, t, state, p):
        #Determine the change in state of the system (differential equation)
        #The state variable is in the following form:
        #state = [theta, omega, vol]
        thetadot = state[1]                 #Change in theta is omega
        omegadot = self.omegadot(state)     #Change in omega is omegadot
        voldot = 0
        for mechanism in self.mechanisms:
            voldot += mechanism.Q(state)    #Change in volume pumped is flow
        #Assemble statedot variable for numerical integrator
        statedot = [thetadot, omegadot, voldot]
        return statedot

class Piston:
    #Piston mechanism for pump simulations
    
    #General piston parameters
    mu = .5
    
    def __init__(self, r_crank, r_rod, r_piston, m_piston, theta):
        #Initialize piston mechanism with 
        self.r_c = r_crank
        self.r_r = r_rod
        self.mass = m_piston
        self.area = r_piston**2*pi
        self.theta0 = theta
    
    def alpha(self, theta):
        #Determine the angle of the pushrod with respect to the piston access
        return asin(self.r_c/self.r_r*sin(theta))

    def xpos(self, state):
        #Calculate the position of the piston
        theta = state[0]+self.theta0    #Include initial theta offset
        comp1 = self.r_c*cos(theta)     #First term   
        comp2 = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2) #Second Term
        xpos = -comp1+comp2             #Sum the terms
        return xpos

    def xdot(self, state):
        #Calculate the velocity of the piston
        theta = state[0]+self.theta0    #Include the initial theta offset
        omega = state[1]                #Pull omega from state
        xdot = -omega*self.xdotcoef(theta)  #Calculate velocity
        return xdot
    
    def xdotcoef(self, theta):
        #Calculate the coefficient for the velocity calculation.
        #Used to split up and simplify the calculation
        comp1 = self.r_c*sin(theta)
        denom = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        comp2 = (self.r_c**2*sin(theta)*cos(theta))/denom
        return comp1-comp2

    def xddotcoef(self, theta):
        #Calculate the coefficient for the acceleration calculation.
        #Used to split up and simplify the calculation
        num = (self.r_c**2*cos(theta)*sin(theta))**2
        denom = self.r_r*sqrt(1-(self.r_c/self.r_r*sin(theta))**2)
        coef = (self.r_c**2/denom+num/(denom**3)+self.r_c*cos(theta))
        return coef

    def Q(self, state):
        #Calculate the flow rate of the piston
        xd = self.xdot(state)   #Get velocity
        if xd > 0:
            #If piston velocity is positive, flow rate is velocity*area
            q = xd*self.area
        else:
            #If velocity of piston is negative or zero, there is no flow
            q = 0
        return q
    
    def power(self, state, backpressure, torque):
        #Calculate the power leaving the system at the piston
        theta = state[0]+self.theta0    #Include initial theta offset
        #Power used to pump water is pressure times area times velocity
        pumpingpower = backpressure*self.area*self.xdot(state)
        #Calcuate tension in the rod for friction calculation
        T_rod = torque/(self.r_c*sin(pi-theta-self.alpha(theta)))
        #Calculate power lost to friction
        power_lost_friction = abs(self.mu*T_rod*sin(theta)*self.xdot(state))
        #Power leaving the system is the sum of lost powers
        power_out = pumpingpower + power_lost_friction
        return power_out

    def de_omegadot_coef(self, state):
        #Determine angular acceleration coefficient for dE/dt expression
        theta = state[0]+self.theta0    #Include initial theta offset
        omega = state[1]                #Pull anglar velocity
        #Calculate the coefficient
        coef = self.mass*omega*self.xdotcoef(theta)**2
        return coef

    def de_offset(self, state):
        #Determine angular acceleration offset for dE/dt expression
        theta = state[0]+self.theta0    #Include initial theta offset
        omega = state[1]                #Pull angular velocity
        #Calculate the coefficient
        offset = self.mass*omega**3*self.xdotcoef(theta)*self.xddotcoef(theta)
        return offset
