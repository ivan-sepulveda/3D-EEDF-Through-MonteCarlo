# Importing Modules
import random
import math
from math import pi as pi
from numpy import *
import numpy as np
from returnOutcome import *
from ReturnAllInterpolatedCSDatFiles import *
from initialConditions import *
import os

def return_iter_evs():
    os.chdir(os.getcwd()+"/CrossSectionsTheoretical")
    i_ev = return_interpolated_ev("xsigma.dat")
    os.chdir(os.getcwd().replace("/CrossSectionsTheoretical", ""))
    return i_ev

sin = lambda x: math.sin(x)
cos = lambda x: math.cos(x)
tan = lambda x: math.tan(x)
atan = lambda x: math.atan(x)

iterEVs = return_iter_evs()
MaxTotalCF = maxx

"""
Abbreviations in this code:
ncf = null collision frequency
tcf_1s = total collision frequency for 1 state
tcf_all = total collision frequency for total of all states
"""

"""
This function was previously called "iround" but I shortened it down to "ir". Python at times has issues processing 
certain floats so this is a makeshift solution. First you use the built-in round() function to round the float down to 
2 decimals. Second we convert that float into a string immediately, forcing python to format it correctly. Lastly we 
turn said string back into a float.

For further elaboration on aforementioned float issues, see below.

a = 0.01
for i in range(100):
    print(a)
    a += 0.01

Normally the function above will give
--------------------------------------------------
0.01
0.02
0.03
0.04
0.05
0.060000000000000005
0.07
0.08
0.09
0.09999999999999999
0.10999999999999999
0.11999999999999998
0.12999999999999998
0.13999999999999999
0.15
0.16

And so forth.

"""

def ir(ener):
    return float(str(round(ener, 2)))

# Previous Version has this angle as theta, I'm just aligning with Kushner notes
def return_beta(Vx, Vy):
    if Vx == 0:
        return pi/2
    else:
        return atan(x=(Vy/Vx))

def return_alpha(Vx, Vy, Vz):
    if Vx == 0 and Vy == 0:
        if Vz > 0:
            return 0
        if Vz < 0:
            return pi
    if Vz == 0:
        return pi/2
    else:
        return atan(sqrt(Vx**2 + Vy**2)/Vz)

def return_dt(energy_inJoules):
    if abs(energy_inJoules*Joules2EV) < 0.01:
        return 100 * (10 ** (-9))
    else:
        return -np.log(1-random.rand())/total_CF(energy_inJoules*Joules2EV, interpolatedCsDict)

class create_electron():
    def __init__(self):
        self.current_time = 0
        self.x_pos = 0
        self.y_pos = 0
        self.z_pos = 0
        self.curr_V_x = 0
        self.curr_V_y = 0
        self.curr_V_z = 0
        self.current_accel_x = 0
        self.current_accel_y = 0
        self.current_accel_z = 0
        self.current_energy_J = 0
        self.list_t = [0]
        self.list_x = [0]
        self.list_y = [0]
        self.list_z = [0]
        self.list_e = [0]
        self.list_vx = [0]
        self.list_vy = [0]
        self.list_vz = [0]
        self.list_v = [0]
        self.collisions = dict()
        self.moves_left = 0
        for key in interpolatedCsDict: self.collisions[key] = []
    def update_data_lists(self, x, y, z, t, e_eV, vx, vy, vz, v):
        self.list_x.append(x)
        self.list_y.append(y)
        self.list_z.append(z)
        self.list_t.append(t)
        if not isnan(e_eV):
            self.list_e.append(e_eV)
        self.list_vx.append(vx)
        self.list_vy.append(vy)
        self.list_vz.append(vz)
        self.list_v.append(v)
    def printstats(self, xyzt = True, vel = True, energy = True, accel = True):
        # Below are all commands the electron's stats
        if xyzt:
            print("x = {0}".format(self.x_pos))
            print("y = {0}".format(self.y_pos))
            print("z = {0}".format(self.z_pos))
            print("t = {0} ns".format(self.current_time*(10**9)))
        if vel:
            print("Vx = {0}".format(self.curr_V_x))
            print("Vy = {0}".format(self.curr_V_y))
            print("Vz = {0}".format(self.curr_V_z))
            print("V = {0}".format(sqrt(self.curr_V_x**2 + self.curr_V_y**2 + self.curr_V_z**2)))
        if energy:
            print("Energy (ev) = {0}".format(self.current_energy_J*Joules2EV))
        if accel:
            print("Ax = {0}".format(self.current_accel_x))
            print("Ay = {0}".format(self.current_accel_y))
            print("Az = {0}".format(self.current_accel_z))
            print("A = {0}".format(sqrt(self.current_accel_x**2 + self.current_accel_y**2 + self.current_accel_z**2)))
    def check_for_collision(self):
        V_i = sqrt(self.curr_V_x**2 + self.curr_V_y**2 + self.curr_V_z**2)
        beta = return_beta(self.curr_V_x, self.curr_V_y)
        alpha = return_alpha(self.curr_V_x, self.curr_V_y, self.curr_V_z)
        # To see updated incident angles as they occur, uncomment lines bellow
        # print("Alpha: {0}".format(rad2deg(alpha)))
        # print("Beta: {0}".format(rad2deg(beta)))
        E_ev = self.current_energy_J*Joules2EV
        process = return_collision_type(E_ev)
        if process != "null":
            # First let's generate the 2 random scattering angles
            theta = pi * (random.random() - 0.5)
            phi = pi * (random.random() - 0.5)
            # Now we can create the new Vx, Vy, and Vz ratios of the final velocity
            x_ratio = cos(beta)*cos(alpha)*sin(theta)*cos(phi) + cos(beta)*sin(alpha)*cos(theta) - sin(beta)*sin(theta)*sin(phi)
            y_ratio = sin(beta)*cos(alpha)*sin(theta)*cos(phi) + sin(beta)*sin(alpha)*cos(theta) - cos(beta)*sin(theta)*sin(phi)
            z_ratio = cos(alpha)*cos(theta) - sin(alpha)*sin(theta)*cos(phi)
            # To see updated collision angles as they occur, uncomment lines bellow
            # print("Theta: {0} (deg)".format(rad2deg(theta)))
            # print("Phi: {0} (deg)".format(rad2deg(theta)))
            DeltaE_J = e_loss[process]*EV2J
            V_f = sqrt(abs(V_i)**2 - (2*DeltaE_J/m_e))
            self.curr_V_x = V_f*x_ratio
            self.curr_V_y = V_f*y_ratio
            self.curr_V_z = V_f*z_ratio
            # Updating Energy
            self.current_energy_J = 0.5*m_e*(V_f**2)
            # Updating data lists
            self.update_data_lists(self.x_pos, self.y_pos, self.z_pos,  self.current_time, self.current_energy_J*Joules2EV, self.curr_V_x, self.curr_V_y, self.curr_V_z, V_f)
            # Reduce allowed propagations by 1
            self.moves_left -= 1
            return self
        if process == "null":
            # Uncomment line below if you would like see when 'null' collisions happen
            # print("No collision. Moving on to next propagation")
            pass
        return self

    def propagate_electron(self):
        # First we must decide the time interval of the propagation
        dt = return_dt(self.current_energy_J)
        # Propagating Electron
        self.x_pos += (self.curr_V_x * dt) + (0.5 * self.current_accel_x * (dt ** 2))
        self.y_pos += (self.curr_V_y * dt) + (0.5 * self.current_accel_y * (dt ** 2))
        self.z_pos += (self.curr_V_z * dt) + (0.5 * self.current_accel_z * (dt ** 2))
        # Updating velocity components because electron E_field accelerates electron
        self.curr_V_x = self.curr_V_x + (self.current_accel_x * dt)
        self.curr_V_y = self.curr_V_y + (self.current_accel_y * dt)
        self.curr_V_z = self.curr_V_z + (self.current_accel_z * dt)
        # Vf will be the velocity post propagation
        Vf = sqrt(self.curr_V_x ** 2 + self.curr_V_y ** 2 + self.curr_V_z ** 2)
        # Updating Enegy
        self.current_energy_J = 0.5 * m_e * (Vf ** 2)
        # Updating time
        self.current_time += dt
        self.update_data_lists(self.x_pos, self.y_pos, self.z_pos, self.current_time,
                               self.current_energy_J * Joules2EV, self.curr_V_x, self.curr_V_y, self.curr_V_z, Vf)
        # Reduce allowed propagations by 1
        self.moves_left -= 1
        if self.current_energy_J*Joules2EV < 99.99 and self.moves_left >= 1:
            # Uncomment below if you would like to see every time the code checks for an electron collision post prop.
            # print("Now checking for collision")
            self.check_for_collision()
        return self

def run_simuation3(number_of_motions=MPE):
    eci = create_electron()
    eci.moves_left += number_of_motions
    eci.current_accel_x = ((electronCharge * E_x) / m_e)
    eci.current_accel_y = ((electronCharge * E_y) / m_e)
    eci.current_accel_z = ((electronCharge * E_z) / m_e)

    while eci.moves_left > 0:
        if eci.current_energy_J*Joules2EV <= 100.00:
            # To print electron status as it propagates, uncomment line below.
            # print("\n\nPropagation: {0}".format(number_of_motions - eci.moves_left + 1))
            eci = eci.propagate_electron()
        else:
            # Below we are deleting the outlier point
            eci.list_e = eci.list_e[:-3]
            eci.list_x = eci.list_x[:-3]
            eci.list_y = eci.list_y[:-3]
            eci.list_z = eci.list_z[:-3]
            eci.list_t = eci.list_t[:-3]
            eci.list_v = eci.list_v[:-3]
            eci.list_vx = eci.list_vx[:-3]
            eci.list_vy = eci.list_vy[:-3]
            eci.list_vz = eci.list_vz[:-3]
            print("Ceasing propagation: electron has over 100 eV")
            break
    return eci