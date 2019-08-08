import random # To generate random values for Monte Carlo method
import math
import numpy as np
from determine_collision import *
from return_collision_frequencies import *
from initial_conditions import *
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

status_format = "{0}x = {1} | {0}y = {2} | {0}z = {3}"

iterEVs = return_iter_evs()
MaxTotalCF = maxx

def return_beta(Vx, Vy):
    if Vx == 0:
        return math.pi/2
    return atan(x=(Vy/Vx))

def return_alpha(Vx, Vy, Vz):
    if Vx == 0 and Vy == 0:
        if Vz > 0:
            return 0
        if Vz < 0:
            return math.pi
    if Vz == 0:
        return math.pi/2
    return atan((np.sqrt(Vx**2 + Vy**2)/Vz))

def return_dt(energy_inJoules):
    if abs(energy_inJoules*Joules2EV) < 0.01:
        return 100 * (10 ** (-9))
    return -np.log(1-random.rand())/total_CF(energy_inJoules*Joules2EV, interpolatedCsDict)

class Electron:
    def __init__(self):
        self.current_time, self.current_energy_J = 0, 0
        self.x_pos, self.y_pos, self.z_pos = 0, 0, 0
        self.current_vx, self.current_vy, self.current_vz = 0, 0, 0
        self.acceleration_x, self.acceleration_y, self.acceleration_z = 0, 0, 0
        self.list_x, self.list_y, self.list_z = [0], [0], [0]
        self.list_t, self.list_e, self.list_v = [0], [0], [0]
        self.list_vx, self.list_vy, self.list_vz = [0], [0], [0]
        self.collisions = dict()
        self.moves_left = 0
        for key in interpolatedCsDict: self.collisions[key] = []
    def update_data_lists(self, x, y, z, t, e_eV, vx, vy, vz, v):
        self.list_x.append(x)
        self.list_y.append(y)
        self.list_z.append(z)
        self.list_t.append(t)
        if (not isnan(e_eV)) and (e_eV <= energy_ev_cutoff):
            self.list_e.append(e_eV)
        self.list_vx.append(vx)
        self.list_vy.append(vy)
        self.list_vz.append(vz)
        self.list_v.append(v)
    def okay_to_propagate(self):
        if self.moves_left >= 1 and self.current_energy_J*Joules2EV <= energy_ev_cutoff:
            return True
        else:
            return False
    def print_electron_status(self, position = True, velocity = True, energy = True, acceleration = True):
        print("Electron {0} at {1} ns".format(id(self), round(self.current_time*(10**9), 2)))
        if position:
            print(status_format.format("", self.x_pos, self.y_pos, self.z_pos))
        if velocity:
            print(status_format.format("V", self.current_vx, self.current_vy, self.current_vz))
            print("V = {0}".format(np.sqrt(self.current_vx**2 + self.current_vy**2 + self.current_vz**2)))
        if energy:
            print("Energy (Joules) = {0}".format(self.current_energy_J))
            print("Energy (ev) = {0}".format(self.current_energy_J*Joules2EV))
        if acceleration:
            print(status_format.format("A", self.acceleration_x, self.acceleration_y, self.acceleration_z))
            print("A = {0}".format(np.sqrt(self.acceleration_x**2 + self.acceleration_y**2 + self.acceleration_z**2)))
        print()
    def check_for_collision(self):
        V_i = np.sqrt(abs(self.current_vx**2 + self.current_vy**2 + self.current_vz**2))
        alpha = return_alpha(self.current_vx, self.current_vy, self.current_vz)
        beta = return_beta(self.current_vx, self.current_vy)
        E_ev = self.current_energy_J*Joules2EV
        process = return_collision_type(E_ev)
        if process != "null":
            theta = math.pi * (random.random() - 0.5) # Random scattering angle 1/2
            phi = math.pi * (random.random() - 0.5) # Random scattering angle 2/2
            # Now we can create the new Vx, Vy, and Vz ratios of the final velocity
            x_ratio = sin(alpha)*cos(beta)*cos(theta) \
                      - sin(beta)*sin(phi)*sin(theta)\
                      + sin(theta)*cos(alpha)*cos(beta)*cos(phi)

            y_ratio = sin(alpha)*sin(beta)*cos(theta)\
                      + sin(beta)*sin(theta)*cos(alpha)*cos(phi)\
                      + sin(theta)*cos(beta)*sin(phi)

            z_ratio = cos(alpha)*cos(theta)\
                      - sin(alpha)*sin(theta)*cos(phi)

            DeltaE_J = energy_loss_dictionary[process]*EV2J
            V_f = (m_e*V_i*cos(theta))/(argonMass*np.sqrt(1 + (m_e/argonMass)))
            V_f += np.sqrt(V_i**2 - ((V_i**2)*m_e/argonMass) - (2*DeltaE_J/m_e))
            V_f = V_f/np.sqrt(1 + (m_e/argonMass))
            self.current_vx, self.current_vy, self.current_vz = V_f*x_ratio, V_f*y_ratio, V_f*z_ratio
            self.current_energy_J = 0.5*m_e*(V_f**2) # Updating Energy
            self.update_data_lists(self.x_pos, self.y_pos, self.z_pos, # Updating data lists
                                   self.current_time, self.current_energy_J*Joules2EV, # Updating data lists
                                   self.current_vx, self.current_vy, self.current_vz, V_f) # Updating data lists
            self.moves_left -= 1 # Reduce allowed propagations by 1

    def propagate_electron(self):
        self.moves_left -= 1 # Reduce remaining propagations by 1
        dt = return_dt(self.current_energy_J) # Determining the time interval of the propagation
        self.x_pos += (self.current_vx * dt) + (0.5 * self.acceleration_x * (dt ** 2)) # Propagating Electron
        self.y_pos += (self.current_vy * dt) + (0.5 * self.acceleration_y * (dt ** 2)) # Propagating Electron
        self.z_pos += (self.current_vz * dt) + (0.5 * self.acceleration_z * (dt ** 2)) # Propagating Electron
        self.current_vx = self.current_vx + (self.acceleration_x * dt) # Updating velocity (acceleration from E-Field)
        self.current_vy = self.current_vy + (self.acceleration_y * dt) # Updating velocity (acceleration from E-Field)
        self.current_vz = self.current_vz + (self.acceleration_z * dt) # Updating velocity (acceleration from E-Field)
        Vf = np.sqrt(self.current_vx ** 2 + self.current_vy ** 2 + self.current_vz ** 2) # Speed post propagation
        self.current_energy_J = 0.5 * m_e * (Vf**2) # Updating Energy
        self.current_time += dt # Updating time
        self.update_data_lists(self.x_pos, self.y_pos, self.z_pos,
                               self.current_time, self.current_energy_J*Joules2EV,
                               self.current_vx, self.current_vy, self.current_vz, Vf)
        if self.moves_left >= 1 and self.current_energy_J*Joules2EV <= energy_ev_cutoff:
            self.check_for_collision()
        return self

def send_electron(print_status=False):
    """Sends the current electron instance through the simulated plasma.

    Args:
    	  print_status: If True, prints the coordinates, velocity, accelerations, and energy of the electron as it moves.

    Returns:
        electron: The same electron instance after it has either completed all its propagations or hit the energy limit
    """
    electron = Electron()
    electron.moves_left += motions_per_electron
    electron.acceleration_x = ((electronCharge * E_x) / m_e)
    electron.acceleration_y = ((electronCharge * E_y) / m_e)
    electron.acceleration_z = ((electronCharge * E_z) / m_e)

    while electron.okay_to_propagate():
      if print_status:
        electron.print_electron_status() # Uncomment to print status as electron propagates space

      electron = electron.propagate_electron()

    return electron


def run_electrons():
    """When completed, this function will return a list of electron class instances. Each instance has accessible lists
    of their respective positions, velocities, accelerations, and energies experienced over a corresponding time.
    The run_electrons function below takes 2 arguments.

    Returns:
        electrons: list of electorn class instances
    """
    electrons = [] # List of Electron class instances we will append to
    for i in range(number_of_electrons):
        electrons.append(send_electron())

    return electrons

