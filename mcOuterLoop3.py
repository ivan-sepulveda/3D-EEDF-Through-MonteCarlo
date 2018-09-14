from CreatePropagateCollide import run_simuation3
import plot as sim_plot
import shutil
import time
from initialConditions import *
import os
from ReturnAllInterpolatedCSDatFiles import goHomeDirectory
# Uncomment line below if you want to time how long it takes the code to run
starttime = time.time()




"""
The run_electrons function below takes 5 arguments.
1) number_of_electrons: number of electrons that will be created and propagated through the simulate Ar Plasma
2) motions_per_elec: number of occurrences the electron will be allowed to propagate within the plasma for an allotted
time 'dt' (also known as time of flight)
3/4/5) Strength of the electric field you are applying in the x, y, & z directions in unites of Volts per meter

For simplicity, I've predefined the function to run 100 electrons, through 100 propagations each in a electric field
with a magnitude of 50 Volts/meter purely in the -x direction.

When completed, this function will return a list of elctron class instances (ECI's). Each ECI has accessible lists of
their respective positions, velocities, accelerations, and energies experienced over a corresponding time.
"""

def run_electrons(number_of_electrons=NOE, motions_per_elec=MPE):
    # All the electrons are python class instances. The list 'electrons' created in the line below stores all these
    # class instances.
    electrons = []
    # Below we are using a 'for loop' to run the desired number of electrons through a set number of propagations
    for i in range(number_of_electrons):
        # Uncomment line below to updated status of which electron is being propagated.
        # print("Electron {0}".format(i+1))
        # The run_simulation3 function used below returns an electron class instance (ECI) of an electron already
        # propagated through the number of motions desired in as inputted into the run_electrons function. We use the
        # built-in append function to add this electron to our list.
        electrons.append(run_simuation3(motions_per_elec))
    return electrons


""" The function run_trials defined below pulls the same 5 parameters requested in run_electrons because the run_electrons
function is used for the output. In Layman's terms, run_electrons creates all the data while run_trials exports the
results into excel files and graphs with in addition to creating a plot that compares our results to those derrived
from similar studies listed below

The function returns the list of ECI's and a list of all energies experienced by all electron at any given point. 

Bolsig Reference
Semi-Emperical Study Reference """
def run_trials(number_of_electrons=NOE, steps_per_electron=MPE, Efieldx=E_x, Efieldy=E_y, Efieldz=E_z):
    goHomeDirectory()
    # This is the directory of all the results.
    results_directory = os.getcwd() + "/Results"
    # If it does not currently exist, the next two lines will create it using an if-statement to check if it exists.
    if not os.path.exists(results_directory):
            os.makedirs(results_directory)
    # Now we change our current directory to the results_directory
    os.chdir(results_directory)
    """ The SpecifyDir string below will be the name of the subdirectory within the results directory.
    We name the folder according to the Electric Field Conditions, Number of electrons propagated, and motions per
    Electron just to assist the user in finding their results based off the initial conditions they set. 
    
    Like before, if the path does not exist we create the path. If it does exists we delete the previous files and
    write the new results"""
    SpecifyDir = results_directory+"/Ey={3}Ex={0}Electrons={1}MotionsPer={2}".format(Efieldx, number_of_electrons, steps_per_electron, Efieldy)
    if os.path.exists(SpecifyDir):
        shutil.rmtree(SpecifyDir)
    os.makedirs(SpecifyDir)
    os.chdir(SpecifyDir)
    # Below will be the master list of all electrons
    list_electrons = run_electrons(number_of_electrons, steps_per_electron)
    # Using the function CombineLists_PullEnergy, we pull every energy experienced by every single electron at any point
    all_E = sim_plot.CombineLists_PullEnergy(list_electrons)
    # Now create the Electron Energy Distribution Function Spreadsheet
    sim_plot.NumElecVsEnergySpreadsheet(all_E)
    # If you would like to create the Electron Energy distribution spreadsheet, uncomment following two lines
    print("Should be creating the spreadsheet now")
    sim_plot.NumElecVsEnergySpreadsheet(all_E)
    # If you would like to plot the Electron Energy distribution graph, uncomment the following two lines
    print("Should be creating the plot now")
    print("all_E:" , all_E)
    sim_plot.plot_eV_vs_num_elec(all_E)
    # If you would like to plot the electrons paths, uncomment the following line
    sim_plot.plot_paths(list_electrons)
    # Now we return to our original directory
    goHomeDirectory()
    return all_E, list_electrons


print("Initiating: run_trials")
Energy, Electrons = run_trials(number_of_electrons=25, Efieldx=10)
# EnergyA, ElectronsA = run_trials(number_of_electrons=25, Efieldx=-25)
# EnergyB, ElectronsB = run_trials(number_of_electrons=250, Efieldx=-10)
# EnergyC, ElectronsC = run_trials(number_of_electrons=250, Efieldx=-25)
print("Complete")
# Uncomment below if you want to time the code.
print('It took', time.time()-starttime, 'seconds.')
SystemExit()
