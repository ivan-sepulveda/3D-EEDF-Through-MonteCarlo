from src.propagate_electrons import send_electron, run_electrons
from src.return_collision_frequencies import silent_create, silent_remove, return_to_start_directory
import src.plot_simulation_data
from src.initial_conditions import *
import src.propagate_electrons
import os

def run_trials():
    """ In Layman's terms, run_electrons creates all the data while run_trials exports the results into excel files and
     graphs with and creates a plot that compares our results to those from similar studies (cited in the report).

    Returns:
         list_electrons: list of electron instances
         all_energies_eV: list of all energies experienced by all electrons at any given point.
    """
    return_to_start_directory()
    all_results_directory = os.getcwd() + "/Results" # Directory of all the results.
    silent_create(all_results_directory)
    os.chdir(all_results_directory) # Now we change our current directory to the results_directory
    simulation_conditions_string = "/Ex={0} Ey={1} Ez={2} NOE={3} MPE={4}"\
        .format(E_x, E_y, E_z, NOE, MPE)
    current_simulation_results_directory = all_results_directory + simulation_conditions_string
    silent_remove(current_simulation_results_directory, is_directory=True)
    os.makedirs(current_simulation_results_directory) # Creating results sub-directory specific to current conditions
    os.chdir(current_simulation_results_directory) # Changing
    list_electrons = run_electrons() # Master list of all electrons
    all_energies_eV = src.plot_simulation_data.consolidate_energy_lists(list_electrons) # Energies electrons experienced
    src.plot_simulation_data.write_eedf_spreadsheet(all_energies_eV) # Electron Energy Distribution Function Spreadsheet
    src.plot_simulation_data.write_normalized_eedf_spreadsheet(all_energies_eV) # Normalized EEDF Spreadsheet
    src.plot_simulation_data.plot_energy_distribution_function(all_energies_eV) # PDF of EEDF
    # plot.plot_paths(list_electrons) # Plots the electrons paths in 3D (computationally expensive)
    return_to_start_directory() # Returning to original directory
    return all_energies_eV, list_electrons

print("staring to run simulation")
energies, electrons = run_trials()
print("done running simulation")

