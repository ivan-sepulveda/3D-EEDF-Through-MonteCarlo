# from run_simulation import run_trials
import numpy as np
# import initial_conditions

example_energies = np.asarray([0.5, 1.5, 2.5, 3.5])
example_eedf_values = np.asarray([0,112, 34, 41, 89 ]) # EEDF = Electron Energy Distribution Function
example_normalized_eedf_values = np.divide(example_eedf_values, sum(example_eedf_values))

def create_electron_collision_frequencies(energy_values, cross_section_values):
	pass


create_electron_cross_sections(0, 0)