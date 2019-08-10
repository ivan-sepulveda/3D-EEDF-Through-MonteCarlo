import numpy as np
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from bolos import parser, grid, solver
# from initial_conditions import *
#
# gr = grid.LinearGrid(0, 100., 200)
#
# fx = np.linspace(0, max_energy_ev, max_energy_ev/0.5 + 1)
# # print(fx)
#
# binzzzz = [(fx[i+1] + fx[i])/2 for i in range(len(fx) - 1)]
#
# # print(binzzzz)
# sqrt_binz = np.sqrt(binzzzz)
# # print(len(fx))
# # print(sqrt_binz)
# # print(len(binzzzz))
#
#
# boltzmann = solver.BoltzmannSolver(gr)
#
#
#
# with open("LXCat-June2013.txt") as fp:
# 	processes = parser.parse(fp)
# 	boltzmann.load_collisions(processes)
#
#
# boltzmann.kT = 300 * solver.KB / solver.ELECTRONVOLT
# boltzmann.EN = 25 * solver.TOWNSEND
# boltzmann.init()
#
# boltzmann.set_density(species='Ar', density=1.0)
#
# fMaxwell = boltzmann.maxwell(2.0)
#
# # print("Mean Energy:", boltzmann.mean_energy(fMaxwell))
#
#
# finterp = boltzmann.grid.interpolate(fMaxwell, gr)
# boltzmann.init()
# f = boltzmann.converge(finterp, maxn=100, rtol=1e-5)
# boltzmann.init()
# # print("Mean Energy:", boltzmann.mean_energy(f))
# # print(f)
# # print(len(f))
#
# plt.plot([0] + binzzzz, [0] + np.multiply(f, sqrt_binz).tolist())
# plt.xlim(0, eedf_x_lim)
# plt.show()