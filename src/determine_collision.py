from numpy import random
from src.return_collision_frequencies import *
import operator
import math
from numpy import isnan
import csv
from src.initial_conditions import *
import pandas as pd
from collections import OrderedDict


interpolatedCsDict = returnIterCsDictionary()
collision_frequencies = create_cf_dicts()



"""
There is no mathematical equation for returning energy lost for a given inelastic collision so we define them manually
by adding them to a dictionary. We define said dictionary within the function return_e_loss_from_outcome below. The 
function takes two parameters, collision_type and return dictionary.

By default, the function will return a quantized energy to be lost by our incident electron based on the collision type
entered into the first parameter. But if for any reason it turns out to be more convenient to create the dictionary once
and call the dictionary rather than the function below when a value is needed, just enter false for the second parameter.

In retrospect, that does sound faster since every time we are calling this function we are creating the dictionary all
over again. 
"""

def prob_for_this_state(energy_eV, cs_process):
    """Determines the probability of a process (i.e. elastic, ionization, gp1, etc.) based on the electrons current
    energy (in electron volts.

    Args:
        energy_eV: The energy of the incident electron in electron volts (up to 100.00 eV and rounded to two decimal
        places, i.e. 1.25 eV, 29.42 eV, 100.00 eV).
        cs_process: The type of collision experienced by the incident electron.

    Returns:
        prob: A floating point number representing probability of the 'cs_process' occurring.
    """
    # eV is our energy in eV, but we run it through the ir() function to filter out NaN values and repeating decimals
    # from python representation errors. For more info see https://docs.python.org/2/tutorial/floatingpoint.html
    eV = ir(energy_eV)
    # We will be iterating through this function many times to find many CF's for various different processes
    if cs_process == "null":
        current_process_cf = NCF4Energy[ir(energy_eV)]
    else:
        current_process_cf = collision_frequencies[cs_process][eV]
    prob = current_process_cf / maxx
    return prob

def return_collision_type(energy_in_eV):

    # e is our energy in eV, but we run it through the ir() function to filter out NaN values and repeating decimals
    # from python representation errors. For more info see https://docs.python.org/2/tutorial/floatingpoint.html
    e = ir(energy_in_eV)
    # r is our randomly generated float between 0 and 1
    r = random.rand()
    # process_list will be a list of lists. Each list will consist of
    # ([process, start_probability, end_probability])
    process_list = []

    current = 0 # End of the probability bar. If P(elastic)=25, after adding elastic to the process list current = 25.
    for process in energy_loss_dictionary:
        if energy_in_eV >= energy_loss_dictionary[process]:
            prob = prob_for_this_state(e, process)
            start_prob = current
            end_prob = start_prob + prob
            process_list.append([process, start_prob, end_prob])
            current += prob

    # First let's add the probability for no collision, AKA 'null'
    # Just using elastic as a
    prob = NCF4Energy[e]/maxx
    start_prob = current
    end_prob = start_prob + prob
    process_list.append(["null", start_prob, end_prob])
    current += prob
    for t in range(len(process_list)):
        current_process = process_list[t]
        if float_in_range(r, current_process[1], current_process[2]):
            return current_process[0]
    return "null"

def sum_excitation_collision_frequencies_abc(all_collision_frequencies):
    collision_frequencies = dict(all_collision_frequencies)
    SumExciteCS = dict()
    del collision_frequencies["elastic"]
    del collision_frequencies["ionization"]
    for currentEnergy in interpolated_ev:
        currentSum = 0
        for specific_excitation in collision_frequencies:
            currentSum += collision_frequencies[specific_excitation][currentEnergy]
        SumExciteCS[currentEnergy] = currentSum
    return SumExciteCS



def condense_collision_frequency_dictionary_abc(all_collision_frequencies):
    """Because we have twenty or so different types of excitation, I wrote this function to sum each individual
    excitation collision frequency into one value keyed into the dictionary as 'sumExcite'.

    Args:
        AllCS: All energy values experienced at any step.

    Returns:
        condensed_collision_frequencies: Dictionary where all excitation collision frequencies are summed into one
        total value per discrete energy (in electron volts).

    """
    condensed_collision_frequencies = dict()
    condensed_collision_frequencies["elastic"] = all_collision_frequencies["elastic"]
    condensed_collision_frequencies["ionization"] = all_collision_frequencies["ionization"]
    condensed_collision_frequencies["sumExcite"] = sum_excitation_collision_frequencies_abc(all_collision_frequencies)
    return condensed_collision_frequencies

