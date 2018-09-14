from numpy import random
from ReturnAllInterpolatedCSDatFiles import *
import operator
import math
from numpy import isnan
import csv
from initialConditions import *

interpolatedCsDict = returnIterCsDictionary()
cf_dictionary__ = create_cf_dicts()



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

    # eV is our energy in eV, but we run it through the ir() function to filter out NaN values and repeating decimals
    # from python representation errors. For more info see https://docs.python.org/2/tutorial/floatingpoint.html
    eV = ir(energy_eV)
    # We will be iterating through this function many times to find many CF's for various different processes
    if cs_process == "null":
        current_process_cf = NCF4Energy[ir(energy_eV)]
    else:
        current_process_cf = cf_dictionary__[cs_process][eV]
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

    possible_processes = []
    eLoss_dict = e_loss
    for pro in eLoss_dict:
        if energy_in_eV >= eLoss_dict[pro]:
            possible_processes.append(pro)
    """
    'current' is the start of the probability bar. It starts at 0
    If the probability for elastic equals 25, after adding elastic to the process list
    'current' will equal 20. If the probability for ionization equals 10, 'current' will then
    equal 30 and so on.
    """
    current = 0
    for process in possible_processes:
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

def return_top_cfs(interpolatedCsDict_, howManyTopCFs = 1):
    # First we make a copy of the original Cross Section Dictionary
    copyIterCS = dict(interpolatedCsDict_)
    # Next we make this list below where our results will be returned
    TopCFs = []
    """We use a for loop to iterate through a copy of interpolatedCsDict to find the process with the highest CF.
    After it finds the top CF it deletes that process from the dictionary to find 2nd place. We then go through the
    new shorter version of aforementioned copy to find third place, and so forth."""
    # For loop in range of how ever many top Collision Frequencies we want to search for
    for k in range(howManyTopCFs):
        # Contenders will be a list of all the Max TCFs this round.
        contenders = list()
        # The current Max Collision Frequency is 0, and we will see it increase as we parse the dictionaries
        MaxCF = 0
        # For every process in our copy, we will compute the CF for every energy (hence the double for loop)
        for proces in copyIterCS:
            for ener in copyIterCS[proces]:
                # The collision energy for the current process parsed is calculated by calling the c_f function
                currentCollisionFrequency = cf_dictionary__[proces][ener]
                # If it us larger than the previous CF, it is the new maximum.
                if currentCollisionFrequency > MaxCF:
                    MaxCF = currentCollisionFrequency
                    # We also add it to the list of top CF's.
                    # The contenders list will be a list of lists. The length of the contenders list will be
                    # equal the the input parameter howManyTopCFs. Each list in contenders will describe everything
                    # we need to know about the top collision frequencies.
                    contenders.append([proces, ener, currentCollisionFrequency])
                # Otherwise we move on to check the next energy value.
        TopCFs.append([proces, ener, currentCollisionFrequency])
        # How we delete the previous Top CF Process from our dictionary. We cannot de
        del copyIterCS[contenders[-1][0]]
    return TopCFs