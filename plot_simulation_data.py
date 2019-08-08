from matplotlib import pylab # Plotting
import matplotlib.pyplot as plt # Plotting
import matplotlib.patches as mpatches # Plotting
import matplotlib.legend as lgnd # Plotting
from matplotlib.font_manager import FontProperties # Plotting
import numpy as np # Math
from math import *  # Math
import scipy  # Math
from decimal import Decimal # Math
from scipy.interpolate import interp1d # Math
from scipy import integrate  # Math
from statistics import mean # Math
from return_collision_frequencies import * # Project Related Data/Functions
from determine_collision import * # Project Related Data/Functions
from propagate_electrons import run_electrons # Project Related Data/Functions
from os import * # Switching directories
import xlsxwriter # Write to Excel
import errno # Error handling
from collections import OrderedDict # Functionality
from textwrap import wrap # No idea
import pandas as pd # Working on it
import glob # No idea

def totalcfs_and_max2(allCS, excludeElasIon = False, returnIncidentEnergy=False):
    """Sums up the collision frequencies from the dictionaries in allCS (all Collision Frequencies) 
    into a dictionary of their sums.

    Args:
        allCS: Dictionary of Collision Frequencies.

    Returns:
        EnergyTotalCF: Dictionary where keys are energy values in electron volts that correspond to a total collision
        frequency for an electron at that energy value.
    """
    iterCS = dict(allCS)
    if excludeElasIon:
        del iterCS["elastic"]
        del iterCS["ionization"]
    EnergyTotalCF = dict() # Creating Dictionary
    allTot = [] # Creating empty list for purpose of finding list max using built in python max()
    interpolated_energy_list = [] # Creating corresponding energy list so that index of energy matches index of TCF
    now = 0.00
    while now < 100.01:
        t = total_CF(ir(now), iterCS)
        interpolated_energy_list.append(now)
        allTot.append(t)
        EnergyTotalCF[ir(now)] = t
        now += 0.01
    index, value = max(enumerate(allTot), key=operator.itemgetter(1))
    if returnIncidentEnergy:
        return r.ir(interpolated_energy_list[index])
    return EnergyTotalCF


def plot_paths(list_of_electrons):
    """Plots paths of all the electrons using Mayavi 3D image package

    Args:
      list_of_electrons: A list/array of electron class instances.
    """
    N = len(list_of_electrons[0].list_t) # number of points per line
    from mayavi import mlab
    mlab.figure(1, size=(400, 400), bgcolor=(0, 0, 0))
    mlab.clf()
    x_positions, y_positions, z_positions = list(), list(), list() # Create lists of coordinate positions
    timestamps, connections = list(), list() # Create lists of coordinate connections
    index = 0 # The index of the current point in the total amount of points
    for i in range(len(list_of_electrons)): # Create each line one after the other
        x_positions.append(np.asarray(list_of_electrons[i].list_x))
        y_positions.append(np.asarray(list_of_electrons[i].list_y))
        z_positions.append(np.asarray(list_of_electrons[i].list_z))
        timestamps.append(np.asarray(list_of_electrons[i].list_t))
        """This is the tricky part: in a line, each point is connected
        to the one following it. We have to express this with the indices
        of the final set of points once all lines have been combined
        together, this is why we need to keep track of the total number of
        points already created (index). """
        connections.append(np.vstack( # Collapse them in one array before plotting
            [np.arange(index, index + N - 1.5), np.arange(index + 1, index + N - .5)]).T)
        index += N

    # Now collapse all positions, scalars and connections in big arrays
    x_positions, y_positions, z_positions = np.hstack(x_positions), np.hstack(y_positions), np.hstack(z_positions)
    timestamps, connections = np.hstack(timestamps), np.vstack(connections)
    src = mlab.pipeline.scalar_scatter(x_positions, y_positions, z_positions, timestamps) # Create the points
    src.mlab_source.dataset.lines = connections # Connect them
    src.update() # Connect them
    lines= mlab.pipeline.tube(src, tube_radius=0.005, tube_sides=6) # The stripper filter cleans up connected lines
    mlab.pipeline.surface(lines, colormap='Accent', line_width=1, opacity=.4) # Finally, display the set of lines
    mlab.view(elevation=0) # And choose a nice view
    mlab.roll(125) # Again, choose a nice view
    mlab.show()

def consolidate_collision_dicts(list_of_electrons):
    """I believe I wrote this function to view the overall number of collisions that happend for each collision type.
    If memory serves me write, I did verify the expected ratios of collisions. However as of now, the function has
    served it's purpose and is no longer needed by the overall script. If any future user would wish to use this
    function, they'd have to figure out how to place unique markers using mayavi, or downgrade the simulation to 2D
    and use the trusty matplotlib markers.

    Args:
        list_of_electrons: The list of all the electrons propagated.

    Returns:
        complete: A dictionary of electrons that experienced that type of collision at any point in their path.
    """
    complete = dict.fromkeys(my_csv_dict.keys(), [])
    dict.fromkeys(my_csv_dict.keys(), [])
    for elec in list_of_electrons: # For electron in the list of electrons
        current = elec.collisions # Current dictionary equals electron.collisions_dictionary
        for key in current: # Each 'key' is a type of collision that electron experienced (besides elastic or null)
            if key != "elastic" and key != "null":
                list_of_collisions = current[key]
                for i in range(len(new_list)):
                    complete[key].append((list_of_collisions[i][0], list_of_collisions[i][1]))
    return complete

def consolidate_energy_lists(list_of_electrons):
    """Combines the energies from all the electrons popagated into one list.

    Args:
        list_of_electrons: Is the list of all electrons propagated.

    Returns:
        EndArray: Array with every single energy value experienced at every electron step.
    """
    EndArray = []
    for elec in list_of_electrons:
        EndArray += elec.list_e[:]
    return EndArray

def find_corresponding_bin(eV_value, SizeOfBin):
    """Returns a tuple that represents the bin. 
    If our bin size is 0.5, a value of 0.75 would correspond to the 0.5-1.0 bin.

    Args:
        eV_value: The value we are looking to bin.
        SizeOfBin: The size of the bins we are using.

    Returns:
        A tuple that represents the start and end points of the bin.
    """
    current = floor(eV_value)
    while current + SizeOfBin <= eV_value:
        current += SizeOfBin
    return ir(current), ir(current+SizeOfBin)

def bin_energy_values(list_all_eV_values, bin_size = BinSize, x_limit=eedf_x_limit):
    """Bins energy values (in electron volts) from a list into a dictionary. The key is the 
    in-between value. For example, if the bin size is 1, the key to the dictionary for the 0.0 to 1.0
    bin will be 0.5. Bins will be tuples where the 2nd value (index = 1) is non-inclusive. Example below
    
    example_list_all_eV_values = [0.5, 1.1, 1.2, 1.3, 2.0, 2.1]
    Bins for this would be as follows for bin_size = 1.0
    (0, 1): 1
    (1, 2): 3
    (2, 3): 2

    So the (0, 1) bin is really from 0 to 0.999. The (1, 2) bin is really from 1 to 1.999. And so on. We use maximum
    energy that occurred in the simulation as a bin cut-off point. For simplicity, we round to the nearest whole number.

    Args:
        list_all_eV_values: List of all the energy values recorded.
        bin_size: Size of the bins. 
        x_limit: If we want to cut off binning at a certain value (because of negligible magnitudes, use this parameter).

    Returns:
        bins: A dictionary that reflects the electron energy distribution.
    """

    MaxWholeNum = ceil(max(list_all_eV_values))
    unseparated_bins = np.arange(0, MaxWholeNum+bin_size, bin_size) # Creates list of floats separated by bin size

    """The dictionary key-value pairs below will hold our separated bins. We'll name this dictionary 'bins'.
    Now we write a for loop that will go through the list of un-separated bins and turns adjacent values into bins.
    For example [0, 0.5, 1, 1.5, 2] would return a dictionary who's keys look like tuples and correspond to a value 0.
    bins[(0, 0.5)] = 0, bins[(0.5, 1.0)] = 0, and so forth. For now these bins are empty but we will add to them.
    Once again, we run into the python issue of processing floats. It even gets slightly in the way when binning our
    values. So what we're going to do is create the dirty_bins dictionary using a one-line for loop that parses
    through our unseparated_bins bins numpy array. Then we create the bins dictionary that parses through dirty_bins
    and uses the ir() function to clean up those pesky run-on float values."""

    dirty_bins = dict(((unseparated_bins[a], unseparated_bins[a+1]), 0) for a in range(len(unseparated_bins)-1))
    bins = dict(((ir(keyy[0]), ir(keyy[1])), 0) for keyy in dirty_bins)
    # Now that we've created all these empty bins we add to them by using a for loop to parse though all the eV values.
    for energy in list_all_eV_values:
        # The next if statement just filter sout and adjust our data for technicalities.
        # If there is a NaN in our data, we just skip is, as we do not know if it is 0 or infinity.
        if not math.isnan(energy):
            bins[find_corresponding_bin(energy, bin_size)] += 1
    return bins

def plot_energy_distribution_function(list_all_eV_values, x_limit = eedf_x_limit):
    """Plots the Electron Energy Distribution function

    Args:
        list_all_eV_values: All energy values experienced at any step.
        x_limit: Energy cutoff point (b/c negligibility of values).
    
    """
    # Creating figure
    fig = plt.figure()
    energy_bins = bin_energy_values(list_all_eV_values)
    list_keys = list()
    list_vals = list()
    for key in energy_bins:
        inBetween = mean([key[0], key[1]])
        if (inBetween < x_limit) and (energy_bins[key] > 0):
            list_keys.append(inBetween)
            list_vals.append(energy_bins[key])

    # plt.scatter(list_keys, list_vals, s=5)
    plt.plot(list_keys, list_vals)
    plt.ylabel("# of Electrons w/ that energy at any point")
    plt.xlabel("Energy (eV)")
    fig.savefig("NumElectronsvsEnergy.pdf")

def write_eedf_spreadsheet(all_eV_values):
    """Writes values from bins into a Microsoft Excel format spreadsheet.

    Args:
        all_eV_values: All energy values experienced at all steps by all electrons.
    """
    silent_remove('EnergyDistribution.xlsx') # Delete old file if it exists.
    workbook = xlsxwriter.Workbook('EnergyDistribution.xlsx') # Create the workbook
    worksheet = workbook.add_worksheet() # Create the worksheet
    CurrentRow = 0 # Current row being written to.
    ener_bins = bin_energy_values(all_eV_values)
    for key in ener_bins:
        InBetween = (key[0] + key[1])*0.5 # Midpoint value of a bin is used to represent said bin
        worksheet.write(CurrentRow, 0, InBetween)
        worksheet.write(CurrentRow, 1, ener_bins[key])
        CurrentRow += 1
    workbook.close()

def write_normalized_eedf_spreadsheet(all_eV_Values):
    """Writes values from bins into a Microsoft Excel format spreadsheet
    First we find the sum of all the electrons in every bin"""
    sum = 0
    binSize = 0.5
    ener_bins = bin_energy_values(all_eV_Values)
    for key in ener_bins:
        sum += ener_bins[key]
    # Now we check if the files already exists. If so we delete it. (If you do not delete it Python will continue
    # to write to an old file.
    silent_remove('NormalizedEnergyDistribution.xlsx')
    # Now we create the workbook
    workbook = xlsxwriter.Workbook('NormalizedEnergyDistribution.xlsx')
    # Next we create the worksheet
    worksheet = workbook.add_worksheet()
    # CurrentRow is the row the script is currently writing into. Like most things in programming, we start at 0.
    CurrentRow = 0
    for key in ener_bins:
        # The variable below, InBetween, is just what it sounds like. To simplify graphing and spreadsheet writing,
        # the midpoint value of our bin is used to represent said bin. For example, 0.5 represents (0, 1.0),
        # 0.25 represents (0.0, 0.5), 0.15 represents (0.10, 0.20), and so on.
        InBetween = (key[0] + key[1])*0.5
        worksheet.write(CurrentRow, 0, InBetween)
        worksheet.write(CurrentRow, 1, ener_bins[key]/sum)
        CurrentRow += 1
    workbook.close()

def sum_excitation_collision_frequencies(all_collision_frequencies):
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

def condense_collision_frequency_dictionary(all_collision_frequencies):
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
    condensed_collision_frequencies["sumExcite"] = sum_excitation_collision_frequencies(all_collision_frequencies)
    return condensed_collision_frequencies

def CSvE2(iterALL, exciteLabel=False, SumExcite = True):
    # For simplicity purposes, going to home directory.
    return_to_start_directory()
    # Creating Figure
    fig = plt.figure()
    # iterEV's are all the interpolated electron volt values from 0.00 to 100.00 in intervals of 0.00.
    iterEV = return_iter_evs()

    ionLST = [iterALL["ionization"][elem] for elem in iterEV]
    elasLST = [iterALL["elastic"][elem] for elem in iterEV]
    if SumExcite:
        SumExciteDict = sum_excitation_collision_frequencies(iterALL)
        SumExciteLST = [SumExciteDict[elem] for elem in iterEV]
        plt.ylim(0, max(ionLST + elasLST + SumExciteLST))
    plt.plot(iterEV, elasLST, label="Elastic",linestyle="-")
    if type(exciteLabel) == type(str("string")):
        exciteLST = [iterALL[exciteLabel][elem] for elem in iterEV]
        plt.plot(iterEV, exciteLST, label=exciteLabel, linestyle=":")
    if SumExcite:
        plt.plot(iterEV, SumExciteLST, label="Sum of Excited States", linestyle=":")
    plt.plot(iterEV, ionLST, label="Ionization")
    plt.xlim(0, eedf_x_limit)
    plt.ylabel("Cross Section ($m^2$)")
    plt.xlabel("Energy (eV)")
    lgd = plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    # plt.title("Cross Section vs Energy")
    return_to_start_directory()
    # This is the directory of all the results.
    results_directory = os.getcwd() + "/Results"
    # If it does not currently exist, the next two lines will create it using an if-statement to check if it exists.
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    os.chdir(results_directory)
    fig.savefig("Cross Section vs Energy.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    return_to_start_directory()
    return

def nonItered_CSvE(iterALL, exciteLabel=False, SumExcite = True):
    return_to_start_directory() # For simplicity purposes, going to home directory.
    fig = plt.figure()     # Creating Figure
    iterEV = return_iter_evs() # Interpolated electron volt values from 0.00 to 100.00 in intervals of 0.00
    evs = return_evs().tolist()
    ni_ion = [iterALL["ionization"][elem] for elem in evs]
    ni_elastic = [iterALL["elastic"][elem] for elem in evs]
    ionLST = [iterALL["ionization"][elem] for elem in iterEV]
    elasLST = [iterALL["elastic"][elem] for elem in iterEV]
    MaxY = max([max(ionLST), max(elasLST)])
    if SumExcite:
        SumExciteDict = sum_excitation_collision_frequencies(iterALL)
        SumExciteLST = [SumExciteDict[elem] for elem in iterEV]
        plt.ylim(0, max(ionLST + elasLST + SumExciteLST))

    plt.scatter(evs, ni_elastic , label="Elastic", c='#1f77b4', marker='x',s=20)
    plt.plot(evs, ni_elastic, c='#1f77b4')
    if type(exciteLabel) == type(str("string")):
        exciteLST = [iterALL[exciteLabel][elem] for elem in iterEV]
        ni_excite = [iterALL[exciteLabel][elem] for elem in evs]
        plt.scatter(evs, ni_excite, label=exciteLabel, marker='x',c='#2ca02c')
        plt.plot(evs, ni_excite, c='#2ca02c')

        MaxY = max([max(ionLST), max(exciteLST), max(exciteLST)])

    if SumExcite:
        plt.scatter(iterEV, SumExciteLST, label="Sum of Excited States")
        MaxY = max([max(ionLST), max(exciteLST), max(SumExciteLST)])

    plt.scatter(evs, ni_ion , label="Ionization", c='orange', marker='x',s=20)
    plt.plot(evs, ni_ion , c='orange', linewidth=1.0)
    plt.ylim(0, MaxY * 1.025)
    plt.xlim(0, eedf_x_limit)
    plt.ylabel("Cross Section ($m^2$)")
    plt.xlabel("Energy (eV)")
    lgd = plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    # plt.title("Cross Section vs Energy")
    return_to_start_directory()
    results_directory = os.getcwd() + "/Results" # This is the directory of all the results.
    # If it does not currently exist, the next two lines will create it using an if-statement to check if it exists.
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    os.chdir(results_directory)
    fig.savefig("Cross Section vs Energy.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    return_to_start_directory()
    plt.show()
    return

def CFvE2(allCF = create_cf_dicts(), SumExcite = True, IndivExcite = False, Null = True):
    return_to_start_directory() # Going back to home directory
    ionLST = [allCF["ionization"][elem] for elem in interpolated_ev]
    elasLST = [allCF["elastic"][elem] for elem in interpolated_ev]
    if type(IndivExcite) ==  type(str("string")):
        exciteLST = [exciteCfDict[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, exciteLST, label=IndivExcite, linestyle=":")
    fig = plt.figure() # Creating Figure
    plt.plot(interpolated_ev, elasLST, label="Elastic", linestyle="--")
    plt.plot(interpolated_ev, ionLST, label="Ionization", linestyle="-.")
    TOTLST = [TotalCF4ThisEnergy[elem] for elem in interpolated_ev]
    MaxY = max(TOTLST)
    plt.plot(interpolated_ev, TOTLST, label="Total Collision Frequency", linestyle="-", c="m")
    if Null:
        # The following 2 lines will plot the null collision frequency as a function of energy
        nullLST = [NCF4Energy[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, nullLST, label="Null", linestyle="-", c="black")
        # The following 4 lines will plot the total + null collision frequency as a function of energy
        TotWithNullLST = [(TotalCF4ThisEnergy[elem] + NCF4Energy[elem]) for elem in interpolated_ev]
        plt.plot(interpolated_ev, TotWithNullLST, label="Total w/ Null", linestyle="-", c="r")
        MaxY = max(TotWithNullLST)
    if SumExcite:
        SumExciteCfDict = sum_excitation_collision_frequencies(allCF)
        SumExciteList = [SumExciteCfDict[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, SumExciteList, label="Excitation", linestyle="-", color="g")
    plt.ylim(0, MaxY * 1.025)
    plt.xlim(0, 40)
    plt.ylabel('Collision Frequency ($s^{-1}) $')
    plt.xlabel("Energy (eV)")
    # fontP = FontProperties()
    # fontP.set_size('small')
    lgd = plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    # plt.title("Collision Frequency vs Energy")
    chdir(getcwd() + "/Results")
    if Null:
        fig.savefig("CFvEnergy with Null.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        fig.savefig("CFvEnergy.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    return_to_start_directory()

def computeSpecificCFs(All_CF_Dict, TotCF_Dict, writeSpreadSheet = True):
    return_to_start_directory()
    results_directory = os.getcwd() + "/Results" # This is the directory of all the results.
    silent_create(results_directory)
    silent_remove("SpecificCollisionFrequencies.xlsx")
    SixCollisionFrequencies = condense_collision_frequency_dictionary(All_CF_Dict) # We can have up to six dictionaries
    SixCollisionFrequencies["total"] = totalcfs_and_max2(AllIterCS)
    SixCollisionFrequencies["null"] = NCF4Energy
    TotWithNullLST = [(TotalCF4ThisEnergy[elem] + NCF4Energy[elem]) for elem in interpolated_ev]
    SixCollisionFrequencies["tot+null"] = dict(zip(interpolated_ev, TotWithNullLST))
    CollisionColumn = dict()
    CollisionColumn["elastic"] = "B"
    CollisionColumn["sumExcite"] = "C"
    CollisionColumn["ionization"] = "D"
    CollisionColumn["total"] = "E"
    CollisionColumn["null"] = "F"
    CollisionColumn["tot+null"] = "G"

    CheckValues = [10, 25] # Some eV values of interest to us
    if (writeSpreadSheet):
        workbook = xlsxwriter.Workbook("SpecificCollisionFrequencies.xlsx") # This will be the excel OrigValuesSheet
        ScaledValuesSheet = workbook.add_worksheet("ScaledValues")
        ScaledValuesSheet.write("A1", "Energy (eV)")
        for cType  in CollisionColumn:
            ScaledValuesSheet.write("{0}1".format(CollisionColumn[cType]), cType)
        ValueRow = dict()
        for i in range(len(CheckValues)):
            ScaledValuesSheet.write("A"+str(i+3), CheckValues[i])
            ValueRow[CheckValues[i]] = str(i+3)
        scale = 10**(-8)
        ScaledValuesSheet.write("B2", r"($\times10^8 s^{-1}$)")
        ScaledValuesSheet.write("C2", r"($\times10^8 s^{-1}$)")
        ScaledValuesSheet.write("D2", r"($\times10^8 s^{-1}$)")
        ScaledValuesSheet.write("E2", r"($\times10^8 s^{-1}$)")
        ScaledValuesSheet.write("F2", r"($\times10^8 s^{-1}$)")
        ScaledValuesSheet.write("G2", r"($\times10^8 s^{-1}$)")

        for value in CheckValues: # The outer layer of the for loop iterate through the eV values we desire to check. .
            totCF = TotCF_Dict[value] # Calculating total CF for this energy
            ScaledValuesSheet.write("F"+ValueRow[value], totCF*scale) # Writing it to OrigValuesSheet
            # The next layer of the for loop will iterate through all the collision types. In other words, the columns.
            for spec_cf in SixCollisionFrequencies:
                ScaledCF = SixCollisionFrequencies[spec_cf][value]*scale
                ScaledValuesSheet.write(str(CollisionColumn[spec_cf] + ValueRow[value]), float(format(ScaledCF, '.3f')))
        workbook.close()


AllIterCS = returnIterCsDictionary()
sum_excitation_collision_frequencies(AllIterCS)
# CSvE2(AllIterCS, exciteLabel="s3p10", SumExcite=False)
# nonItered_CSvE(AllIterCS, exciteLabel="s3p10", SumExcite=False)
# nonItered_CSvE(AllIterCS, exciteLabel="s3p10",SumExcite=False)
# CFvE2()
# CFvE2(Null=True)
tcf_dict = totalcfs_and_max2(AllIterCS)
AllIterCF = create_cf_dicts()
computeSpecificCFs(AllIterCF, tcf_dict, writeSpreadSheet=False)

def probabilities(allCF = create_cf_dicts(), SumExcite = True, IndivExcite = False, Null = True):
    return_to_start_directory() # Going back to home directory
    ionLST = [allCF["ionization"][elem] for elem in interpolated_ev]
    elasLST = [allCF["elastic"][elem] for elem in interpolated_ev]
    if Null:
        TotWithNullLST = [(TotalCF4ThisEnergy[elem] + NCF4Energy[elem]) for elem in interpolated_ev]
        MaxY = max(TotWithNullLST)
        TotWithNullLST = np.true_divide(np.asarray(TotWithNullLST), MaxY)
        ionLST = np.true_divide(np.asarray(ionLST), MaxY)
        elasLST = np.true_divide(np.asarray(elasLST), MaxY)
        nullLST = np.true_divide(np.asarray([NCF4Energy[elem] for elem in interpolated_ev]), MaxY)

    if type(IndivExcite) == type(str("string")):
        exciteLST = [exciteCfDict[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, exciteLST, label=IndivExcite, linestyle="-")

    fig = plt.figure()# Creating Figure

    plt.plot(interpolated_ev, elasLST, label="Elastic", linestyle="-")
    plt.plot(interpolated_ev, ionLST, label="Ionization", linestyle="-", linewidth=0.75)

    TOTLST = [TotalCF4ThisEnergy[elem] for elem in interpolated_ev]
    TOTLST = np.asarray(TOTLST)
    TOTLST = np.true_divide(TOTLST, MaxY)
    if Null:
        plt.plot(interpolated_ev, TotWithNullLST, label="Total", linestyle="-", c="r") # Total (including NCF) vs Energy
        plt.fill_between(interpolated_ev, elasLST, TotWithNullLST, color='#f08b8d', alpha=0.5)
    if SumExcite:
        SumExciteCfDict = sum_excitation_collision_frequencies(allCF)
        SumExciteList = [SumExciteCfDict[elem] for elem in interpolated_ev]
        SumExciteList = np.asarray(SumExciteList)
        SumExciteList = np.true_divide(SumExciteList, MaxY)
        plt.plot(interpolated_ev, SumExciteList, label="Excitation", linestyle="-", color="g", linewidth=0.75)
        plt.fill_between(interpolated_ev, np.maximum(ionLST, SumExciteList), elasLST, color='#87ceeb', alpha=0.5)
        plt.fill_between(interpolated_ev, SumExciteList, np.maximum(SumExciteList, ionLST), color='#f2d099', alpha=0.5)
        plt.fill_between(interpolated_ev, np.minimum(ionLST, SumExciteList), SumExciteList, color='grey', alpha=0.5)
        plt.fill_between(interpolated_ev, np.zeros(len(interpolated_ev)), np.minimum(SumExciteList, ionLST),
                         color='lightgreen', alpha=0.5)

    p_null = mpatches.Patch(color='#f08b8d', alpha=0.5, label='Null')
    p_elastic = mpatches.Patch(color='#87ceeb', label='Elastic')
    p_ion = mpatches.Patch(color='#f2d099', label='Ionization')
    p_excite = mpatches.Patch(color='lightgreen', label='Excitation')

    lgd = plt.legend(handles=[p_null, p_elastic, p_ion, p_excite], loc="upper left", bbox_to_anchor=(1, 1))
    plt.ylim(0, 1.0085)
    plt.xlim(0, 27.5)
    plt.ylabel('Cumulative Collision Probability')
    plt.xlabel("Energy (eV)")
    chdir(getcwd() + "/Results")
    if Null:
        fig.savefig("Probabilities including Null.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        fig.savefig("Probabilities excluding Null.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    return_to_start_directory()


# def reduced_ev_vs_electron_temp(list_reduced_evs, num_electrons=NOE, num_motions=MPE, bin_size=0.5):
#     fig = plt.figure()
#     ev_values = np.arange(bin_size/2, 100, bin_size)
#     eedfs = pd.DataFrame(np.array(np.arange(0.25, 100, 0.50)), columns = ['Energy_in_eV'])
#     energy_values = eedfs['Energy_in_eV'].values
#
#     eedfs = pd.DataFrame((np.arange(0.25, 100, 0.5), columns='energy'))
#     for eV in list_reduced_evs:
#         eedfs[str(eV) + "_Td"] = pd.Series(np.zeros(len(eedfs.index)), index=eedfs.index)
#         current_energies = consolidate_energy_lists(run_electrons(num_electrons, num_motions))
#         ener_bins = bin_energy_values(current_energies)
#         for key in ener_bins:
#             current_energy = key[1] - bin_size*0.5
#             eedfs.loc[eedfs['Energy_in_eV'] == current_energy, str(eV) + "_Td"] \
#                 = ener_bins[(current_energy-bin_size*0.5, current_energy+bin_size*0.5)]
#
#         eedf_values = eedfs[str(eV) + "_Td"].values
#         simpson = integrate.simps(2 * energy_values * eedf_values / (3 * sum(eedf_values) * 0.5), energy_values)
#         plt.scatter(eV, simpson, color='red')
#
#
#     plt.xlabel("Reduced Electric Field (Td)")
#     plt.ylabel("Electron Temperature (eV)")
#     plt.savefig("Reduced Electric Field vs Electron Temperature.png")
#     plt.show()

def draw_probability_graph():
    probabilities(Null=True)
    plt.show()

draw_probability_graph()