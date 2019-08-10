import os
import ntpath as nt
import numpy, math
from numpy import loadtxt
import csv
import errno
import time
from src.initial_conditions import *
import operator
import shutil


start_directory = os.getcwd()
def return_to_start_directory():
    os.chdir(start_directory)

def corresponding_atom_density(collision_type):
    """
    Returns the corresponding atom density from a collision type.

    The atom density for ground state argon atoms is 10^4 times larger than that of atoms already in an excited state.
    So when determining the collision frequency, we must use the corresponding atom volume density.

    Args:
        collision_type: The type of collision that happened (i.e. gp1, elastic, ionization, s5p10, etc.).

    Returns:
        Corresponding atom density as a float.
    """
    if collision_type[0] == "s" or collision_type == "electron":
        return metastable_state_atom_density
    return ground_state_atom_density

def collision_frequency(ener_eV, cs_dictionary, collision_type):
    """Returns collision frequency for a given energy value. If the given energy is greater than 100.00, eV, it will
    be rounded down to 100.00 eV (also sets given NaN values as 100.00 eV).

    Args:
        ener_eV: The energy of the electron in Electron Volts.
        cs_dictionary: Dictionary of all electron cross section values for all possible collisions.
        collision_type: The type of collision that occurred (i.e. elastic, ionization, gp1, etc.).

    Returns:
        cf: The corresponding collision frequency as a float.

     """
    if numpy.isnan(ener_eV) or ener_eV > 100.00:
        ener_eV = 100.00
    cf = corresponding_atom_density(collision_type)*(cs_dictionary[ir(ener_eV)])*math.sqrt(2*ir(ener_eV)*EV2J/m_e)
    return cf


def silent_remove(path, is_directory=False):
    """Removes a file without raising errno.ENOENT error if the file doesn't exist

    Args:
        path: The file or directory path as a string.
        is_directory: Boolean value reflecting if path is a directory. Default set to false, as we assume files.
    Raises:
        Exception: Any error/exception other than 'no such file or directory'.
    """
    try:
        if is_directory and os.path.exists(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


def silent_create(directory_name):
    """Creates a directory if it does not yet exist

      Args:
        directory_name: Intended directory name as astring.

      Returns:
        True if directory was created, False if it already existed beforehand."""
    if os.path.exists(directory_name):
        return False
    os.makedirs(directory_name)
    return True

def returnInterpolatedCsFromDAT(energy_eV_dat_file, ui_cross_section_dat_file):
    """Returns interpolated Cross Sections as a dictionary from a DAT file

    Args:
        energy_eV_dat_file: xsigma.dat
        ui_cross_section_dat_file: Cross section dat file for specific collision type

    Returns:
        The energy-probability dictionary for that collision type
    """
    return_to_start_directory()
    os.chdir(theoretical_cs_directory)
    energy_levels = loadtxt(energy_eV_dat_file)
    uninterpolated_cross_sections = loadtxt(ui_cross_section_dat_file) # ui stands for un-interpolated
    all_eV_levels = []
    interpolated_cross_sections = [] # interpolated_cross_sections stands for interpolated cross sections
    for i in range(len(energy_levels) -1):
        lower_limit = energy_levels[i]
        lower_cs = uninterpolated_cross_sections[i]
        upper_limit = energy_levels[i+1]
        upper_cs = uninterpolated_cross_sections[i+1]
        if i == 0:
            interpolated_energy = [lower_limit] + list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
            interpolated_sigma = [lower_cs] + list(numpy.linspace(lower_cs, upper_cs, 201)[1:200]) + [upper_cs]
        else:
            interpolated_energy = list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
            interpolated_sigma = list(numpy.linspace(lower_cs, upper_cs, 201)[1:200]) + [upper_cs]
        all_eV_levels = all_eV_levels[:] + interpolated_energy
        interpolated_cross_sections = interpolated_cross_sections[:] + interpolated_sigma
    for i in range(len(all_eV_levels)):
        all_eV_levels[i] = ir(all_eV_levels[i])
    energy_probability_dictionary = dict(zip(all_eV_levels, interpolated_cross_sections))
    if ui_cross_section_dat_file == "sigmaelastic.dat":
        elasCSatTwoEv = energy_probability_dictionary[2.00]
        FixedIntialCsValues = numpy.arange(0, elasCSatTwoEv,elasCSatTwoEv/200)
        for enerKey in energy_probability_dictionary:
            if enerKey < 2.00:
                energy_probability_dictionary[enerKey] = FixedIntialCsValues[int(enerKey*100)]
    return_to_start_directory()
    return energy_probability_dictionary

def dat_to_diction(energy_eV_dat_file, ui_cross_section_dat_file):
    """Turns a DAT File in to a dictionary.

    Args:
        energy_eV_dat_file: xsigma.dat
        ui_cross_section_dat_file: Cross section dat file for specific collision type

    Returns:
        The energy-probability dictionary for that collision type
    """
    return_to_start_directory()
    os.chdir(os.getcwd() + "/CrossSectionsTheoretical")
    energy_levels = loadtxt(energy_eV_dat_file)
    uninterpolated_cross_sections = loadtxt(ui_cross_section_dat_file) # ui stands for un-interpolated
    all_eV_levels = []
    interpolated_cross_sections = []
    for i in range(len(energy_levels) -1):
        lower_limit = energy_levels[i]
        lower_cs = uninterpolated_cross_sections[i]
        upper_limit = energy_levels[i+1]
        upper_cs = uninterpolated_cross_sections[i+1]
        if i == 0:
            interpolated_energy = [lower_limit] + list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
            interpolated_sigma = [lower_cs] + list(numpy.linspace(lower_cs, upper_cs, 201)[1:200]) + [upper_cs]
        else:
            interpolated_energy = list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
            interpolated_sigma = list(numpy.linspace(lower_cs, upper_cs, 201)[1:200]) + [upper_cs]
        all_eV_levels = all_eV_levels[:] + interpolated_energy
        interpolated_cross_sections = interpolated_cross_sections[:] + interpolated_sigma
    for i in range(len(all_eV_levels)):
        all_eV_levels[i] = ir(all_eV_levels[i])
    energy_probability_dictionary = dict(zip(all_eV_levels, interpolated_cross_sections))
    if ui_cross_section_dat_file == "sigmaelastic.dat":
        elasCSatTwoEv = energy_probability_dictionary[2.00]
        FixedIntialCsValues = numpy.arange(0, elasCSatTwoEv,elasCSatTwoEv/200)
        for enerKey in energy_probability_dictionary:
            if enerKey < 2.00:
                energy_probability_dictionary[enerKey] = FixedIntialCsValues[int(enerKey*100)]
    return_to_start_directory()
    return energy_probability_dictionary

def create_interpolated_cs_csv_file(collision_type, interpolated_dictionary):
    """Creates interpolated Cross Section CSV File. If you're going to be running this script often, it's 

    Args:
        collision_type: Type of collision.
        interpolated_dictionary: Is the interpolated dicionary we are turning into a CSV file.

    """
    return_to_start_directory()
    if not os.path.exists(os.getcwd()+ "/InterpolatedCrossSections"):
        os.makedirs(os.getcwd()+ "/InterpolatedCrossSections")
        os.chdir(os.getcwd()+ "/InterpolatedCrossSections")
    else:
        os.chdir(os.getcwd()+ "/InterpolatedCrossSections")
    with open('i_{0}.csv'.format(collision_type), 'w', newline='') as csvfile:
        for inter_energy in interpolated_dictionary:
            iwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            iwriter.writerow([inter_energy, interpolated_dictionary[inter_energy]])
    return_to_start_directory()

def all_cs_theoretical():
    """Takes all Cross Sections and combines them into a dictionary of dictionaries.

    Returns:
        dictionary_of_ds: A dictionary of energy-collision frequency dictionaries.
    """
    return_to_start_directory()
    dictionary_of_ds = dict()
    if os.path.exists(os.getcwd() + "/InterpolatedCrossSections"):
        os.chdir(os.getcwd() + "/InterpolatedCrossSections")
    else:
        for datFile in os.listdir(theoretical_cs_directory):
            if (".DS_Store" not in os.path.realpath(datFile)) and ("xsigma.dat" not in os.path.realpath(datFile)):
                collision_type = nt.basename(os.path.realpath(datFile)).replace("sigma", "").replace(".dat", "")
                csDictionary = returnInterpolatedCsFromDAT("xsigma.dat", nt.basename(os.path.realpath(datFile)))
                dictionary_of_ds[collision_type] = csDictionary
    return_to_start_directory()
    return dictionary_of_ds


def ir(ener):
    """Due to issues in Floating Point Arithmetic, I wrote the function ir() in order to help move between two point
    decimal accuracy floats and strings, as both of these formats were passed around frequently. Long story short,
    if you have any float, ir() will round it to two decimal points of accuracy.

    Args:
        ener: The energy in electron volts you are rounding down to two decimal accuracy.

    Returns:
        That energy as a two point accuracy float (or at least it tries its best).
    """
    return float(str(round(float(ener), 2)))

def returnIterCsDictionary():
    """Returns a dictionary regardless of whether or not the cross sections have been interpolated.
    If they have not been interpolated into csv files, the function will parse and interpolate the .dat files into our
    sought after dictionary. If the interpolated csv files already exist, it will read them into the same dictionary.
    On average, former takes about 6.5 seconds while the latter takes between on 0.9 seconds.
    """
    return_to_start_directory()
    if not os.path.exists(os.getcwd() + "/InterpolatedCrossSections"):
        dod = all_cs_theoretical()
        for col_type in dod:
            create_interpolated_cs_csv_file(col_type, dod[col_type])
        # Now that
        return_to_start_directory()
        return returnIterCsDictionary()
    else:
        # If the aforementioned path does not exist, that means the interpolated CS csv files are already created.
        # So now we read all those csv files into one dictionary.
        return_to_start_directory()
        os.chdir(os.getcwd()+'/InterpolatedCrossSections')
        dictionary_of_ds = dict()
        for csvFile in os.listdir(os.getcwd()):
            interpolated_cs_dict = dict()
            with open(csvFile, mode='r') as infile:
                reader = csv.reader(infile)
                for row in reader:
                    # row will be a list of two strings that looks like: ['12.23', '5.48435e-26']
                    # The first item, index = 0, is the interpolated energy
                    # The second item, index = 1, is the interpolated cross section
                    # The former will be the key to the dictionary, while the latter will be the value.
                    interpolated_cs_dict[float(row[0])] = float(row[1])
            dictionary_of_ds[csvFile.replace("i_", "").replace(".csv", "")] = interpolated_cs_dict
        return_to_start_directory()
        return dictionary_of_ds

interpolatedCsDict = returnIterCsDictionary()

def create_cf_files(dictOfDicts=interpolatedCsDict):
    """Creates the Collision Frequency files from a Dictionary of Dictionaries
    Args:
        dictOfDicts: Dictionary of Collision Frequency dictionaries.
    """
    return_to_start_directory()
    if not os.path.exists(os.getcwd()+ "/CollisionFrequencies"):
        os.makedirs(os.getcwd()+ "/CollisionFrequencies")
        os.chdir(os.getcwd() + "/CollisionFrequencies")
    else:
        os.chdir(os.getcwd() + "/CollisionFrequencies")
    cf_dictionary_ofDictionaries = dict()
    for proce in dictOfDicts:
        cf_dictionary_ofDictionaries[proce] = dict()
        with open('cf_{0}.csv'.format(proce), 'w', newline='') as csvfile:
            iwriter = csv.writer(csvfile, delimiter=',',
                                 quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for en in dictOfDicts[proce]:
                current_cf = collision_frequency(en, dictOfDicts[proce], proce)
                iwriter.writerow([en, current_cf])
                cf_dictionary_ofDictionaries[proce.replace("cf_", "")][en] = current_cf
    return_to_start_directory()
    return cf_dictionary_ofDictionaries


def create_cf_dicts():
    """
    Creates collision frequency dictionaries from pre-existing files.

    Returns:
        dictionary_of_ds: Dictionary of dictionaries.
    """
    # If path does not exist, then we have to create the files and can return the dictionary while we're at it
    return_to_start_directory()
    if not os.path.exists(os.getcwd() + "/CollisionFrequencies"):
        return create_cf_files()
    # If path exists returns True that means the interpolated CS csv files are already created.
    # So now we read all those csv files into one dictionary.
    else:
        return_to_start_directory()
        os.chdir(os.getcwd()+'/CollisionFrequencies')
        dictionary_of_ds = dict()
        for csvFile in os.listdir(os.getcwd()):
            interpolated_cs_dict = dict()
            with open(csvFile, mode='r') as infile:
                reader = csv.reader(infile)
                for row in reader:
                    # row will be a list of two strings that looks like: ['12.23', '5.48435e-26']
                    # The first item, index = 0, is the interpolated energy
                    # The second item, index = 1, is the interpolated cross section
                    # The former will be the key to the dictionary, while the latter will be the value.
                    interpolated_cs_dict[float(row[0])] = float(row[1])
            dictionary_of_ds[csvFile.replace("cf_", "").replace(".csv", "")] = interpolated_cs_dict
        return_to_start_directory()
    return dictionary_of_ds


def return_interpolated_ev(energy_eV_dat_file):
    """
    Returns interpolated energy values from xsigma.dat file.

    Args:
        energy_eV_dat_file: Should be xsigma.dat file.

    Returns:
        all_eV_levels: Are all the energies from xsigma.dat, but instead of going by intervals of two, goes by
        intervals of 0.01.
    """
    return_to_start_directory()
    os.chdir(theoretical_cs_directory)
    energy_levels = loadtxt(energy_eV_dat_file)
    all_eV_levels = []
    for i in range(len(energy_levels) -1):
        lower_limit = energy_levels[i]
        upper_limit = energy_levels[i+1]
        if i == 0:
            interpolated_energy = [lower_limit] + list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
        else:
            interpolated_energy = list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
        all_eV_levels = all_eV_levels[:] + interpolated_energy
    for i in range(len(all_eV_levels)):
        all_eV_levels[i] = ir(all_eV_levels[i])
    return_to_start_directory()
    return all_eV_levels

def return_iter_evs():
    """
    Returns interpolated energy values from xsigma.dat file and deals with directory navigation.
    """
    os.chdir(theoretical_cs_directory)
    i_ev = return_interpolated_ev("xsigma.dat")
    os.chdir(os.getcwd().replace("/CrossSectionsTheoretical", ""))
    return i_ev

def return_evs():
    """
    Returns interpolated energy values from xsigma.dat file and deals with directory navigation.
    """
    return_to_start_directory()
    os.chdir(os.getcwd()+"/CrossSectionsTheoretical")
    energy_levels = loadtxt("xsigma.dat")
    os.chdir(os.getcwd().replace("/CrossSectionsTheoretical", ""))
    return energy_levels

interpolated_ev = return_iter_evs()

def float_in_range(flt, start, stop):
    if start <= flt <= stop:
        return True
    return False


# The next function, total_CF will return total collision frequency for a given energy
def total_CF(energy_in_eV, dictionary_of_cs_dictionaries):
    tot = 0
    for process in dictionary_of_cs_dictionaries:
        tot += collision_frequency(energy_in_eV, dictionary_of_cs_dictionaries[process], process)
    return tot

def totalcfs_and_max(returnIncidentEnergy=False):
    """ This function defines two things in one go
    First it finds the Maximum Total Collision Frequency and the value it occurs at
    Second it creates a dictionary where they key is an energy between 0 and 100 and the returned value
    is the the Total Collision Frequency for that energy

    Args:
        returnIncidentEnergy:

    Returns:
        value:
        EnergyTotalCF
    """
    EnergyTotalCF = dict() # Creating Dictionary
    allTot = [] # Empty list for sole purpose of finding list max using built in python max()
    interpolated_energy_list = [] # Corresponding energy list so that index of energy matches index of TCF
    now = 0.00
    while now < 100.01:
        t = total_CF(ir(now), interpolatedCsDict)
        interpolated_energy_list.append(now)
        allTot.append(t)
        EnergyTotalCF[ir(now)] = t
        now += 0.01
    index, value = max(enumerate(allTot), key=operator.itemgetter(1))
    if returnIncidentEnergy:
        return ir(interpolated_energy_list[index])
    return value, EnergyTotalCF

# MTCF = Maximum Total Collision Frequency

def return_ncf_diction():
    """
    Creates and returns null collision frequency dictionary
    """
    allnull = dict()
    now = 0.00
    while now < 100.01:
        allnull[ir(now)] = maxx - TotalCF4ThisEnergy[ir(now)]
        now += 0.01
    return allnull


def allEV():
    """
    Reads in the xsigma.dat file containing the eV intervals/values at which cross sections were measured.
    """
    chdir(getcwd() + '/CrossSectionsTheoretical')
    all_ev = loadtxt("xsigma.dat")
    chdir(getcwd().replace("/CrossSectionsTheoretical", ""))
    return all_ev




TotalCollisionFrequenciesAndMax = totalcfs_and_max()
maxx = TotalCollisionFrequenciesAndMax[0]

TotalCF4ThisEnergy = TotalCollisionFrequenciesAndMax[1]
NCF4Energy = return_ncf_diction()
