import os
import ntpath as nt
import numpy
from numpy import loadtxt
import csv
import errno
import time
from initialConditions import *
import math
import operator




home = os.getcwd()
def goHomeDirectory():
    os.chdir(home)
    pass

def corresponding_atom_density(collision_type):
    # All ground state collisions start with 'g'
    if collision_type[0] == "g":
        return GSAD
    if collision_type == "elastic" or collision_type == "ionization":
        return GSAD
    if collision_type[0] == "s":
        return MSSAD

def collision_frequency(ener_eV, cs_dictionary, collision_type):
    # If the energy somehow 'ran away' and became a nan value
    # of if the energy is over 100, we're going to treat it like 100 since that's the highest value we have data for
    if numpy.isnan(ener_eV) or ener_eV > 100.00:
        ener_eV = 100.00
    cf = corresponding_atom_density(collision_type)*(cs_dictionary[ir(ener_eV)])*math.sqrt(2*ir(ener_eV)*EV2J/m_e)
    return cf


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def silentCreate(directoryName):
    # If the does not exist, we create it using os.makedirs (make-directory)
    if not os.path.exists(directoryName):
        os.makedirs(directoryName)
        # Then we return True since we did indeed need to create this directory
        return True
    # If the path already exists, we return false, as we did not need to create anything
    else:
        return False

def returnInterpolatedCsFromDAT(energy_eV_dat_file, ui_cross_section_dat_file):
    # ui stands for un-interpolated
    goHomeDirectory()
    os.chdir(os.getcwd() + "/CrossSectionsTheoretical")
    energy_levels = loadtxt(energy_eV_dat_file)
    ui_cross_sections = loadtxt(ui_cross_section_dat_file)
    all_eV_levels = []
    # i_cs stands for interpolated cross sections
    i_cs = []
    for i in range(len(energy_levels) -1):
        lower_limit = energy_levels[i]
        lower_cs = ui_cross_sections[i]
        upper_limit = energy_levels[i+1]
        upper_cs = ui_cross_sections[i+1]
        if i == 0:
            interpolated_energy = [lower_limit] + list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
            interpolated_sigma = [lower_cs] + list(numpy.linspace(lower_cs, upper_cs, 201)[1:200]) + [upper_cs]
        else:
            interpolated_energy = list(numpy.linspace(lower_limit, upper_limit, 201)[1:200]) + [upper_limit]
            interpolated_sigma = list(numpy.linspace(lower_cs, upper_cs, 201)[1:200]) + [upper_cs]
        all_eV_levels = all_eV_levels[:] + interpolated_energy
        i_cs = i_cs[:] + interpolated_sigma
    for i in range(len(all_eV_levels)):
        all_eV_levels[i] = ir(all_eV_levels[i])
    energy_probability_dictionary = dict(zip(all_eV_levels, i_cs))
    if ui_cross_section_dat_file == "sigmaelastic.dat":
        print("Fixing elastic cross sections")
        elasCSatTwoEv = energy_probability_dictionary[2.00]
        FixedIntialCsValues = numpy.arange(0, elasCSatTwoEv,elasCSatTwoEv/200)
        for enerKey in energy_probability_dictionary:
            if enerKey < 2.00:
                energy_probability_dictionary[enerKey] = FixedIntialCsValues[int(enerKey*100)]
    goHomeDirectory()
    return energy_probability_dictionary

def create_interpolated_cs_csv_file(collision_type, interpolated_dictionary):
    goHomeDirectory()
    print("About to write cross section CSV file for:", collision_type)

    if not os.path.exists(os.getcwd()+ "/InterpolatedCrossSections"):
        os.makedirs(os.getcwd()+ "/InterpolatedCrossSections")
        os.chdir(os.getcwd()+ "/InterpolatedCrossSections")
    else:
        os.chdir(os.getcwd()+ "/InterpolatedCrossSections")
    with open('i_{0}.csv'.format(collision_type), 'w', newline='') as csvfile:
        for inter_energy in interpolated_dictionary:
            iwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            iwriter.writerow([inter_energy, interpolated_dictionary[inter_energy]])
    goHomeDirectory()
    pass

def all_cs_theoretical():
    print("starting all cs theoretical function")
    goHomeDirectory()
    if os.path.exists(os.getcwd() + "/InterpolatedCrossSections"):
        print("interpolated cross sections does exist")
        os.chdir(os.getcwd() + "/InterpolatedCrossSections")
    else:
        print("interpolated cross sections folder does not exist yet")
        dictionary_of_ds = dict()
        for datFile in os.listdir(os.getcwd()+"/CrossSectionsTheoretical"):
            if ".DS_Store" not in os.path.realpath(datFile):
                if "xsigma.dat" not in os.path.realpath(datFile):
                    collision_type = nt.basename(os.path.realpath(datFile)).replace("sigma", "").replace(".dat", "")
                    csDictionary = returnInterpolatedCsFromDAT("xsigma.dat", nt.basename(os.path.realpath(datFile)))
                    dictionary_of_ds[collision_type] = csDictionary
    goHomeDirectory()
    return dictionary_of_ds


def ir(ener):
    return float(str(round(ener, 2)))

# The function below will return a dictionary regardless of whether or not the cross sections have been interpolated.
# If they have not been interpolated into csv files, the function will parse and interpolate the .dat files into our
# sought after dictionary. If the interpolated csv files already exist, it will read them into the same dictionary.
# On average, former takes about 6.5 seconds while the latter takes between on 0.9 seconds.
def returnIterCsDictionary():
    goHomeDirectory()
    if not os.path.exists(os.getcwd() + "/InterpolatedCrossSections"):
        print("Interpolated CS folder does not exist yet")
        dod = all_cs_theoretical()
        for col_type in dod:
            create_interpolated_cs_csv_file(col_type, dod[col_type])

        # Now that
        goHomeDirectory()
        return returnIterCsDictionary()
    # If silentCreate returns false that means the interpolated CS csv files are already created.
    # So now we read all those csv files into one dictionary.
    else:
        goHomeDirectory()
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

        goHomeDirectory()
        return dictionary_of_ds

interpolatedCsDict = returnIterCsDictionary()

def create_cf_files(dictOfDicts=interpolatedCsDict):
    goHomeDirectory()
    if not os.path.exists(os.getcwd()+ "/CollisionFrequencies"):
        os.makedirs(os.getcwd()+ "/CollisionFrequencies")
        os.chdir(os.getcwd() + "/CollisionFrequencies")
    else:
        os.chdir(os.getcwd() + "/CollisionFrequencies")
    cf_dictionary_ofDictionaries = dict()
    for proce in dictOfDicts:
        print("About to create an interpolated cross section collision frequency csv file for:", proce)
        cf_dictionary_ofDictionaries[proce] = dict()
        with open('cf_{0}.csv'.format(proce), 'w', newline='') as csvfile:
            iwriter = csv.writer(csvfile, delimiter=',',
                                 quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for en in dictOfDicts[proce]:
                current_cf = collision_frequency(en, dictOfDicts[proce], proce)
                iwriter.writerow([en, current_cf])
                cf_dictionary_ofDictionaries[proce.replace("cf_", "")][en] = current_cf
    goHomeDirectory()
    return cf_dictionary_ofDictionaries


def create_cf_dicts():
    # If path does not exist, then we have to create the files and can return the dictionary while we're at it
    goHomeDirectory()
    if not os.path.exists(os.getcwd() + "/CollisionFrequencies"):
        # print("collision frequencies folder does not exist yet, lets create it")
        # os.makedirs(os.getcwd() + "/CollisionFrequencies")
        return create_cf_files()
    # If path exists returns True that means the interpolated CS csv files are already created.
    # So now we read all those csv files into one dictionary.
    else:
        goHomeDirectory()
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
        goHomeDirectory()
    return dictionary_of_ds


def return_interpolated_ev(energy_eV_dat_file):

    goHomeDirectory()
    os.chdir(os.getcwd()+"/CrossSectionsTheoretical")
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
    goHomeDirectory()
    return all_eV_levels

def return_iter_evs():
    os.chdir(os.getcwd()+"/CrossSectionsTheoretical")
    i_ev = return_interpolated_ev("xsigma.dat")
    os.chdir(os.getcwd().replace("/CrossSectionsTheoretical", ""))
    return i_ev

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


"""
The next function defines two things in one go
First it finds the Maximum Total Collision Frequency and the value it occurs at
Second it creates a dictionary where they key is an energy between 0 and 100 and the returned value
is the the Total Collision Frequency for that energy
"""
def totalcfs_and_max(returnIncidentEnergy=False):
    # Creating Dictionary
    EnergyTotalCF = dict()
    # Creating empty list for purpose of finding list max using built in python max()
    allTot = []
    # Creating corresponding energy list so that index of energy matches index of TCF
    interpolated_energy_list = []
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
    allnull = dict()
    now = 0.00
    while now < 100.01:
        allnull[ir(now)] = maxx - TotalCF4ThisEnergy[ir(now)]
        now += 0.01
    return allnull





"The atom density for ground state argon atoms is 10^4 times larger than that of atoms already in an excited state." \
"So when determining the collision frequency, we must use the corresponding atom volume density."""

"""The corresponding_atom_density function below takes an input of the collision type and outputs the corresponding
 atom volume density."""



def allEV():
    chdir(getcwd() + '/CrossSectionsTheoretical')
    all_ev = loadtxt("xsigma.dat")
    chdir(getcwd().replace("/CrossSectionsTheoretical", ""))
    return all_ev




TotalCollisionFrequenciesAndMax = totalcfs_and_max()
maxx = TotalCollisionFrequenciesAndMax[0]
TotalCF4ThisEnergy = TotalCollisionFrequenciesAndMax[1]
# This dictionary will return null frequency for an energy faster than calling the function each time it is needed
NCF4Energy = return_ncf_diction()

