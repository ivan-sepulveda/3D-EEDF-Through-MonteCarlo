import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import pylab
import matplotlib.legend as lgnd
import numpy as np
from ReturnAllInterpolatedCSDatFiles import *
from math import *
import itertools as itr
from matplotlib.font_manager import FontProperties
from os import *
from returnOutcome import *
import xlsxwriter
import errno
from decimal import Decimal
from collections import OrderedDict
from statistics import mean

J_eV_Factor = 6.241509 * (10 ** 18)


list_colors = ['r', 'b', 'c', 'm']
CycleColors = itr.cycle(list_colors)

def totalcfs_and_max2(allCS, excludeElasIon = False, returnIncidentEnergy=False):
    iterCS = dict(allCS)
    if excludeElasIon:
        del iterCS["elastic"]
        del iterCS["ionization"]
    # Creating Dictionary
    EnergyTotalCF = dict()
    # Creating empty list for purpose of finding list max using built in python max()
    allTot = []
    # Creating corresponding energy list so that index of energy matches index of TCF
    interpolated_energy_list = []
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
    # The number of points per line
    N = len(list_of_electrons[0].list_t)

    # The scalar parameter for each line

    from mayavi import mlab
    # Below i
    mlab.figure(1, size=(400, 400), bgcolor=(0, 0, 0))
    mlab.clf()

    # We create a list of positions and connections, each describing a line.
    # We will collapse them in one array before plotting.
    x_ = list()
    y = list()
    z = list()
    s = list()
    connections = list()

    # The index of the current point in the total amount of points
    index = 0


    # Create each line one after the other in a loop
    for i in range(len(list_of_electrons)):
        x_.append(np.asarray(list_of_electrons[i].list_x))
        y.append(np.asarray(list_of_electrons[i].list_y))
        z.append(np.asarray(list_of_electrons[i].list_z))
        s.append(np.asarray(list_of_electrons[i].list_t))
        # This is the tricky part: in a line, each point is connected
        # to the one following it. We have to express this with the indices
        # of the final set of points once all lines have been combined
        # together, this is why we need to keep track of the total number of
        # points already created (index)
        connections.append(np.vstack(
            [np.arange(index, index + N - 1.5), np.arange(index + 1, index + N - .5)]).T)
        index += N


    # Now collapse all positions, scalars and connections in big arrays
    x_ = np.hstack(x_)
    y = np.hstack(y)
    z = np.hstack(z)
    s = np.hstack(s)
    connections = np.vstack(connections)

    # Create the points
    src = mlab.pipeline.scalar_scatter(x_, y, z, s)

    # Connect them
    src.mlab_source.dataset.lines = connections
    src.update()

    # The stripper filter cleans up connected lines
    lines= mlab.pipeline.tube(src, tube_radius=0.005, tube_sides=6)

    # Finally, display the set of lines
    mlab.pipeline.surface(lines, colormap='Accent', line_width=1, opacity=.4)

    # And choose a nice view
    mlab.view(elevation=0)
    # mlab.roll(125)
    mlab.show()







    return

def consolidate_collision_dicts(list_of_electrons):
    complete = dict()
    complete["elastic"] = []
    complete["ionization"] = []
    complete["gp1"] = []
    complete["gp2"] = []
    complete["gp3"] = []
    complete["gp4"] = []
    complete["gp5"] = []
    complete["gp6"] = []
    complete["gp7"] = []
    complete["gp8"] = []
    complete["gp9"] = []
    complete["gp10"] = []
    complete["s3p1"] = []
    complete["s3p2"] = []
    complete["s3p3"] = []
    complete["s3p4"] = []
    complete["s3p5"] = []
    complete["s3p6"] = []
    complete["s3p7"] = []
    complete["s3p8"] = []
    complete["s3p9"] = []
    complete["s3p10"] = []
    complete["s5p1"] = []
    complete["s5p2"] = []
    complete["s5p3"] = []
    complete["s5p4"] = []
    complete["s5p5"] = []
    complete["s5p6"] = []
    complete["s5p7"] = []
    complete["s5p8"] = []
    complete["s5p9"] = []
    complete["s5p10"] = []
    # For electron in the list of electrons
    for elec in list_of_electrons:
        # Current dictionary equals electron.collisions_dictionary
        current = elec.collisions
        # Each 'key' is a type of collision that electron experienced (besides elastic or null)
        for key in current:
            if key != "elastic":
                new_list = current[key]
                for i in range(len(new_list)):
                    complete[key].append((new_list[i][0], new_list[i][1]))
    return complete

def CombineLists_PullEnergy(list_of_electrons):
    EndArray = []
    for elec in list_of_electrons:
        EndArray += elec.list_e[:]
    return EndArray

def find_bin(eV_value, SizeOfBin):
    current = floor(eV_value)
    while current + SizeOfBin <= eV_value:
        current += SizeOfBin
    return ir(current), ir(current+SizeOfBin)

def bin_eV_values_into_dict(list_all_eV_values, bin_size = BinSize, x_limit=eedf_x_limit):
    # First off, bins will be tuples where the 2nd value (index = 1) is non-inclusive. Example below
    # example_list_all_eV_values = [0.5, 1.1, 1.2, 1.3, 2.0, 2.1]
    # Bins for this would be as follows for bin_size = 1.0
    # (0, 1): 1
    # (1, 2): 3
    # (2, 3): 2
    # So the (0, 1) bin is really from 0 to 0.999.
    # So the (1, 2) bin is really from 1 to 1.999.
    # So the (2, 3) bin is really from 2 to 2.999.
    # Now we find the maximum energy that occurred in the simulation so that we know when to cut off bin creation.
    # For simplicity we round up to the nearest whole number.
    MaxWholeNum = ceil(max(list_all_eV_values))
    # np.arrange will create a list of floats all separated by our bin size. Because all of these floats are listed one
    # after another in a list, we call this list "unseparated_bins".
    unseparated_bins = np.arange(0, MaxWholeNum+bin_size, bin_size)
    # The dictionary key-value pairs below will hold our separated bins. We'll name this dictionary 'bins'.
    # Now we write a for loop that will go through the list of un-separated bins and turns adjacent values into bins.
    # For example [0, 0.5, 1, 1.5, 2] would return a dictionary who's keys look like tuples and correspond to a value 0.
    # bins[(0, 0.5)] = 0, bins[(0.5, 1.0)] = 0, and so forth. For now these bins are empty but we will add to them.
    # Once again, we run into the python issue of processing floats. It even gets slightly in the way when binning our
    # values. So what we're going to do is create the dirty_bins dictionary using a one-line for loop that parses
    # through our unseparated_bins bins numpy array. Then we create the bins dictionary that parses through dirty_bins
    # and uses the ir() function to clean up those pesky run-on float values.
    dirty_bins = dict(((unseparated_bins[a], unseparated_bins[a+1]), 0) for a in range(len(unseparated_bins)-1))
    bins = dict(((ir(keyy[0]), ir(keyy[1])), 0) for keyy in dirty_bins)
    # Now that we've created all these empty bins we add to them by using a for loop to parse though all the eV values.
    for energy in list_all_eV_values:
        # The next if statement just filter sout and adjust our data for technicalities.
        # If there is a NaN in our data, we just skip is, as we do not know if it is 0 or infinity.
        if not math.isnan(energy):
            bins[find_bin(energy, bin_size)] += 1
    return bins
def plot_eV_vs_num_elec(list_all_eV_values, x_limit = eedf_x_limit):
    # Creating figure
    fig = plt.figure()
    enerbins = bin_eV_values_into_dict(list_all_eV_values)
    list_keys = list()
    list_vals = list()
    for key in enerbins:
        inBetween = mean([key[0], key[1]])
        if (inBetween < x_limit) and (enerbins[key] > 0):
            list_keys.append(inBetween)
            list_vals.append(enerbins[key])

    # plt.scatter(list_keys, list_vals, s=5)
    plt.plot(list_keys, list_vals)
    plt.ylabel("# of Electrons w/ that energy at any point")
    plt.xlabel("Energy (eV)")
    fig.savefig("NumElectronsvsEnergy.pdf")
    return

def plot_both_EDFs(MCC_eV_values, TCC_eV_values, ne, spe, ex, ey ,limit=15):
    # This function plots two different EEDF's onto one graph.
    fig = plt.figure()
    MCC_data = dict()
    for i in range(limit):
        MCC_data[float(i)] = 0
    for i in range(len(MCC_eV_values)):
        if not isnan(MCC_eV_values[i]):
            if floor(MCC_eV_values[i]) < limit:
                MCC_data[floor(MCC_eV_values[i])] += 1
    MCC_keys = []
    MCC_vals = []
    for key in MCC_data:
        if (key < limit) and (MCC_data[key] > 0):
            MCC_keys.append(key)
            MCC_vals.append(MCC_data[key])
    # Time for the theoretical
    TCC_data = dict()
    for i in range(limit):
        TCC_data[float(i)] = 0
    for i in range(len(TCC_eV_values)):
        if not isnan(TCC_eV_values[i]):
            if floor(TCC_eV_values[i]) < limit:
                TCC_data[floor(TCC_eV_values[i])] += 1
    TCC_keys = []
    TCC_vals = []
    for key in TCC_data:
        if (key < limit) and (TCC_data[key] > 0):
            TCC_keys.append(key)
            TCC_vals.append(TCC_data[key])

    # plt.xticks(np.arange(floor(min(list_all_eV_values)), floor(max(list_all_eV_values)) + 1, 1.0))
    # Plotting Measured Cross Sections
    plt.scatter(MCC_keys, MCC_vals, s=5, c="b")
    # Plotting Theoretical Cross Sections
    plt.scatter(TCC_keys, TCC_vals, s=5, c="r")
    # The below line might connect the scatter plot dots
    measured, = plt.plot(MCC_keys, MCC_vals, c="b", label = "measured")
    measured.set_label("Measured")
    theoretical, = plt.plot(TCC_keys, TCC_vals, c="r", label="theoretical")
    theoretical.set_label("Theoretical")
    plt.legend()
    # plt.title("EDF from Measured and Theoretical Cross Sections")
    plt.ylabel("# of Electrons w/ that energy at any point")
    plt.xlabel("Energy (eV)")
    fig.savefig("EDF_Ex_{0}_Ey_{1}_NE_{2}_SPE_{3}.pdf".format(ex, ey, ne, spe))
    return

def NumElecVsEnergySpreadsheet(all_eV_values):
    # First we check if the files already exists. If so we delete it. (If you do not delete it Python will continue
    # to write to an old file.
    silentremove('EnergyDistribution.xlsx')
    # Now we create the workbook
    workbook = xlsxwriter.Workbook('EnergyDistribution.xlsx')
    # Next we create the worksheet
    worksheet = workbook.add_worksheet()
    # CurrentRow is the row the script is currently writing into. Like most things in programming, we start at 0.
    CurrentRow = 0
    ener_bins = bin_eV_values_into_dict(all_eV_values)
    for key in ener_bins:
        # The variable below, InBetween, is just what it sounds like. To simplify graphing and spreadsheet writing,
        # the midpoint value of our bin is used to represent said bin. For example, 0.5 represents (0, 1.0),
        # 0.25 represents (0.0, 0.5), 0.15 represents (0.10, 0.20), and so on.
        InBetween = (key[0] + key[1])*0.5
        worksheet.write(CurrentRow, 0, InBetween)
        worksheet.write(CurrentRow, 1, ener_bins[key])
        CurrentRow += 1
    workbook.close()
    return

def sumExciteCSCFdict(allCS_CF_Dict):
    copyAllCS = dict(allCS_CF_Dict)
    SumExciteCS = dict()
    del copyAllCS["elastic"]
    del copyAllCS["ionization"]
    for currentEnergy in interpolated_ev:
        currentSum = 0
        for exciteCollisionCS in copyAllCS:
            currentSum += copyAllCS[exciteCollisionCS][currentEnergy]
        SumExciteCS[currentEnergy] = currentSum
    return SumExciteCS

def ElasIonSumExcite(AllCS):
    ElasIonSumAllExcite = dict()
    ElasIonSumAllExcite["elastic"] = AllCS["elastic"]
    ElasIonSumAllExcite["ionization"] = AllCS["ionization"]
    ElasIonSumAllExcite["sumExcite"] = sumExciteCSCFdict(AllCS)
    return ElasIonSumAllExcite

def CSvE2(iterALL, exciteLabel=False, SumExcite = True):
    # For simplicity purposes, going to home directory.
    goHomeDirectory()
    # Creating Figure
    fig = plt.figure()
    # iterEV's are all the interpolated electron volt values from 0.00 to 100.00 in intervals of 0.00.
    iterEV = return_iter_evs()
    ionLST = [iterALL["ionization"][elem] for elem in iterEV]
    elasLST = [iterALL["elastic"][elem] for elem in iterEV]
    if SumExcite:
        SumExciteDict = sumExciteCSCFdict(iterALL)
        SumExciteLST = [SumExciteDict[elem] for elem in iterEV]
    plt.plot(iterEV, elasLST, label="Elastic",linestyle="-")
    if type(exciteLabel) == type(str("string")):
        exciteLST = [iterALL[exciteLabel][elem] for elem in iterEV]
        plt.plot(iterEV, exciteLST, label=exciteLabel, linestyle=":")
    if SumExcite:
        plt.plot(iterEV, SumExciteLST, label="Sum of Excited States", linestyle=":")
    plt.plot(iterEV, ionLST, label="Ionization", linestyle="-.")
    plt.ylim(0, max(ionLST + elasLST + SumExciteLST))
    plt.xlim(0, 40)
    plt.ylabel("Cross Section $m^2$")
    plt.xlabel("Energy (eV)")
    lgd = plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    # plt.title("Cross Section vs Energy")
    goHomeDirectory()
    # This is the directory of all the results.
    results_directory = os.getcwd() + "/Results"
    # If it does not currently exist, the next two lines will create it using an if-statement to check if it exists.
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    os.chdir(results_directory)
    fig.savefig("Cross Section vs Energy.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    goHomeDirectory()
    return

def CFvE2(allCF = create_cf_dicts(), SumExcite = True, IndivExcite = False, Null = True):
    # Going back to home directory
    goHomeDirectory()
    # Changing directory to access Cross Section Files
    ionLST = [allCF["ionization"][elem] for elem in interpolated_ev]
    elasLST = [allCF["elastic"][elem] for elem in interpolated_ev]
    if type(IndivExcite) ==  type(str("string")):
        exciteLST = [exciteCfDict[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, exciteLST, label=IndivExcite, linestyle=":")
    # Creating Figure
    fig = plt.figure()
    plt.plot(interpolated_ev, elasLST, label="Elastic", linestyle="--")
    plt.plot(interpolated_ev, ionLST, label="Ionization", linestyle="-.")
    TOTLST = [TotalCF4ThisEnergy[elem] for elem in interpolated_ev]
    MaxY = max(TOTLST)
    plt.plot(interpolated_ev, TOTLST, label="TCF", linestyle="-", c="m")
    if Null:
        # The following 2 lines will plot the null collision frequency as a function of energy
        nullLST = [NCF4Energy[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, nullLST, label="Null", linestyle="-", c="black")
        # The following 4 lines will plot the total + null collision frequency as a function of energy
        TotWithNullLST = [(TotalCF4ThisEnergy[elem] + NCF4Energy[elem]) for elem in interpolated_ev]
        plt.plot(interpolated_ev, TotWithNullLST, label="Total w/ Null", linestyle="-", c="r")
        MaxY = max(TotWithNullLST)
    if SumExcite:
        SumExciteCfDict = sumExciteCSCFdict(allCF)
        SumExciteList = [SumExciteCfDict[elem] for elem in interpolated_ev]
        plt.plot(interpolated_ev, SumExciteList, label="SumExcite", linestyle="-", color="g")
    plt.ylim(0, MaxY * 1.025)
    plt.xlim(0, 40)
    plt.ylabel('Collision Frequency $s^{-1} $')
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
    goHomeDirectory()
    return

def computeSpecificCFs(All_CF_Dict, TotCF_Dict, writeSpreadSheet = True):
    goHomeDirectory()
    # This is the directory of all the results.
    results_directory = os.getcwd() + "/Results"
    # If it does not currently exist, the next two lines will create it using an if-statement to check if it exists.
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    # If an old SpecificCollisionFrequencies.xlsx file exists, we will remove it with the next
    try:
        os.remove("SpecificCollisionFrequencies.xlsx")
    except OSError:
        pass
    # We will have a total of six dictionaries, the line below only has 3
    SixCollisionFrequencies = ElasIonSumExcite(All_CF_Dict)
    # The next 3 will be Null, Total, and Tot + Null
    SixCollisionFrequencies["total"] = totalcfs_and_max2(AllIterCS)
    SixCollisionFrequencies["null"] = NCF4Energy
    # One more
    TotWithNullLST = [(TotalCF4ThisEnergy[elem] + NCF4Energy[elem]) for elem in interpolated_ev]
    SixCollisionFrequencies["tot+null"] = dict(zip(interpolated_ev, TotWithNullLST))

    # We are now defining a dictionary that will correlate CF-Type with an excel column
    CollisionColumn = dict()
    CollisionColumn["elastic"] = "B"
    CollisionColumn["sumExcite"] = "C"
    CollisionColumn["ionization"] = "D"
    CollisionColumn["total"] = "E"
    CollisionColumn["null"] = "F"
    CollisionColumn["tot+null"] = "G"

    # Below are the eV values of interest to us
    CheckValues = [0, 5, 10, 14, 15, 20, 25]
    # This will be the excel OrigValuesSheet
    workbook = xlsxwriter.Workbook("SpecificCollisionFrequencies.xlsx")
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

    # The outer layer of the for loop iterate through the eV values we desire to check. In other words, the rows.
    for value in CheckValues:
        # Calculating total CF for this energy
        totCF = TotCF_Dict[value]
        # Writing it to OrigValuesSheet
        ScaledValuesSheet.write("F"+ValueRow[value], totCF*scale)
        # The next layer of the for loop will iterate through all the collision types. In other words, the columns.
        for spec_cf in SixCollisionFrequencies:
            ScaledCF = SixCollisionFrequencies[spec_cf][value]*scale
            ScaledValuesSheet.write(str(CollisionColumn[spec_cf] + ValueRow[value]), float(format(ScaledCF, '.3f')))
    workbook.close()
    return


AllIterCS = returnIterCsDictionary()
# sumExciteCSCFdict(AllIterCS)
# CSvE2(AllIterCS)
# CFvE2()
# CFvE2(Null=False)
tcf_dict = totalcfs_and_max2(AllIterCS)
AllIterCF = create_cf_dicts()
computeSpecificCFs(AllIterCF, tcf_dict)