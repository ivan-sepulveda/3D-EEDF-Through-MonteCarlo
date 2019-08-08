import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

boffard_label = "[6] Boﬀard et. al."
hagelaar_label = "[5] Hagelaar et. al."

# --------------------------------------------------------------------------------------------------------------------
# 50td_10000electrons_comparison
# --------------------------------------------------------------------------------------------------------------------

pd.set_option('display.max_columns', 500)

ev_values_mc_50 = pd.read_excel("Results/Ey=0Ex=50Electrons=10000MotionsPer=250/50td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "T")['E (eV)'].astype(str)
ev_values_mc_50 = ev_values_mc_50.iloc[1:]
ev_values_mc_50 = ev_values_mc_50.astype(float)
ev_values_mc_50 = ev_values_mc_50.values

index_mc_50 = np.where(ev_values_mc_50==30.25)[0][0]
ev_values_mc_50 = ev_values_mc_50[:index_mc_50+1]

mc_normalized_50 = pd.read_excel("Results/Ey=0Ex=50Electrons=10000MotionsPer=250/50td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "W")
mc_normalized_50 = mc_normalized_50.dropna()
mc_normalized_50 = mc_normalized_50[1].values[1:]
mc_normalized_50 = mc_normalized_50[:index_mc_50+1]

ev_values_hagelaar_50 = pd.read_excel("Results/Ey=0Ex=50Electrons=10000MotionsPer=250/50td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "Z")
ev_values_hagelaar_50 = ev_values_hagelaar_50['E (eV).1'].astype(str)
ev_values_hagelaar_50 = ev_values_hagelaar_50.iloc[1:]
ev_values_hagelaar_50 = ev_values_hagelaar_50.astype(float)
ev_values_hagelaar_50 = ev_values_hagelaar_50.values

index_hagelaar_50 = np.where(ev_values_hagelaar_50==21.13)[0][0]
ev_values_hagelaar_50 = ev_values_hagelaar_50[:index_hagelaar_50+1]


hagelaar_normalized_50 = pd.read_excel("Results/Ey=0Ex=50Electrons=10000MotionsPer=250/50td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "AC")['Hagelaar'].iloc[1:]

hagelaar_normalized_50 = hagelaar_normalized_50.dropna()
hagelaar_normalized_50 = hagelaar_normalized_50.astype(float)
hagelaar_normalized_50 = hagelaar_normalized_50.values[:index_hagelaar_50+1]


ev_values_boffard_50 = pd.read_excel("Results/Ey=0Ex=50Electrons=10000MotionsPer=250/50td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "A").iloc[3:,0]

ev_values_boffard_50 = ev_values_boffard_50.astype(float).dropna()
ev_values_boffard_50 = ev_values_boffard_50.values


index_boffard_50 = len(ev_values_boffard_50) - 1

boffard_normalized_50 = pd.read_excel("Results/Ey=0Ex=50Electrons=10000MotionsPer=250/50td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "Q").iloc[3:,0]

boffard_normalized_50 = boffard_normalized_50.dropna()
boffard_normalized_50 = boffard_normalized_50.astype(float)
boffard_normalized_50 = boffard_normalized_50.values[:index_boffard_50+1]


results_comparison_50td = plt.figure()
# plt.title('E/V = 50 Td')
plt.plot(ev_values_mc_50, mc_normalized_50, label="Monte Carlo", color='darkred')
plt.plot(ev_values_hagelaar_50, hagelaar_normalized_50, label=hagelaar_label, color='darkblue')
plt.plot(ev_values_boffard_50, boffard_normalized_50, label=boffard_label, color='black')
# plt.xlim(0, min([max(ev_values_mc_50), max(ev_values_hagelaar_50), max(ev_values_hagelaar_50)]))
plt.xlim(0, 15)
plt.ylabel(r'EEDF (eV$^{-1}$)')
plt.xlabel('Energy (eV)')
plt.legend()
plt.show()
results_comparison_50td.savefig("Ex=50 Motions=10k")
plt.close('all')


# --------------------------------------------------------------------------------------------------------------------
# 25td_10000electrons_comparison
# --------------------------------------------------------------------------------------------------------------------


pd.set_option('display.max_columns', 500)

ev_values_mc_25 = pd.read_excel("Results/Ey=0Ex=25Electrons=10000MotionsPer=250/25td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "T")['E (eV)'].astype(str)
#
ev_values_mc_25 = ev_values_mc_25.iloc[1:]
ev_values_mc_25 = ev_values_mc_25.astype(float)
ev_values_mc_25 = ev_values_mc_25.values

index_mc_25 = np.where(ev_values_mc_25==30.25)[0][0]
ev_values_mc_25 = ev_values_mc_25[:index_mc_25+1]

mc_normalized_25 = pd.read_excel("Results/Ey=0Ex=25Electrons=10000MotionsPer=250/25td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "W")
mc_normalized_25 = mc_normalized_25.dropna()
mc_normalized_25 = mc_normalized_25[1].values[1:]
mc_normalized_25 = mc_normalized_25[:index_mc_25+1]

ev_values_hagelaar_25 = pd.read_excel("Results/Ey=0Ex=25Electrons=10000MotionsPer=250/25td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "Z")
ev_values_hagelaar_25 = ev_values_hagelaar_25['E (eV).1'].astype(str)
ev_values_hagelaar_25 = ev_values_hagelaar_25.iloc[1:]
ev_values_hagelaar_25 = ev_values_hagelaar_25.astype(float)
ev_values_hagelaar_25 = ev_values_hagelaar_25.values

index_hagelaar_25 = np.where(ev_values_hagelaar_25==21.13)[0][0]
ev_values_hagelaar_25 = ev_values_hagelaar_25[:index_hagelaar_25+1]


hagelaar_normalized_25 = pd.read_excel("Results/Ey=0Ex=25Electrons=10000MotionsPer=250/25td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "AC")['Hagelaar'].iloc[1:]

hagelaar_normalized_25 = hagelaar_normalized_25.dropna()
hagelaar_normalized_25 = hagelaar_normalized_25.astype(float)
hagelaar_normalized_25 = hagelaar_normalized_25.values[:index_hagelaar_25+1]


ev_values_boffard_25 = pd.read_excel("Results/Ey=0Ex=25Electrons=10000MotionsPer=250/25td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "A").iloc[3:,0]

ev_values_boffard_25 = ev_values_boffard_25.astype(float).dropna()
ev_values_boffard_25 = ev_values_boffard_25.values

index_boffard_25 = len(ev_values_boffard_25) - 1

boffard_normalized_25 = pd.read_excel("Results/Ey=0Ex=25Electrons=10000MotionsPer=250/25td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "Q").iloc[3:,0]

boffard_normalized_25 = boffard_normalized_25.dropna()
boffard_normalized_25 = boffard_normalized_25.astype(float)
boffard_normalized_25 = boffard_normalized_25.values[:index_boffard_25+1]


results_comparison_25td = plt.figure()
# plt.title('E/V = 25 Td')
plt.plot(ev_values_mc_25, mc_normalized_25, label="Monte Carlo", color='darkred')
plt.plot(ev_values_hagelaar_25, hagelaar_normalized_25, label=hagelaar_label, color='darkblue')
plt.plot(ev_values_boffard_25, boffard_normalized_25, label=boffard_label, color='black')
# plt.xlim(0, min([max(ev_values_mc_25), max(ev_values_hagelaar_25), max(ev_values_hagelaar_25)]))
plt.xlim(0, 15)
plt.ylabel(r'EEDF (eV$^{-1}$)')
plt.xlabel('Energy (eV)')
plt.legend()
plt.show()
results_comparison_25td.savefig("Ex=25 Motions=10k")
plt.close('all')




pd.set_option('display.max_columns', 500)

ev_values_mc_10 = pd.read_excel("Results/Ey=0Ex=10Electrons=10000MotionsPer=250/10td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "T")['E (eV)'].astype(str)
#
ev_values_mc_10 = ev_values_mc_10.iloc[1:]
ev_values_mc_10 = ev_values_mc_10.astype(float)
ev_values_mc_10 = ev_values_mc_10.values

index_mc_10 = np.where(ev_values_mc_10==30.25)[0][0]
ev_values_mc_10 = ev_values_mc_10[:index_mc_10+1]

mc_normalized_10 = pd.read_excel("Results/Ey=0Ex=10Electrons=10000MotionsPer=250/10td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "W")
mc_normalized_10 = mc_normalized_10.dropna()
mc_normalized_10 = mc_normalized_10[1].values[1:]
mc_normalized_10 = mc_normalized_10[:index_mc_10+1]

ev_values_hagelaar_10 = pd.read_excel("Results/Ey=0Ex=10Electrons=10000MotionsPer=250/10td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "Z")
ev_values_hagelaar_10 = ev_values_hagelaar_10['E (eV).1'].astype(str)
ev_values_hagelaar_10 = ev_values_hagelaar_10.iloc[1:]
ev_values_hagelaar_10 = ev_values_hagelaar_10.astype(float)
ev_values_hagelaar_10 = ev_values_hagelaar_10.values

index_hagelaar_10 = np.where(ev_values_hagelaar_10==17.12)[0][0]
ev_values_hagelaar_10 = ev_values_hagelaar_10[:index_hagelaar_10+1]


hagelaar_normalized_10 = pd.read_excel("Results/Ey=0Ex=10Electrons=10000MotionsPer=250/10td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "AC")['Hagelaar'].iloc[1:]

hagelaar_normalized_10 = hagelaar_normalized_10.dropna()
hagelaar_normalized_10 = hagelaar_normalized_10.astype(float)
hagelaar_normalized_10 = hagelaar_normalized_10.values[:index_hagelaar_10+1]


ev_values_boffard_10 = pd.read_excel("Results/Ey=0Ex=10Electrons=10000MotionsPer=250/10td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "A").iloc[3:,0]

ev_values_boffard_10 = ev_values_boffard_10.astype(float).dropna()
ev_values_boffard_10 = ev_values_boffard_10.values

index_boffard_10 = len(ev_values_boffard_10) - 1

boffard_normalized_10 = pd.read_excel("Results/Ey=0Ex=10Electrons=10000MotionsPer=250/10td_10000electrons_comparison.xlsx",
                   sheet_name='Sheet1', usecols = "Q").iloc[3:,0]

boffard_normalized_10 = boffard_normalized_10.dropna()
boffard_normalized_10 = boffard_normalized_10.astype(float)
boffard_normalized_10 = boffard_normalized_10.values[:index_boffard_10+1]


results_comparison_10td = plt.figure()
# plt.title('E/V = 10 Td')
plt.plot(ev_values_mc_10, mc_normalized_10, label="Monte Carlo", color='darkred')
plt.plot(ev_values_hagelaar_10, hagelaar_normalized_10, label=hagelaar_label, color='darkblue')
plt.plot(ev_values_boffard_10, boffard_normalized_10, label=boffard_label, color='black')
# plt.xlim(0, min([max(ev_values_mc_10), max(ev_values_hagelaar_10), max(ev_values_hagelaar_10)]))
plt.xlim(0, 15)
plt.ylabel(r'EEDF (eV$^{-1}$)')
plt.xlabel('Energy (eV)')
plt.legend()
plt.show()
results_comparison_10td.savefig("Ex=10 Motions=10k")
plt.close('all')