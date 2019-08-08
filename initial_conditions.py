
"""Now that all our functions and subroutines are defined, we can actually run electrons and conduct data analysis.
 The variables below will be all of our initial conditions."""

E_x = -25 # Electric Field x-component (in cartesian coordinate system)
E_y = 0 # Electric Field y-component (in cartesian coordinate system)
E_z = 0 # Electric Field z-component (in cartesian coordinate system)
NOE = 2 # Number of Electrons
number_of_electrons = 2
MPE = 3 # Motions Per Electron
motions_per_electron = 3
ground_state_atom_density = 1 * 10 ** 21 # Ground State Atom Density in atoms per meter cubed
metastable_state_atom_density = 1 * 10 ** 17 # MetaStable State Atom Density in atoms per meter cubed.
m_e = 9.1 * (10 ** (-31)) # Mass of an electron
electronCharge = 1.6 * (10 ** (-19)) * (-1) # Elementary Charge (negative)
argonMass = 6.6335209 * (10 ** (-26)) # Mass of an Argon atom
BinSize = 1 # Size of our bins
Joules2EV = 6.241509 * (10 ** 18) # Scale Factor: Multiply energy in Joules to convert to electron volts
EV2J = 1/Joules2EV # Scale Factor: Multiply energy in electron volts to convert to Joules
eedf_x_limit = 40
eedf_x_lim = 20
max_energy_ev = 100
energy_ev_cutoff = 50



# Below we define energy lost in electron volts from a collision
energy_loss_dictionary = dict()
energy_loss_dictionary["elastic"], energy_loss_dictionary["ionization"] = 0.0, 15.8 # Non-Excite Collisions with Argon
energy_loss_dictionary["null"] = 0.0
energy_loss_dictionary["gp1"] = 13.48 # Ground to 2nd excited state
energy_loss_dictionary["gp2"] = 13.33 # Ground to 2nd excited state
energy_loss_dictionary["gp3"] = 13.30 # Ground to 2nd excited state
energy_loss_dictionary["gp4"] = 13.28 # Ground to 2nd excited state
energy_loss_dictionary["gp5"] = 13.27 # Ground to 2nd excited state
energy_loss_dictionary["gp6"] = 13.17 # Ground to 2nd excited state
energy_loss_dictionary["gp7"] = 13.15 # Ground to 2nd excited state
energy_loss_dictionary["gp8"] = 13.09 # Ground to 2nd excited state
energy_loss_dictionary["gp9"] = 13.08 # Ground to 2nd excited state
energy_loss_dictionary["gp10"] = 12.91 # Ground to 2nd excited state
energy_loss_dictionary["s3p1"] = 1.76 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p2"] = 1.61 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p3"] = 1.58 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p4"] = 1.56 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p5"] = 1.55 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p6"] = 1.45 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p7"] = 1.43 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p8"] = 1.37 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p9"] = 1.36 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s3p10"] = 1.19 # 1st Excited (s3 metastable) to 2nd excited state
energy_loss_dictionary["s5p1"] = 1.93 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p2"] = 1.78 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p3"] = 1.75 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p4"] = 1.73 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p5"] = 1.72 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p6"] = 1.62 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p7"] = 1.60 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p8"] = 1.54 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p9"] = 1.53 # 1st Excited (s5 metastable) to 2nd excited state
energy_loss_dictionary["s5p10"] = 1.36 # 1st Excited (s5 metastable) to 2nd excited state
