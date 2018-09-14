
"""Now that all our functions and subroutines are defined, we can actually run electrons and conduct data analysis.
 The variables below will the all of our initial conditions.

E_x, E_y, and E_z: electric fields applied in Volts/Meter.
NOE: number of electrons
MPE: motions per electron
GSAD: Ground State Atom Density in atoms per meter cubed.
MSSAD: MetaStable State Atom Density in atoms per meter cubed.

"""

E_x = -25
E_y = 0
E_z = 0
NOE = 50
MPE = 50
GSAD = 1*10**21
MSSAD = 1*10**17
m_e = 9.1 * (10 ** (-31))
electronCharge = 1.6 * (10 ** (-19)) * (-1)
argonMass = float(6.6335209 * (10 ** (-26)))
BinSize = 0.5
Joules2EV = 6.241509 * (10 ** 18)
EV2J = 1/Joules2EV
eedf_x_limit = 40



# Below we define energy lost in electron volts from a collision
e_loss = dict()
e_loss["elastic"] = 0.0
e_loss["ionization"] = 15.8
e_loss["gp1"] = 13.48
e_loss["gp2"] = 13.33
e_loss["gp3"] = 13.30
e_loss["gp4"] = 13.28
e_loss["gp5"] = 13.27
e_loss["gp6"] = 13.17
e_loss["gp7"] = 13.15
e_loss["gp8"] = 13.09
e_loss["gp9"] = 13.08
e_loss["gp10"] = 12.91
e_loss["s3p1"] = 1.76
e_loss["s3p2"] = 1.61
e_loss["s3p3"] = 1.58
e_loss["s3p4"] = 1.56
e_loss["s3p5"] = 1.55
e_loss["s3p6"] = 1.45
e_loss["s3p7"] = 1.43
e_loss["s3p8"] = 1.37
e_loss["s3p9"] = 1.36
e_loss["s3p10"] = 1.19
e_loss["s5p1"] = 1.93
e_loss["s5p2"] = 1.78
e_loss["s5p3"] = 1.75
e_loss["s5p4"] = 1.73
e_loss["s5p5"] = 1.72
e_loss["s5p6"] = 1.62
e_loss["s5p7"] = 1.60
e_loss["s5p8"] = 1.54
e_loss["s5p9"] = 1.53
e_loss["s5p10"] = 1.36
e_loss["null"] = 0.0
