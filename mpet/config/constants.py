# global constants

#: Reference temperature, K
T_ref = 298.
#: Boltzmann constant, J/(K particle)
k = 1.381e-23
#: Electron charge, C
e = 1.602e-19
#: Avogadro constant, particle / mol
N_A = 6.022e23
#: TODO: this parameter needs a name, C/mol
F = e * N_A
#: General particle classification (1 var)
two_var_types = ["diffn2", "CHR2", "homog2", "homog2_sdn"]
#: General particle classification (2 var)
one_var_types = ["ACR", "diffn", "CHR", "homog", "homog_sdn"]
