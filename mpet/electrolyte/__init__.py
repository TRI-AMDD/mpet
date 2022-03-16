r"""
This module provides functions defining properties of the ion-conducting
phase -- the electrolyte Manage functions for the parameters involved in
Stefan-Maxwell based concentrated electrolyte transport theory for
binary electrolytes.

Each electrolyte set must output functions for the following as a
function of c (electrolyte concentration, M)
 - Dchem [m^2/s] = the prefactor for grad(c) in species conservation
 - sigma [S/m] = the conductivity
 - (1 + dln(f_\pm)/dln(c)) = the "thermodynamic factor"
 - t_+^0 = the transference number of the cations
T in these equations is nondimensionalized wrt 298K
"""
