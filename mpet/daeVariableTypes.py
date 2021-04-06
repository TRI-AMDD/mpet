"""This module defines custom daetools variable types"""
import daetools.pyDAE as dae

mole_frac_t = dae.daeVariableType(
    name="mole_frac_t", units=dae.unit(), lowerBound=0,
    upperBound=1, initialGuess=.25, absTolerance=1.e-6)
conc_t = dae.daeVariableType(
    name="conc_t", units=dae.unit(), lowerBound=0,
    upperBound=1e20, initialGuess=1.00, absTolerance=1.e-6)
elec_pot_t = dae.daeVariableType(
    name="elec_pot_t", units=dae.unit(), lowerBound=-1e20,
    upperBound=1e20, initialGuess=0, absTolerance=1.e-6)
