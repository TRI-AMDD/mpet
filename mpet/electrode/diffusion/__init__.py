"""This folder contains diffusion functions that return
 the filling-fraction dependent variation of
the transport coefficient, D(y), such that
Flux = -D_ref*D(y)*grad(y) for solid solution transport or
Flux = -D_ref*D(y)*grad(mu) for thermo-based transport
where y here is the filling fraction, D_ref has dimensions of
length^2/time, D(y) is dimensionless, and mu, the chemical
potential, has been scaled to the thermal energy, k*Tref. For more
details on the two forms, see muRfuncs which defines both solid
solution (_ss) and materials based on simpler thermodynamic models.
"""
