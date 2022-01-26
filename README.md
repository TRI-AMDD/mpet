Development: [![Coverage Status](https://coveralls.io/repos/github/TRI-AMDD/mpet-dev/badge.svg?branch=development)](https://coveralls.io/github/TRI-AMDD/mpet-dev?branch=development)
Master: [![Coverage Status](https://coveralls.io/repos/github/TRI-AMDD/mpet-dev/badge.svg?branch=master)](https://coveralls.io/github/TRI-AMDD/mpet-dev?branch=master)
# MPET -- Multiphase Porous Electrode Theory

This software is designed to run simulations of batteries with porous electrodes using porous electrode theory, which is a volume-averaged, multiscale approach to capture the coupled behavior of electrolyte and active material within electrodes. As a result, with physical parameter inputs and run protocols (specified current or voltage profiles), it makes predictions about the internal dynamics within a battery (electrolyte concentration and potential, solid phase concentrations, reaction rates, etc.) and also macroscopic, easily measurable electrochemical quantities such as total current and voltage. In this way, it is similar to the [`dualfoil`](http://www.cchem.berkeley.edu/jsngrp/fortran.html) code released by Newman and coworkers from Berkeley. This software has much of the functionality contained in `dualfoil` (it is currently missing, e.g., temperature dependence). However, beyond the standard porous electrode theory simulations, this software can also simulate electrodes in which the active materials phase separate using non-equilibrium thermodynamics within a phase field modeling framework. Such behavior is common in widely used electrode materials, including graphite and LiFePO4.

If you use this software in academic work, please cite the relevant references detailing its development as presented in the `LICENSE` file. For more details on the theory implemeneted in the code, see:

Smith, R. B., and Bazant M. Z., Multiphase Porous Electrode Theory, [Journal of the Electrochemical Society](https://doi.org/10.1149/2.0171711jes), 2017, 164 (11) E3291-E3310, [arXiv preprint](https://arxiv.org/abs/1702.08432).

## Documentation

Documentation is available here ([https://mpet.readthedocs.io](https://mpet.readthedocs.io)) for installing, running, and analyzing results with mpet.


## Troubleshooting

Please use the Issues section of the Bitbucket repository (https://bitbucket.org/bazantgroup/mpet/issues) to file issues and/or bug reports with the software.
