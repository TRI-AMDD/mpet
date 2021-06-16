Theory
===========================================


This software is designed to run simulations of batteries with porous electrodes using porous electrode theory,
which is a volume-averaged, multiscale approach to capture the coupled behavior of electrolyte and active material
within electrodes. As a result, with physical parameter inputs and run protocols (specified current or voltage
profiles), it makes predictions about the internal dynamics within a battery (electrolyte concentration and potential,
solid phase concentrations, reaction rates, etc.) and also macroscopic, easily measurable electrochemical quantities
such as total current and voltage. In this way, it is similar to the dualfoil code released by Newman and coworkers
from Berkeley. This software has much of the functionality contained in dualfoil (it is currently missing, e.g.,
temperature dependence). However, beyond the standard porous electrode theory simulations, this software can also
simulate electrodes in which the active materials phase separate using non-equilibrium thermodynamics within a phase
field modeling framework. Such behavior is common in widely used electrode materials, including graphite and LiFePO4.

If you use this software in academic work, please cite the relevant references detailing
its development as presented in the LICENSE file.

Reference
------------------------------------------------------------

For more details on the theory implemeneted in the code, see: ::

    Smith, R. B., and Bazant M. Z.,
    Multiphase Porous Electrode Theory,
    Journal of the Electrochemical Society, 2017, 164 (11) E3291-E3310
    https://iopscience.iop.org/article/10.1149/2.0171711jes



Model
------------------------------

[To Do] In that section we could quickly introduce the theoretical background. Equation can easily be displayed :

.. math::

    \frac{\partial (\epsilon c_{l,i})}{\partial t} = -\nabla \cdot \mathbf{F}_{l,i} + R_{V,i}

Multiscale Porous Electrode Theory (MPET) is a multiscale model based off porous electrode theory, but in a thermodynamically consistent framework to be able to study phase separating materials. It is able to capture particle scale phase separation, but also reflect the results in a porous electrode scale framework.

The particle scale model consists of multiple particles in the same electrolyte volume bath at the chemical potential :math:`{\phi}_l`. The solid scale reaction equation can be modeled using a simple diffusion model for a non-phase separating model, an Allen-Cahn equation for a depth-averaged phase separating model, and the Cahn-Hilliard equation for a full diffusion-reaction phase separating model.

The Cahn-Hilliard and diffusion equations are modeled as

.. math::

    \frac{\partial c_i}{\partial t} = - \nabla \cdot \mathbf{F}_i

where :math:`\mathbf{F}_i = - \frac{D_i c_i}{k_B T} \nabla {\mu}_i` is the flux. Reaction boundary conditions of :math:`\mathbf{n} \cdot \mathbf{F}_i = -R_{p,i}` are applied for non phase separating models, while natural boundary condition from the surface free energy are applied for phase separating models, with the surface energy :math:`{\gamma}_s` as :math:`\mathbf{n} \cdot \frac{\partial g}{\partial \nabla c_i} = \mathbf{n} \cdot {\kappa} \nabla c_i = \frac{\partial {\gamma}_s}{\partial c_i}`.

The Allen-Cahn equation and homogenuous reactions both assume homogenuous reaction. The Allen-Cahn/depth averaged model has only been verified for lithium iron phosphate models, but both can be written as

.. math::

    \frac{\partial c_i}{\partial t} = \frac{A_p}{V_p} j_p

where :math:`A_p` is the area of the particle, :math:`V_p` is the volume of the particle, and :math:`j_p` is the reaction flux without accounting for diffusive transport.

Multiple reaction rates can be used in the single particle models, including Butler Volmer reactions with different prefactors, Marcus reactions, Marcus-Hush-Chisdey reactions, and coupled ion electron transfer. The Butler Volmer reaction is modeled as

.. math::

    i = i_0 \left( \exp{\left( -\frac{{\alpha}e{\eta}_{eff}}{k_B T}\right)} - \exp{\left( \frac{\left( 1-{\alpha}\right)e{\eta}_{eff}}{k_B T}\right)}\right)

where :math:`{\eta}_{eff}` is the effective overpotential that accounts for film resistance, :math:`{\eta}_{eff} = {\eta} + iR_{film}` and :math:`i_0` is the exchange current density.
Marcus and Marcus-Hush-Chidsey (MHC) kinetics are similar models, but one assumes a localized electron density and the other assumes a delocalized electron density, corresponding to insulating and metallic states respectively, where only the density of states :math:`\rho(z)` integral differs between the two. The reduction and oxidation currents can be used as 

.. math::

    i_{red} = ek^0c_0\int_{-\infty}^{\infty} {\rho}(z)n_e(z)W_{red}dz \\
    i_{ox} = ek^0c_R\int_{-\infty}^{\infty} {\rho}(z)\left(1-n_e(z)\right)W_{ox}dz

where :math:`k^0` is a prefactor, :math:`z` is the energy level of the electrons, and :math:`n_e(z)` is the Fermi function. The transition probabilities :math:`W` are found to be 

.. math::

    W_{red} = k_w \exp{\left( - \frac{w_O}{k_B T}\right)} \exp{\left( - \frac{\left( {\lambda} + e{\eta}_f\right)^2}{4{\lambda} k_B T}\right)} \\
    W_{ox} = k_w \exp{\left( - \frac{w_R}{k_B T}\right)} \exp{\left( - \frac{\left( {\lambda} - e{\eta}_f\right)^2}{4{\lambda} k_B T}\right)}

where :math:`{\lambda}` is the reorganization energy. The driving force of the system, the modified overpotential, is defined as :math:`e{\eta}_f = e{\eta} + k_B T \ln{\left( \frac{c_O}{c_R}\right)}`.

Meanwhile, for the Marcus Hush Chidsey model, the reaction model can be written as 

.. math::

    i_{red/ox} = k_M c_{O/R}\exp{\left( -\frac{w_{O/R}}{k_B T}\right)} \int_{-\infty}^{\infty} \exp{\left( - \frac{\left( z-{\lambda}\mp e{\eta}_f\right)^2}{4 {\lambda} k_B T}\right)} \frac{dz}{1+\exp{(z/k_BT)}}

Following Zeng et al. J. Electroanal. Chem., 2014, we see that 

.. math::
    k_{red/ox} \approx \frac{\sqrt{{\pi}\tilde{{\lambda}}}}{1+\exp{\left( \pm \tilde{{\eta}}_f\right)}} erfc \left( \frac{\tilde{{\lambda}} - \sqrt{1+ \sqrt{\tilde{{\lambda}}}+\tilde{{\eta}}_f^2}}{2\sqrt{\tilde{{\lambda}}}}\right)

is true as a simplification to the MHC model.
