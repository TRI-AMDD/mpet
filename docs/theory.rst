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

For more details on the theory implemented in the code, see: ::

    Smith, R. B., and Bazant M. Z.,
    Multiphase Porous Electrode Theory,
    Journal of the Electrochemical Society, 2017, 164 (11) E3291-E3310
    https://iopscience.iop.org/article/10.1149/2.0171711jes



Model
------------------------------

Multiscale Porous Electrode Theory (MPET) is a multiscale model whose original form is porous electrode theory, but in a thermodynamically consistent framework to be able to study phase separating materials. It is able to capture particle scale phase separation, but also reflect the results in a porous electrode scale framework in an efficient manner to to computationally efficient. We include a simple description of MPET here, but more can be read in Smith & Bazant, J. Electrochem. Soc., 2017; and Bazant, Acc. Chem. Res., 2013. 

There are two scales in such a model to efficiently simulate a porous electrode. The particle scale model consists of multiple particles with the same volume scale properties, e.g. in the same electrolyte volume bath at the chemical potential :math:`{\phi}_l` and with the same solid potential :math:`{\phi}_s`. The solid scale reaction equation can be modeled using a simple diffusion model for a non-phase separating model, an Allen-Cahn equation for a depth-averaged phase separating model, and the Cahn-Hilliard equation for a full diffusion-reaction phase separating model.

The Cahn-Hilliard and diffusion equations are modeled as

.. math::

    \frac{\partial c_i}{\partial t} = - \nabla \cdot \mathbf{F}_i

where :math:`\mathbf{F}_i = - \frac{D_i c_i}{k_B T} \nabla {\mu}_i` is the flux from the chemical potential (:math:`{\mu}` gradient). Reaction boundary conditions of :math:`\mathbf{n} \cdot \mathbf{F}_i = -R_{p,i}` are applied for non phase separating models, while natural boundary condition from the surface free energy are applied for phase separating models, with the surface energy :math:`{\gamma}_s` as :math:`\mathbf{n} \cdot \frac{\partial g}{\partial \nabla c_i} = \mathbf{n} \cdot {\kappa} \nabla c_i = \frac{\partial {\gamma}_s}{\partial c_i}`.

The difference between the reaction models is that the Allen-Cahn equation and homogenuous reactions both assume homogenuous reaction, compared to boundary reactions as in the diffusion/Cahn-Hilliard model. The Allen-Cahn/depth averaged model has only been verified to be an accurate model for lithium iron phosphate, but the Allen-Cahn model and the homogenuous model can both be written as

.. math::

    \frac{\partial c_i}{\partial t} = \frac{A_p}{V_p} j_p

where :math:`A_p` is the area of the particle, :math:`V_p` is the volume of the particle, and :math:`j_p` is the reaction flux without accounting for diffusive transport.

Multiple reaction rates can be used in the single particle models, including Butler Volmer reactions with different prefactors, Marcus reactions, Marcus-Hush-Chisdey reactions, and coupled ion electron transfer, as long as the reactions are thermodynamically consistent (satisfy De Donder's equation). The most commonly used equation in electrochemical modeling, the Butler Volmer reaction, is modeled as

.. math::

    i = i_0 \left( \exp{\left( -\frac{{\alpha}e{\eta}_{eff}}{k_B T}\right)} - \exp{\left( \frac{\left( 1-{\alpha}\right)e{\eta}_{eff}}{k_B T}\right)}\right)

where :math:`{\eta}_{eff}` is the effective overpotential that accounts for film resistance, :math:`{\eta}_{eff} = {\eta} + iR_{film}` and :math:`i_0` is the exchange current density. In the BV equation and other models, :math:`{alpha}` is the charge transfer coefficient, which is 0.5 for symmetric reactions. 
In contrast with the traditional BV equations, there are reaction models that assume that electron transfer is the dominant reaction mechanism, all Marcus-type models. Marcus and Marcus-Hush-Chidsey (MHC) kinetics are similar models, but one assumes a localized electron density and the other assumes a delocalized electron density, corresponding to insulating and metallic states respectively, where only the density of states :math:`\rho(z)` integral differs between the two. In the Marcus model, the density of states is a delta function, while for the MHC model it is delocalized. The reduction and oxidation currents can be modeled as 

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

is true as a simplification to the MHC model, to prevent having to evaluate the integral over the density of states, a computationally expensive procedure.

The coupled-ion electron transfer model was developed by our group as a way of modeling the coupling of electron and lithium ion transfer in a reaction. In this model, electron transfer is modeled quantum mechanically but ion transfer is modeled classically. The model can be expressed as

.. math::
    i_{red} = \frac{e k^*_0 c_O}{{\gamma}_{TS}} n_e \exp{\left( - \frac{\left( {\lambda}+e{\eta}_f-x_{\varepsilon}\right)^2}{4 {\lambda}k_B T}\right)} \\
    i_{ox} = \frac{e k^*_0 c_R}{{\gamma}_{TS}} (1-n_e) \exp{\left( - \frac{\left( {\lambda}-e{\eta}_f-x_{\varepsilon}\right)^2}{4 {\lambda}k_B T}\right)} \\
    i = \int_{-\infty}^{\infty} \left( i_{red} - i_{ox}\right){\rho} d{\varepsilon}

where :math:`x_{\varepsilon} = E_f - {\varepsilon}` and :math:`e{\eta_f} = e{\eta} - k_B T \ln{c_O/c_R}` and :math:`k^*_0` is the reaction state prefactor, :math:`{\gamma}_{TS}` is the transition state overpotential.

Now we return to our volume scale model to reflect the behavior of the porous electrode. System scale transport can be modeled using transport equations through multiple volumes. The electrolyte transport equations can be modeled

.. math::
    \frac{\partial \epsilon c_{l,i}}{\partial t} = - \nabla \cdot \mathbf{F_{l,i}} + R_{V,i}

where :math:`R_{V,i}` is the total volumetric reaction rate in the volume modeled.

For a charge balance, we assume that quasineutrality is also satisified over the length scales, giving

.. math::
    \frac{\partial \epsilon \rho_e }{\partial t} \approx 0 = - \nabla \cdot \mathbf{i_l} + \sum_i z_i eR_{V,i}

for the current equations. No flux boundary conditions for concentration and current are applied at the current collector as :math:`\mathbf{n} \cdot \mathbf{F_l} = 0` and :math:`\mathbf{n} \cdot \mathbf{i_l} = 0`.

For electrolyte transport, dilute or concentrated (Stefan-Maxwell) electrolyte models can be used for electrolyte transport in MPET. Since we defined the chemical potential as :math:`{\mu}_{l,i} = k_B T \ln{a_{l,i}} + {\mu}_{l,i}^0 + z_i e {\phi}_l`, where :math:`{\mu}_{l,i}^0` is the reference state chemical potential and :math:`{\phi}_l` is the electrostatic potential, in a dilute solution model transport can simply be modeled as 

.. math::
    \mathbf{F}_{l,i} = -\left( D_{l,chem,i} \nabla c_{l,i} + \frac{D_{l,i}c_{l,i}z_i}{k_B T} \nabla {\phi}_l\right)

with :math:`D_{l,chem,i} = D_{l,i}\left( 1 + \frac{\partial \ln{{\gamma}_{l,i}}}{\partial \ln{c_{l,i}}}\right)` and :math:`{\gamma}_{l,i}` as the activity coefficient of species :math:`i`. 

Meanwhile, the Stefan-Maxwell concentrated solution model includes larger concentration gradients, which requires that gradients of species :math:`i` also affect the transport of species :math:`j` and thus is more commonly used in battery modeling because of its higher accuracy. Newman simplified the model for a binary electrolyte in his classic textbook to be


.. math::
    \mathbf{F}_{l,+} = - \frac{{\nu}_{+} {\epsilon}}{{\tau} }D_l \nabla c_l + \frac{t^0_{+} \mathbf{i}_l}{z_{+} e} \\
    \mathbf{F}_{l,-} = - \frac{{\nu}_{-} {\epsilon}}{{\tau} }D_l \nabla c_l + \frac{t^0_{-} \mathbf{i}_l}{z_{-} e}


where 


.. math::
    D_l = \mathbf{D} \frac{c_T}{c_{l,0}}\left( 1 + \frac{\partial \ln{{\gamma}_{l,\pm}}}{\partial \ln{c_l}}\right) \\
    \mathbf{D} = \frac{\mathbf{D}_{0+} \mathbf{D}_{0-} (z_+ - z_{-})}{z_+\mathbf{D}_{0+} - z_{-}\mathbf{D}_{0-}} \\
    t_{+}^0 = 1 - t_{-}^0 = \frac{z_+ \mathbf{D}_{0+}}{z_+ \mathbf{D}_{0+} - z_{-} \mathbf{D}_{0-}}
 

The effective diffusivity of the cation/anion in these models is defined as :math:`c_{l,i} \nabla {\mu}_{l,i} = k_B T \sum_j \frac{c_{l,i}c_{l,j}}{c_T \mathbf{D}_{ij}}\left( \mathbf{v}_j - \mathbf{v}_i\right)`. 

In addition to electrolyte transport, the solid microstructure model for charge is also defined using conservation of charge, with :math:`0 = -\nabla \cdot \mathbf{i}_2 - \sum_i z_i e R_{V,i}`, where the current density is assumed to be modeled with an Ohm's law equation :math:`\mathbf{i}_s = - \frac{1-{\epsilon}}{{\tau}} {\sigma}_s \nabla {\phi}_s`, where :math:`{\sigma}_s` is the solid conductivity of the system.

The resistance between particles can be modeled with parallel contact resistances as :math:`G_{j,k}\left( {\phi}_{j,k}-{\phi}_{j,k+1} = I_{j,k}\right)`, where :math:`G_{j,k}` is the conductance and :math:`I_{j,k}` is the current between particles. Charge conservation requires that :math:`I_{j,k} - I_{j,k+1} = \int_{A_{k+1}} j_{j,k+1}dA`, where :math:`j_{j,k+1}` is the intercalation rate into particle k+1 with surface area :math:`A_{k+1}`. 

.. The coupling between the particle scale and electrode model is achieved through "replicating" the set of particles we are simulating based on the amount of active material in each volume of the battery electrode, and through the current or voltage constraints throughout the battery. The total volumetric reaction term is defined as :math:`R_{V,i} = - \left( 1-{\epsilon}\right)P_L \sum_p \frac{V_p}{V_u} \frac{\partial \bar{c}_{p,i}}{\partial t}`, where :math:`V_u = \sum_p V_p` is the sum over particle areas and :math:`P_L` is the loading of active material in the solid phase.
.. The current constraint or voltage constraint equations can be written as 

.. math::
    i_{cell} = \sum_i \int_{L_a} z_i e R_{V,i} dL = - \sum_i \int_{L_c} z_i e R_{V,i} dL

over either the anode (left side) or cathode (right side).
The overall cell voltage is defined as 

.. math::
    \Delta {\phi}_{cell} = \Delta {\phi}_{appl} - i_{cell} R_{ser}

where :math:`R_{ser}` is the resistance of the cell per area, which closes our porous electrode model.
