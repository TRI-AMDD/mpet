UIhelp='''********************************** MPET user interface *********************************
mpetUI provides a user interface to setup and run MPET simulations.

COMMANDS:

	****** To quit mpetUI ******
	$mpetUI --quit 

	****** To display help for mpetUI ******
	$mpetUI --help 

	****** To run the user interface in GUI mode ******
	$mpetUI --gui 

	****** To run the user interface in text mode ******
	The text mode allows UI to input required simulation parameters as arguments.
	$mpetUI --paramater1=value1
	$mpetUI --paramater1=value1 --parameter2=value2 ..

PARAMETERS:

  Description of parameters for mpetUI.

  ******** Write or replace the parameter configuration file ***********

	--Anode_config_file   : Specify the Anode parameter configuration file
	--Cathode_config_file : Specify the Cathode parameter configuration file
	--System_config_file  : Specify the System parameter configuration file

  ******** Anode Particle Parameters ********

	--anode_particle 				(legacy mode: "type") [-]
	  Type of solid particles(ACR,CHR,CHR2,homog,homog_sdn,homog2,diffn) for anode

	--anode_particle_discretization (legacy mode: "discretization") [m]
	  Discretization of solid, Used for non-homogeneous type

	--anode_particle_shape 			(legacy mode: "shape") [-]
	  Shape of the particle (sphere,cylinder,c3)

	--anode_particle_thickness      (legacy mode: "thickness") [m]
	  C3 a-axis length or cylinder particle thickness,

  ******** Anode Material Parameters ******** 

	--anode_material_muRfunc 		(legacy mode: "muRfunc") [-]
	  The chemical potential of the reduced state as a function of filling fraction.

	--anode_material_Logpad 		(legacy mode: "logpad") [-]
	  Log pad (option for muRfuncs containing log terms -- not solid solution)

	--anode_material_noise 			(legacy mode: "noise") [-]
	  Add Langevin noise to the solid dynamics to simulate some random thermal fluctuations

	--anode_material_mag_noise 		(legacy mode: "noise_prefac")  [-]
	  The standard deviation of the normal noise added to the solid dynamics

	--anode_material_num_noise 		(legacy mode: "numnoise")  [-]
	  Number of segments over which the noise is re-randomized

	--anode_material_oma 			(legacy mode: "Omega_a") [J/Li]
	  Regular solution parameter

	--anode_material_omb 			(legacy mode: "Omega_b") [J/Li]
	  Secondary interaction terms for 2-variable models

	--anode_material_omc 			(legacy mode: "Omega_c") [J/Li] 
	  Secondary interaction terms for 2-variable models

	--anode_material_grad_pen 		(legacy mode: "kappa")  [J/m]
	  Gradient penalty,

	--anode_material_stress_coeff 	(legacy mode: "B")  [Pa]
	  Stress coefficient,  

	--anode_material_EvdW 			(legacy mode: "EvdW")  [J/Li]
	  Van der Waals-like interaction energy for graphite models

	--anode_material_rho_s 			(legacy mode: "rho_s")  [Li/m^3]
	  Total Li site density within the solid

	--anode_material_D 				(legacy mode: "D") [m^2/s]
	  Dimensional diffusivity prefactor (used if type = diffn[2], CHR[2])

	--anode_material_Dfunc 			(legacy mode: "Dfunc")  [-]
	  Dimensionless diffusivity function,

	--anode_material_dgammadc 		(legacy mode: "dgammadc")  [J*m/Li]
	  Change in surface energy with respect to solid concentration,

	--anode_material_cwet 			(legacy mode: "cwet")  [-]
	  ACR dimensionless surface wetting (used with LiFePO4, type=ACR)

  ******** Anode Reaction Parameters ********

	--anode_reaction_type 			(legacy mode: "rxnType") [-]
	  Reaction type (BV,BV_gMod01,BV_raw,BV_mod01,BV_mod02,Marcus,MHC)

	--anode_reaction_ecd 			(legacy mode: "k0")  [A/m^2]
	  Exchange current density, (0.1 for h2/Pt)

	--anode_reaction_charget 		(legacy mode: "alpha")  [-]
	  Charge transfer coefficient (only used for BV)

	--anode_reaction_Rorg_energy 	(legacy mode: "lambda")  [J/Li]
	  Reorganizational energy, (used for Marcus and MHC)

	--anode_reaction_filmres 		(legacy mode: "Rfilm")  [Ohm m^2]
	  Film resistance,


  ******** Cathode Particle Parameters ********

	--cathode_particle 				(legacy mode: "type") [-]
	  Type of solid particles(ACR,CHR,CHR2,homog,homog_sdn,homog2,diffn) for cathode

	--cathode_particle_discretization (legacy mode: "discretization") [m]
	  Discretization of solid,  Used for non-homogeneous type

	--cathode_particle_shape 			(legacy mode: "shape") [m]
	  Shape of the particle (sphere,cylinder,c3)

	--cathode_particle_thickness      (legacy mode: "thickness") [m] 
	  C3 a-axis length or cylinder particle thickness,

  ******** Cathode Material Parameters ******** 

	--cathode_material_muRfunc 		(legacy mode: "muRfunc") [-]
	  The chemical potential of the reduced state as a function of filling fraction.

	--cathode_material_Logpad 		(legacy mode: "logpad") [-]
	  Log pad (option for muRfuncs containing log terms -- not solid solution)

	--cathode_material_noise 			(legacy mode: "noise") [-]
	  Add Langevin noise to the solid dynamics to simulate some random thermal fluctuations

	--cathode_material_mag_noise 		(legacy mode: "noise_prefac")  [-]
	  The standard deviation of the normal noise added to the solid dynamics

	--cathode_material_num_noise 		(legacy mode: "numnoise")  [-]
	  Number of segments over which the noise is re-randomized

	--cathode_material_oma 			(legacy mode: "Omega_a") [J/Li]
	  Regular solution parameter,

	--cathode_material_omb 			(legacy mode: "Omega_b")  [J/Li]
	  Secondary interaction terms for 2-variable models

	--cathode_material_omc 			(legacy mode: "Omega_c")  [J/Li]
	  Secondary interaction terms for 2-variable models

	--cathode_material_grad_pen 		(legacy mode: "kappa")  [J/m]
	  Gradient penalty, J/m

	--cathode_material_stress_coeff 	(legacy mode: "B")  [Pa]
	  Stress coefficient, Pa 

	--cathode_material_EvdW 			(legacy mode: "EvdW")  [J/Li]
	  Van der Waals-like interaction energy for graphite models

	--cathode_material_rho_s 			(legacy mode: "rho_s")  [Li/m^3]
	  Total Li site density within the solid, 

	--cathode_material_D 				(legacy mode: "D") [m^2/s]
	  Dimensional diffusivity prefactor (used if type = diffn[2], CHR[2]), 

	--cathode_material_Dfunc 			(legacy mode: "Dfunc")  [-]
	  Dimensionless diffusivity function,

	--cathode_material_dgammadc 		(legacy mode: "dgammadc")  [J*m/Li]
	  Change in surface energy with respect to solid concentration, 

	--cathode_material_cwet 			(legacy mode: "cwet")  [-]
	  ACR dimensionless surface wetting (used with LiFePO4, type=ACR)

  ******** Cathode Reaction Parameters ********

	--cathode_reaction_type 			(legacy mode: "rxnType") [-]
	  reaction type (BV,BV_gMod01,BV_raw,BV_mod01,BV_mod02,Marcus,MHC)

	--cathode_reaction_ecd 			(legacy mode: "k0")  [A/m^2]
	  Exchange current density,  (0.1 for h2/Pt)

	--cathode_reaction_charget 		(legacy mode: "alpha")  [-]
	  Charge transfer coefficient (only used for BV)

	--cathode_reaction_Rorg_energy 	(legacy mode: "lambda")  [J/Li]
	  Reorganizational energy,  (used for Marcus and MHC)

	--cathode_reaction_filmres 		(legacy mode: "Rfilm")  [Ohm m^2]
	  Film resistance, 


  ******** Simulation Parameters ********

	--sim_pars_prof 				(legacy mode: "profileType") [-]
	  Constant voltage or current or segments of one of them

	--sim_pars_crate 				(legacy mode: "Crate") [number of capacities / hr]
	  Battery (dis)charge c-rate (only used for CC), 

	--sim_pars_Vcutmax 				(legacy mode: "Vmax")   [V]
	  Maximum Voltage cutoff, V

	--sim_pars_Vcutmin				(legacy mode: "Vmin")   [V]
	  Minimum Voltage cutoff, V

	--sim_pars_capfrac 				(legacy mode: "capFrac")  [-]
	  Fraction of full (dis)charge to simulate (only used for CC)

	--sim_pars_applv 				(legacy mode: "Vset")  [V]
	  Battery applied voltage (only used for CV), V

	--sim_pars_segments 			(legacy mode: "segments") [-]
	  CC/CV segments defining profile for profileType = CCsegments or CVsegments

	--sim_pars_rampt 				(legacy mode: "tramp") [s]
	  Ramp time to linearly transition to each new segment for CCsegments/CVsegments,

	--sim_pars_fintime 				(legacy mode: "tend")  [s]
	  Final time (only used for CV), 

	--sim_pars_numdisc 				(legacy mode: "tsteps")  [-]
	  The output will have all simulation values at a linear spacing between initial 
	  and final times with tsteps total outputs.

	--sim_pars_reltol 				(legacy mode: "relTol")  [-]
	  Relative Tolerance

	--sim_pars_abstol 				(legacy mode: "absTol")   [-]
	  Absolute Tolerance

	--sim_pars_temp 				(legacy mode: "T")  [K]
	  Temp, 

	--sim_pars_randseed 			(legacy mode: "randomSeed")  [-]
	  Random seed. Set to true to give a random seed in the simulation 
	  (affects noise, particle size distribution)

	--sim_pars_seed 				(legacy mode: "seed")  [-]
	  Value of the random seed, must be an integer

	--sim_pars_Rseries 				(legacy mode: "Rser")	[Ohm m^2]
	  Series resistance, 

	--sim_pars_cat_numdisc 			(legacy mode: "Nvol_c")	[-]
	--sim_pars_sep_numdisc 			(legacy mode: "Nvol_s")	[-]
	--sim_pars_an_numdisc 			(legacy mode: "Nvol_a")	[-]
	  Cathode, seperator, and anode numer disc. in x direction(volumes in electrodes)
	  Nvol_c must be >= 1.If Nvol_c = 1 & Nvol_a = 0 & Nvol_s = 0, simulate a single 
	  volume with no separator and infinitely fast counter electrode. If Nvol_a = 0, 
	  simulate a Li foil electrode

	--sim_pars_cat_nPartVol 		(legacy mode: "Npart_c")	[-]
	  Number of particles per volume for cathode

	--sim_pars_an_nPartVol 			(legacy mode: "Npart_a")	[-]
	  Number of particles per volume for anode

  ******** Electrode Parameters ********

	--sim_elect_rc_Lifoil 			(legacy mode: "k0_foil")	[A/m^2]
	  Rate constant of the Li foil electrode, Used only if Nvol_a = 0

	--sim_elect_Lifoil_filmres 		(legacy mode: "Rfilm_foil")	[Ohm m^2]
	  Film resistance on the Li foil, 


  ******** Particle Distribution Parameters ********

	--sim_PartDist_cat_mean 		(legacy mode: "mean_c")		[m]
	--sim_PartDist_cat_std 			(legacy mode: "stddev_c")	[m]
	--sim_PartDist_an_mean 			(legacy mode: "mean_a")		[m]
	--sim_PartDist_an_std 			(legacy mode: "stddev_a")	[m]
	  Electrode particle size distribution info,  C3 -- size along [100] direction
	  sphere or cylinder -- radius. If using stddev = 0, set Npart = 1 for that 
	  electrode to avoid wasting computational effort.

	--sim_PartDist_cat_fillfrac 	(legacy mode: "cs0_c")		[-]
	--sim_PartDist_an_fillfrac 		(legacy mode: "cs0_a")		[-]
	  Initial electrode filling fractions (for disch, anode starts full, 
	  cathode starts empty)

  ******** Conductivity Parameters ********

	--sim_cond_cat_bulk_par 		(legacy mode: "simBulkCond_c")	[-]
	--sim_cond_an_bulk_par 			(legacy mode: "simBulkCond_a")	[-]
	  Simulate bulk cathode conductivity (Ohms Law)(True/False) 

	--sim_cond_cat_bulk 			(legacy mode: "sigma_s_c")	[S/m]
	--sim_cond_an_bulk 				(legacy mode: "sigma_s_a")	[S/m]
	  Dimensional conductivity (used if simBulkCond = true), 

	--sim_cond_cat_connect_par 		(legacy mode: "simPartCond_c")	[-]
	--sim_cond_an_connect_par 		(legacy mode: "simPartCond_a")	[-]
	  Simulate particle connectivity losses (Ohms Law)(True/False)

	--sim_condpart_cat_mean 		(legacy mode: "G_mean_c")	[1/Ohm]
	--sim_condpart_cat_std 			(legacy mode: "G_stddev_c")	[1/Ohm]
	--sim_condpart_an_mean 			(legacy mode: "G_mean_a")	[1/Ohm]	
	--sim_condpart_an_std 			(legacy mode: "G_stddev_a")	[1/Ohm]
	  Conductance between particles,

  ******** Geometry Parameters ********

	--sim_geo_cat_thickness 		(legacy mode: "L_c")	[m]
	--sim_geo_an_thickness 			(legacy mode: "L_a")	[m]
	--sim_geo_sep_thickness 		(legacy mode: "L_s")	[m]
	  Thicknesses,  (cathode, anode and seperator)

	--sim_geo_cat_volload 			(legacy mode: "P_L_c")	[-]
	--sim_geo_an_volload 			(legacy mode: "P_L_a")	[-]
	  Volume loading percents of active material (volume fraction of solid
	  that is active material)

	--sim_geo_cat_por 				(legacy mode: "poros_c")	[-]
	--sim_geo_an_por 				(legacy mode: "poros_a")	[-]
	--sim_geo_sep_por 				(legacy mode: "poros_s")	[-]
	  Porosities (liquid volume fraction in each region)

	--sim_geo_cat_Brug 				(legacy mode: "BruggExp_c")	[-]
	--sim_geo_an_Brug 				(legacy mode: "BruggExp_a")	[-]
	--sim_geo_sep_Brug 				(legacy mode: "BruggExp_s")	[-]
	  Bruggeman exponent (tortuosity = porosity^bruggExp)

  ******** Electrolyte properties ********

	--sim_elprop_InConc 			(legacy mode: "c0")	[mol/m^3]
	  Initial electrolyte conc., 

	--sim_elprop_cation_chargeNum 	(legacy mode: "zp")	[-]
	--sim_elprop_anion_chargeNum 	(legacy mode: "zm")	[-]
	  Cation/anion charge number (e.g. 2, -1 for CaCl_2)

	--sim_elprop_disNup 			(legacy mode: "nup")	[-]
	--sim_elprop_disNum 			(legacy mode: "num")	[-]
	  Cation/anion dissociation number (e.g. 1, 2 for CaCl_2)

	--sim_elprop_model 				(legacy mode: "elyteModelType")	[-]
	  Electrolyte model,dilute: assume dilute, binary electrolyte model; phi in
	  electrolyte is an "inner potential" or an "quasi-electrostatic
	  potential" SM: use Stefan-Maxwell model for electrolyte transport; phi in
	  electrolyte is that measured with a Li metal reference electrode
	  relative to some fixed position in solution.

	--sim_elprop_SMpropset 			(legacy mode: "SMset")	[-]
	  Stefan-Maxwell property set (test1, LiClO4_PC, valoen_bernardi)
	  test1: parameter set for testing
	  LiClO4_PC: electrolyte/solvent used in Fuller, Doyle, Newman 1994
	             conductivity taken from dualfoil5.2.f
	  valoen_bernardi: LiPF6 in carbonates as in Bernardi and Go 2011

	--sim_elprop_electrans 			(legacy mode: "n")	[-]
	  Reference electrode (defining the electrolyte potential) information:
	  number of electrons transfered in the reaction, 1 for Li/Li+

	--sim_elprop_Stoicof 			(legacy mode: "sp")	[-]
	  Stoichiometric coefficient of cation, -1 for Li/Li+

	--sim_elprop_catdiff 			(legacy mode: "Dp")		[m^2/s]
	--sim_elprop_andiff 			(legacy mode: "Dm")		[m^2/s]
	  Dilute solution properties (used only if elyteModelType = "dilute")
	  Cation/anion diff, 
****************************************************************************************'''
