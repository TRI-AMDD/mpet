#!/usr/bin/env python3

#Import VKML tools
from VKML.gui.tk import Parameter, Parser, Menu
import sys, os
import shutil
#import subprocess
#import mpet.main as main
toolHomeDir = os.path.dirname(os.path.abspath(__file__))
from UI import params as mpetUIparams
from UI import uihelp
from UI import parfile
param_keyDict=mpetUIparams.mpet_paramsUI_dict

def paramfiles():
	for param in p.params:
		value=param.out
		name=param.name
		if name in param_keyDict.keys():
			#print(name,'=',value)
			ns=name.split('_')
			mpetUIparams.findstrNrep(ns[0],param_keyDict[name],value)	

def RunMPETgui():
	paramfiles()
	if (sim_elect_cat_overpar()=='true'):
		shutil.copy('params_c.cfg','params_c_GUIbackup.cfg')
		shutil.copy(Cathode_config_file.filename,'params_c.cfg')
	if (sim_elect_an_overpar()=='true'):
		shutil.copy('params_a.cfg','params_a_GUIbackup.cfg')
		shutil.copy(Anode_config_file.filename,'params_a.cfg')
	if (sim_elect_sim_overpar()=='true'):
		shutil.copy('params_system.cfg','params_system_GUIbackup.cfg')
		shutil.copy(System_config_file.filename,'params_system.cfg')
	#subprocess.call(["python","mpetrun.py","params_system.cfg"])
	main.main('params_system.cfg')

def RunMPETtext():
	for arg in sys.argv:
		if arg.startswith('--') and '=' in arg:
			arg = arg[2:]
			name = arg.split('=')[0]
			value = arg.split('=')[-1]
			if name in param_keyDict.keys():
				for param in p.params:
					if name == param.name:
						ns=name.split('_')
						mpetUIparams.findstrNrep(ns[0],param_keyDict[name],value)
			elif name in ['Cathode_config_file','Anode_config_file','System_config_file']:
				if name=='Cathode_config_file':
					shutil.copy('params_c.cfg','params_c_UIjustused.cfg')
					shutil.copy(value,'params_c.cfg')
				elif name=='Anode_config_file':
					shutil.copy('params_a.cfg','params_a_UIjustused.cfg')
					print(value)
					shutil.copy(value,'params_a.cfg')
				elif name=='System_config_file':
					shutil.copy('params_system.cfg','params_system_UIjustused.cfg')
					shutil.copy(value,'params_system.cfg')
			else:
				print('Wrong Parameter Name: Check parameter command name @',name)
				sys.exit()
		elif (arg not in ['--help','--quit','--gui']) and (arg.startswith('--') and '=' not in arg):
			print('Wrong option seleted. Please refer to help using "--help" option')
			sys.exit()
	#subprocess.call(["python","mpetrun.py","params_system.cfg"])
	main.main('params_system.cfg')
#Create functions
#def callback():
def anpm1():
	anode_part_pars=anode_particle.out
	p1=Parser(title='MPET simulation',interactive=1)
	mAp = Menu(title='Anode Particle Parameters', parser=p1)
	text_particle_shape='sphere -- spherical particle, All surface area is available for reaction\nC3 -- rectangular particle with b surface available for reaction\ncylinder -- cylindrical particle, reactions on surfaces with radial normal'
	Parameter(name='', menu=mAp, variable='label')
	Parameter(name='Anode Particle Parameters', menu=mAp, variable='label')
	Parameter(name='', menu=mAp, variable='label')
	anode_particle_discretization=Parameter(name='anode_particle_discretization', display_name='Discretization of solid (m):',variable=str, menu=mAp, default='50.0e-8',tooltip='Used for non-homogeneous particle type \n When using phase separating particles, must be less than characteristic interface width,\nlambda_b = (kappa/(rho_s*Omega_a))^(0.5)) [Bazant 2013]')
	if anode_part_pars=='CHR2':
		ANDps_OPT = ('cylinder','sphere')
		anode_particle_shape=Parameter(name='anode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	elif anode_part_pars=='ACR':
		ANDps_OPT = ('C3','')
		anode_particle_shape=Parameter(name='anode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	elif anode_part_pars=='CHR':
		ANDps_OPT = ('cylinder','sphere')
		anode_particle_shape=Parameter(name='anode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	elif anode_part_pars in ['homog','homog2','diffn','homog_sdn']:
		ANDps_OPT = ('sphere','cylinder','C3')
		anode_particle_shape=Parameter(name='anode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	anode_particle_thickness=Parameter(name='anode_particle_thickness', display_name='Particle thickness (m):',variable=str, menu=mAp, default='20.0e-9',tooltip='Cylindrical particle thickness or Length of c3( rectangular) particle')   
	Parameter(name='', menu=mAp, variable='label')
	run_func=Parameter(name=' Set Parameters ', variable='function', default=exit,menu=mAp,button_text='DONE',tooltip='Set the parameters')
	p1()
	return p1
def func1():
    q=anpm1()    
    anode_particle_discretization.out=q.params[2]()
    anode_particle_shape.out=q.params[3]()
    anode_particle_thickness.out=q.params[4]()
def anpm2():
	anode_mat_pars=anode_material_muRfunc.out
	p1=Parser(title='MPET simulation',interactive=1)
	mAm = Menu(title='Anode Material Parameters', parser=p1)
	Parameter(name='', menu=mAm, variable='label')
	Parameter(name='Anode Material Parameters', menu=mAm, variable='label')
	Parameter(name='', menu=mAm, variable='label')
	muRfunc_list=['LiC6_coke_ss2','LiC6','LiMn2O4_ss','LiMn2O4_ss2','LiC6_coke_ss','LiC6_ss','LiC6_ss2','LiC6_2step_ss','Li_ss','NCA_ss1','NCA_ss2','testIS_ss','testRS_ss','ideal_sln','reg_sln','graphite_2param_homog','graphite_1param_homog','graphite_1param_homog_2','non_homog_rect_fixed_csurf','non_homog_round_wetting','general_non_homog','LiFePO4','LiC6_1param','testRS','testRS_ps']
	if anode_mat_pars in muRfunc_list:
		AND_tf=('false','true')
		anode_material_Logpad=Parameter(name='anode_material_Logpad', display_name='Log pad (option for muRfunc containing log terms -- not solid solution) "logpad":', variable=tuple,menu=mAm, default=(AND_tf),tooltip='Consider turning this on if you experience errors with log(negative)\n It extends the ideal solution into regions outside (0, 1) by a \n (very steep) linear extension instead of a singularity.')
		anode_material_noise=Parameter(name='anode_material_noise', display_name='Add Langevin noise to the solid dynamics "noise":', variable=tuple,menu=mAm, default=(AND_tf),tooltip='add Langevin noise to the solid dynamics to simulate some\n random thermal fluctuations. This can help to, e.g., trigger a \n spinodal decomposition but slows simulations so should not be \n enabled unless needed.')
		anode_material_mag_noise=Parameter(name='anode_material_mag_noise', display_name='Magnitude of the noise(if noise=TRUE)"noise_prefac":', variable=str,menu=mAm, default='1e-6',tooltip='Magnitude of the noise (used if noise = true)\n The standard deviation of the normal noise added to the solid\n dynamics')
		anode_material_num_noise=Parameter(name='anode_material_num_noise', display_name='Number of segments over which the noise is re-randomized "numnoise":', variable=str,menu=mAm, default='200',tooltip='Number of segments over which the noise is re-randomized')
		anode_material_oma=Parameter(name='anode_material_oma', display_name='Regular solution parameter (J/Li) "Omega_a":', variable=str,menu=mAm, default='1.8560e-20',tooltip='# The following parameters are used for some of the materials defined \nby muRfunc. All must be defined in this input file. For solid \nsolution materials (containing "_ss"), none of the following are used \n except D \n Note k*298 K = 4.115e-21 J/particle \n Regular solution parameter, J/Li \nDefault = 4.51 (k*298 K)')	
		anode_material_omb=Parameter(name='anode_material_omb', display_name='Secondary interaction terms for 2-variable models "Omega_b":', variable=str,menu=mAm, default='5.761532e-21',tooltip=' Secondary interaction terms for 2-variable models\n See, e.g., Ferguson and Bazant 2014')
		anode_material_omc=Parameter(name='anode_material_omc', display_name='Secondary interaction terms for 2-variable models "Omega_c":', variable=str,menu=mAm, default='8.23076e-20',tooltip=' Secondary interaction terms for 2-variable models\n See, e.g., Ferguson and Bazant 2014')	
		anode_material_grad_pen=Parameter(name='anode_material_grad_pen', display_name='Gradient penality (J/m) "kappa":', variable=str,menu=mAm, default='5.0148e-10',tooltip='Gradient penalty, J/m')
		anode_material_stress_coeff=Parameter(name='anode_material_stress_coeff', display_name='Stress Coeffecient (Pa) "B":', variable=str,menu=mAm, default='0.1916e9',tooltip='Stress coefficient, Pa')	
		anode_material_EvdW=Parameter(name='anode_material_EvdW', display_name='Van der Waals-like interaction energy for graphite models "EvdW":', variable=str,menu=mAm, default='0.0',tooltip='Van der Waals-like interaction energy for graphite models \ndefault: 0.0')		
		anode_material_rho_s=Parameter(name='anode_material_rho_s', display_name='Total Li site density within the solid (Li/m^3) "rho_s":', variable=str,menu=mAm, default='1.58981e28',tooltip='Total Li site density within the solid, Li/m^3 \n Even for 2-parameter models, this should be the full site density \n(not that of an individual layer/lattice)')	
		anode_material_D=Parameter(name='anode_material_D', display_name='Dimensional diffusivity prefactor (used if type = diffn[2], CHR[2]) (m^2/s) "D":', variable=str,menu=mAm, default='5e-13',tooltip='Dimensional diffusivity prefactor (used if type = diffn[2], CHR[2]), m^2/s')	
		AND_dfunc=('constant','lattice')
		anode_material_Dfunc=Parameter(name='anode_material_Dfunc', display_name='Dimensionless diffusivity function "Dfunc":', variable=tuple,menu=mAm, default=(AND_dfunc),tooltip='# Dimensionless diffusivity function, (see props_am.Dfuncs class) \n This is a function of solid filling fraction, y, such that \n Flux = -D*Dfunc(y)*grad(c) for solid solution materials, and \nFlux = -D*Dfunc(y)*grad(mu) for thermo-based transport models')
		anode_material_dgammadc=Parameter(name='anode_material_dgammadc', display_name='Surface wetting (used if type = CHR[2]) "dgammadc":', variable=str,menu=mAm, default='0e-30',tooltip='Surface wetting (used if type = CHR[2]).\n Change in surface energy with respect to solid concentration, J*m/Li \n set greater than 0 for surface wetting, less than 0 for dewetting \n default: 0e-30')	
		anode_material_cwet=Parameter(name='anode_material_cwet', display_name='ACR dimensionless surface wetting (used with LiFePO4, type=ACR) "cwet":', variable=str,menu=mAm, default='0.98',tooltip='ACR dimensionless surface wetting (used with LiFePO4, type=ACR)')	
	run_func=Parameter(name=' Set Parameters ', variable='function', default=exit,menu=mAm,button_text='DONE',tooltip='Set the parameters')
	p1()
	return p1
def func2():
	q=anpm2()
	anode_material_Logpad.out=q.params[2]()
	anode_material_noise.out=q.params[3]()
	anode_material_mag_noise.out=q.params[4]()
	anode_material_num_noise.out=q.params[5]()
	anode_material_oma.out=q.params[6]()
	anode_material_omb.out=q.params[7]()
	anode_material_omc.out=q.params[8]()
	anode_material_grad_pen.out=q.params[9]()
	anode_material_stress_coeff.out=q.params[10]()
	anode_material_EvdW.out=q.params[11]()
	anode_material_rho_s.out=q.params[12]()
	anode_material_D.out=q.params[13]()
	anode_material_Dfunc.out=q.params[14]()
	anode_material_dgammadc.out=q.params[15]()
	anode_material_cwet.out=q.params[16]()

def anpm3():
	anode_rxn_pars=anode_reaction_type.out
	p1=Parser(title='MPET simulation',interactive=1)
	mAr = Menu(title='Anode Particle Parameters', parser=p1)
	Parameter(name='', menu=mAr, variable='label')
	Parameter(name='Anode Reaction Parameters', menu=mAr, variable='label')
	Parameter(name='', menu=mAr, variable='label')
	rxn_list=['BV_mod02','BV_gMod01','BV','BV_raw','BV_mod01','Marcus','MHC']
	if anode_rxn_pars in rxn_list:
		anode_reaction_ecd=Parameter(name='anode_reaction_ecd', display_name='Exchange current density (A/m^2) "k0":',variable=str, menu=mAr, default='8.3038')
		anode_reaction_charget=Parameter(name='anode_reaction_charget', display_name='Charge transfer coefficient (only used for BV) "alpha":',variable=str, menu=mAr, default='0.5')
		anode_reaction_Rorg_energy=Parameter(name='anode_reaction_Rorg_energy', display_name='Reorganizational energy (J/Li) (used for Marcus and MHC) "lambda":',variable=str, menu=mAr, default='6.26e-20')
		anode_reaction_filmres=Parameter(name='anode_reaction_filmres', display_name='Film resistance (Ohm m^2) "Rfilm":',variable=str, menu=mAr, default='0e-0')
	Parameter(name='', menu=mAr, variable='label')
	run_func=Parameter(name=' Set Parameters ', variable='function', default=exit,menu=mAr,button_text='DONE')
	p1()
	return p1
def func3():
	q=anpm3()
	anode_reaction_ecd.out=q.params[2]()
	anode_reaction_charget.out=q.params[3]()
	anode_reaction_Rorg_energy.out=q.params[4]()
	anode_reaction_filmres.out=q.params[5]()   

def catpm1():
	cathode_part_pars=cathode_particle.out
	p1=Parser(title='MPET simulation',interactive=1)
	mAp = Menu(title='Cathode Particle Parameters', parser=p1)
	text_particle_shape='sphere -- spherical particle, All surface area is available for reaction\nC3 -- rectangular particle with b surface available for reaction\ncylinder -- cylindrical particle, reactions on surfaces with radial normal'
	Parameter(name='', menu=mAp, variable='label')
	Parameter(name='Cathode Particle Parameters', menu=mAp, variable='label')
	Parameter(name='', menu=mAp, variable='label')
	cathode_particle_discretization=Parameter(name='cathode_particle_discretization', display_name='Discretization of solid (m):',variable=str, menu=mAp, default='50.0e-8',tooltip='Used for non-homogeneous particle type \n When using phase separating particles, must be less than characteristic interface width,\nlambda_b = (kappa/(rho_s*Omega_a))^(0.5)) [Bazant 2013]')
	if cathode_part_pars=='CHR2':
		ANDps_OPT = ('cylinder','sphere')
		cathode_particle_shape=Parameter(name='cathode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	elif cathode_part_pars=='ACR':
		ANDps_OPT = ('C3','')
		cathode_particle_shape=Parameter(name='cathode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	elif cathode_part_pars=='CHR':
		ANDps_OPT = ('cylinder','sphere')
		cathode_particle_shape=Parameter(name='cathode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	elif cathode_part_pars in ['homog','homog2','diffn','homog_sdn']:
		ANDps_OPT = ('sphere','cylinder','C3')
		cathode_particle_shape=Parameter(name='cathode_particle_shape', display_name='Choose Particle Shape:',variable=tuple, menu=mAp, default=(ANDps_OPT),tooltip=text_particle_shape)
	cathode_particle_thickness=Parameter(name='cathode_particle_thickness', display_name='Particle thickness (m):',variable=str, menu=mAp, default='20.0e-9',tooltip='Cylindrical particle thickness or Length of c3( rectangular) particle')   
	Parameter(name='', menu=mAp, variable='label')
	run_func=Parameter(name=' Set Parameters ', variable='function', default=exit,menu=mAp,button_text='DONE',tooltip='Set the parameters')
	p1()
	return p1
def func4():
    q=catpm1()    
    cathode_particle_discretization.out=q.params[2]()
    cathode_particle_shape.out=q.params[3]()
    cathode_particle_thickness.out=q.params[4]()
def catpm2():
	cathode_mat_pars=cathode_material_muRfunc.out
	p1=Parser(title='MPET simulation',interactive=1)
	mAm = Menu(title='Cathode Material Parameters', parser=p1)
	Parameter(name='', menu=mAm, variable='label')
	Parameter(name='Cathode Material Parameters', menu=mAm, variable='label')
	Parameter(name='', menu=mAm, variable='label')
	muRfunc_list=['LiC6_coke_ss2','LiC6','LiMn2O4_ss','LiMn2O4_ss2','LiC6_coke_ss','LiC6_ss','LiC6_ss2','LiC6_2step_ss','Li_ss','NCA_ss1','NCA_ss2','testIS_ss','testRS_ss','ideal_sln','reg_sln','graphite_2param_homog','graphite_1param_homog','graphite_1param_homog_2','non_homog_rect_fixed_csurf','non_homog_round_wetting','general_non_homog','LiFePO4','LiC6_1param','testRS','testRS_ps']
	if cathode_mat_pars in muRfunc_list:
		AND_tf=('false','true')
		cathode_material_Logpad=Parameter(name='cathode_material_Logpad', display_name='Log pad (option for muRfunc containing log terms -- not solid solution) "logpad":', variable=tuple,menu=mAm, default=(AND_tf),tooltip='Consider turning this on if you experience errors with log(negative)\n It extends the ideal solution into regions outside (0, 1) by a \n (very steep) linear extension instead of a singularity.')
		cathode_material_noise=Parameter(name='cathode_material_noise', display_name='Add Langevin noise to the solid dynamics "noise":', variable=tuple,menu=mAm, default=(AND_tf),tooltip='add Langevin noise to the solid dynamics to simulate some\n random thermal fluctuations. This can help to, e.g., trigger a \n spinodal decomposition but slows simulations so should not be \n enabled unless needed.')
		cathode_material_mag_noise=Parameter(name='cathode_material_mag_noise', display_name='Magnitude of the noise(if noise=TRUE)"noise_prefac":', variable=str,menu=mAm, default='1e-6',tooltip='Magnitude of the noise (used if noise = true)\n The standard deviation of the normal noise added to the solid\n dynamics')
		cathode_material_num_noise=Parameter(name='cathode_material_num_noise', display_name='Number of segments over which the noise is re-randomized "numnoise":', variable=str,menu=mAm, default='200',tooltip='Number of segments over which the noise is re-randomized')
		cathode_material_oma=Parameter(name='cathode_material_oma', display_name='Regular solution parameter (J/Li) "Omega_a":', variable=str,menu=mAm, default='1.8560e-20',tooltip='# The following parameters are used for some of the materials defined \nby muRfunc. All must be defined in this input file. For solid \nsolution materials (containing "_ss"), none of the following are used \n except D \n Note k*298 K = 4.115e-21 J/particle \n Regular solution parameter, J/Li \nDefault = 4.51 (k*298 K)')	
		cathode_material_omb=Parameter(name='cathode_material_omb', display_name='Secondary interaction terms for 2-variable models "Omega_b":', variable=str,menu=mAm, default='5.761532e-21',tooltip=' Secondary interaction terms for 2-variable models\n See, e.g., Ferguson and Bazant 2014')
		cathode_material_omc=Parameter(name='cathode_material_omc', display_name='Secondary interaction terms for 2-variable models "Omega_c":', variable=str,menu=mAm, default='8.23076e-20',tooltip=' Secondary interaction terms for 2-variable models\n See, e.g., Ferguson and Bazant 2014')	
		cathode_material_grad_pen=Parameter(name='cathode_material_grad_pen', display_name='Gradient penality (J/m) "kappa":', variable=str,menu=mAm, default='5.0148e-10',tooltip='Gradient penalty, J/m')
		cathode_material_stress_coeff=Parameter(name='cathode_material_stress_coeff', display_name='Stress Coeffecient (Pa) "B":', variable=str,menu=mAm, default='0.1916e9',tooltip='Stress coefficient, Pa')	
		cathode_material_EvdW=Parameter(name='cathode_material_EvdW', display_name='Van der Waals-like interaction energy for graphite models "EvdW":', variable=str,menu=mAm, default='0.0',tooltip='Van der Waals-like interaction energy for graphite models \ndefault: 0.0')		
		cathode_material_rho_s=Parameter(name='cathode_material_rho_s', display_name='Total Li site density within the solid (Li/m^3) "rho_s":', variable=str,menu=mAm, default='1.58981e28',tooltip='Total Li site density within the solid, Li/m^3 \n Even for 2-parameter models, this should be the full site density \n(not that of an individual layer/lattice)')	
		cathode_material_D=Parameter(name='cathode_material_D', display_name='Dimensional diffusivity prefactor (used if type = diffn[2], CHR[2]) (m^2/s) "D":', variable=str,menu=mAm, default='5e-13',tooltip='Dimensional diffusivity prefactor (used if type = diffn[2], CHR[2]), m^2/s')	
		AND_dfunc=('constant','lattice')
		cathode_material_Dfunc=Parameter(name='cathode_material_Dfunc', display_name='Dimensionless diffusivity function "Dfunc":', variable=tuple,menu=mAm, default=(AND_dfunc),tooltip='# Dimensionless diffusivity function, (see props_am.Dfuncs class) \n This is a function of solid filling fraction, y, such that \n Flux = -D*Dfunc(y)*grad(c) for solid solution materials, and \nFlux = -D*Dfunc(y)*grad(mu) for thermo-based transport models')
		cathode_material_dgammadc=Parameter(name='cathode_material_dgammadc', display_name='Surface wetting (used if type = CHR[2]) "dgammadc":', variable=str,menu=mAm, default='0e-30',tooltip='Surface wetting (used if type = CHR[2]).\n Change in surface energy with respect to solid concentration, J*m/Li \n set greater than 0 for surface wetting, less than 0 for dewetting \n default: 0e-30')	
		cathode_material_cwet=Parameter(name='cathode_material_cwet', display_name='ACR dimensionless surface wetting (used with LiFePO4, type=ACR) "cwet":', variable=str,menu=mAm, default='0.98',tooltip='ACR dimensionless surface wetting (used with LiFePO4, type=ACR)')	
	run_func=Parameter(name=' Set Parameters ', variable='function', default=exit,menu=mAm,button_text='DONE',tooltip='Set the parameters')
	p1()
	return p1
def func5():
	q=catpm2()
	cathode_material_Logpad.out=q.params[2]()
	cathode_material_noise.out=q.params[3]()
	cathode_material_mag_noise.out=q.params[4]()
	cathode_material_num_noise.out=q.params[5]()
	cathode_material_oma.out=q.params[6]()
	cathode_material_omb.out=q.params[7]()
	cathode_material_omc.out=q.params[8]()
	cathode_material_grad_pen.out=q.params[9]()
	cathode_material_stress_coeff.out=q.params[10]()
	cathode_material_EvdW.out=q.params[11]()
	cathode_material_rho_s.out=q.params[12]()
	cathode_material_D.out=q.params[13]()
	cathode_material_Dfunc.out=q.params[14]()
	cathode_material_dgammadc.out=q.params[15]()
	cathode_material_cwet.out=q.params[16]()

def catpm3():
	cathode_rxn_pars=cathode_reaction_type.out
	p1=Parser(title='MPET simulation',interactive=1)
	mAr = Menu(title='Cathode Particle Parameters', parser=p1)
	Parameter(name='', menu=mAr, variable='label')
	Parameter(name='Cathode Reaction Parameters', menu=mAr, variable='label')
	Parameter(name='', menu=mAr, variable='label')
	rxn_list=['BV_mod02','BV_gMod01','BV','BV_raw','BV_mod01','Marcus','MHC']
	if cathode_rxn_pars in rxn_list:
		cathode_reaction_ecd=Parameter(name='cathode_reaction_ecd', display_name='Exchange current density (A/m^2) "k0":',variable=str, menu=mAr, default='8.3038')
		cathode_reaction_charget=Parameter(name='cathode_reaction_charget', display_name='Charge transfer coefficient (only used for BV) "alpha":',variable=str, menu=mAr, default='0.5')
		cathode_reaction_Rorg_energy=Parameter(name='cathode_reaction_Rorg_energy', display_name='Reorganizational energy (J/Li) (used for Marcus and MHC) "lambda":',variable=str, menu=mAr, default='6.26e-20')
		cathode_reaction_filmres=Parameter(name='cathode_reaction_filmres', display_name='Film resistance (Ohm m^2) "Rfilm":',variable=str, menu=mAr, default='0e-0')
	Parameter(name='', menu=mAr, variable='label')
	run_func=Parameter(name=' Set Parameters ', variable='function', default=exit,menu=mAr,button_text='DONE')
	p1()
	return p1
def func6():
	q=catpm3()
	cathode_reaction_ecd.out=q.params[2]()
	cathode_reaction_charget.out=q.params[3]()
	cathode_reaction_Rorg_energy.out=q.params[4]()
	cathode_reaction_filmres.out=q.params[5]() 

#Create parser
p=Parser(title='MPET simulation',help_text=uihelp.UIhelp)

#Create menus
#welcome = Menu(title='Welcome', parser=p)
mAp = Menu(title='Anode Particle Parameters', parser=p)
mAm = Menu(title='Anode Material Parameters', parser=p)
mAr = Menu(title='Anode Reaction Parameters', parser=p)

mCp = Menu(title='Cathode Particle Parameters', parser=p)
mCm = Menu(title='Cathode Material Parameters', parser=p)
mCr = Menu(title='Cathode Reaction Parameters', parser=p)


mSp1 = Menu(title='Simulation Parameters', parser=p)
mSp2 = Menu(title='Simulation Parameters contd.', parser=p)
mSoe = Menu(title='Simulation Electrode Parameters',parser=p)
mSpd = Menu(title='Simulation Particles Distribution',parser=p)
mSc = Menu(title='Simulation Conductivities',parser=p)
mSg = Menu(title='Simulation Geometry',parser=p)
mSel = Menu(title='Simulation Electrolyte Properties',parser=p)
mSlaunch = Menu(title='MPET simulation launch',parser=p)
#PARAMETERS FOR WELCOME MENU

#Parameters for 'Anode Particle Properties'
Parameter(name='', menu=mAp, variable='label')
Parameter(name='Anode Particle Parameters', menu=mAp, variable='label')
Parameter(name='', menu=mAp, variable='label')
ANDspt_OPT_txt='CHR2: 2-concentration field Cahn-Hilliard-Reaction\n ACR: Allen-Cahn-Reaction particles\nCHR: Cahn-Hilliard-Reaction\n homog: homogeneously filling particles\n homog_sdn_homogeneously filling particles with size-dependent nucleation\nhomog2: 2-concentration field model with each layer treated as homogeneous.\n diffn: Particles with rxn at surface and internal solid diffusion'
ANDspt_OPT_pars = ('diffn','CHR2','ACR','CHR','homog','homog_sdn','homog2')
anode_particle=Parameter(name='anode_particle', display_name='Choose types of solid particles:', variable=tuple,menu=mAp, default=(ANDspt_OPT_pars),tooltip=ANDspt_OPT_txt)
anode_particle_discretization=Parameter(name='anode_particle_discretization',display_name=' ',variable=str, menu=mAp, default='50.0e-8',tooltip='Used for non-homogeneous particle type \n When using phase separating particles, must be less than characteristic interface width,\nlambda_b = (kappa/(rho_s*Omega_a))^(0.5)) [Bazant 2013]',show_par='False')
run_function=Parameter(name='  Set Parameters ', variable='function', default=func1,menu=mAp,button_text='SET',tooltip='Open new menu containing respective params')
anode_particle_shape=Parameter(name='anode_particle_shape',display_name=' ',variable=str, menu=mAp, default='sphere',tooltip='text_particle_shape',show_par='False')
anode_particle_thickness=Parameter(name='anode_particle_thickness',display_name=' ',variable=str, menu=mAp, default='20.0e-9',show_par='False')
#Parameters for 'Anode Material Properties'

Parameter(name='', menu=mAm, variable='label')
Parameter(name='Anode Material Parameters', menu=mAm, variable='label')
Parameter(name='', menu=mAm, variable='label')

anode_material_Logpad=Parameter(name='anode_material_Logpad',display_name=' ',  variable=str,\
              menu=mAm, default='false',show_par=False)
anode_material_noise=Parameter(name='anode_material_noise',display_name=' ', variable=str,\
              menu=mAm, default='false',show_par=False)
anode_material_mag_noise=Parameter(name='anode_material_mag_noise',display_name=' ', variable=str,\
              menu=mAm, default='1e-6',show_par=False)
anode_material_num_noise=Parameter(name='anode_material_num_noise',display_name=' ', variable=str,\
              menu=mAm, default='200',show_par=False)
anode_material_oma=Parameter(name='anode_material_oma',display_name=' ', variable=str,\
              menu=mAm, default='1.8560e-20',show_par=False)
AND_murf = ('LiC6_coke_ss2','LiC6','LiMn2O4_ss','LiMn2O4_ss2','LiC6_coke_ss','LiC6_ss','LiC6_ss2','LiC6_2step_ss','Li_ss','NCA_ss1','NCA_ss2','testIS_ss','testRS_ss','ideal_sln','reg_sln','graphite_2param_homog','graphite_1param_homog','graphite_1param_homog_2','non_homog_rect_fixed_csurf','non_homog_round_wetting','general_non_homog','LiFePO4','LiC6_1param','testRS','testRS_ps')
anode_material_muRfunc=Parameter(name='anode_material_muRfunc', display_name='Choose muRfunctions :', variable=tuple,\
              menu=mAm, default=(AND_murf),tooltip='muRfunc defines the chemical potential of the reduced state\n(the electron conducting phase, usually solid) as a function of\nfilling fraction.')	
anode_material_omb=Parameter(name='anode_material_omb',display_name=' ',  variable=str,\
              menu=mAm, default='5.761532e-21',show_par=False)
anode_material_omc=Parameter(name='anode_material_omc',display_name=' ',  variable=str,\
              menu=mAm, default='8.23076e-20',show_par=False)	
anode_material_grad_pen=Parameter(name='anode_material_grad_pen',display_name=' ',  variable=str,\
              menu=mAm, default='5.0148e-10',show_par=False)
anode_material_stress_coeff=Parameter(name='anode_material_stress_coeff',display_name=' ',  variable=str,\
              menu=mAm, default='0.1916e9',show_par=False)	
anode_material_EvdW=Parameter(name='anode_material_EvdW',display_name=' ',  variable=str,\
              menu=mAm, default='0.0',show_par=False)
run_function=Parameter(name='  Set Parameters ', variable='function', default=func2,menu=mAm,button_text='SET',tooltip='Open new menu containing respective params') 		
anode_material_rho_s=Parameter(name='anode_material_rho_s',display_name=' ',  variable=str,\
              menu=mAm, default='1.58981e28',show_par=False)	
anode_material_D=Parameter(name='anode_material_D',display_name=' ', variable=str,\
              menu=mAm, default='5e-13',show_par=False)	

anode_material_Dfunc=Parameter(name='anode_material_Dfunc',display_name=' ',  variable=str,\
              menu=mAm, default='constant',show_par=False)
anode_material_dgammadc=Parameter(name='anode_material_dgammadc',display_name=' ', variable=str,\
              menu=mAm, default='0e-30',show_par=False)	
anode_material_cwet=Parameter(name='anode_material_cwet',display_name=' ',  variable=str,\
              menu=mAm, default='0.98',show_par=False)
#Parameters for 'Anode Reaction Properties'

Parameter(name='', menu=mAr, variable='label')
Parameter(name='Anode Reaction Parameters', menu=mAr, variable='label')
Parameter(name='', menu=mAr, variable='label')
ANDrxn_OPT = ('BV_mod02','BV_gMod01','BV','BV_raw','BV_mod01','Marcus','MHC')
ANDrxn_txt='#   BV: Butler-Volmer from Bazant Acc. Chem. Res. 2013 with gamma_TS = 1/(1-c)\n#   BV_gMod01: Like BV with gamma_TS = 1/(c*(1-c))\n#   BV_raw: BV with ecd = k0\n#   BV_mod01: As in Fuller, Doyle, Newman 1994 cathode\n#   BV_mod02: As in Fuller, Doyle, Newman 1994 anode\n#   Marcus: As in Bazant Acc. Chem. Res. 2013 \n#   MHC: See Zeng, Smith, Bai, Bazant, 2014'
anode_reaction_type=Parameter(name='anode_reaction_type', display_name='Choose Reaction type:', variable=tuple,\
              menu=mAr, default=(ANDrxn_OPT),tooltip=ANDrxn_txt)
anode_reaction_ecd=Parameter(name='anode_reaction_ecd',display_name=' ',\
              variable=str, menu=mAr, default='8.3038',show_par=False)
anode_reaction_charget=Parameter(name='anode_reaction_charget',display_name=' ', \
              variable=str, menu=mAr, default='0.5',show_par=False)
run_function=Parameter(name='  Set Parameters ', variable='function', default=func3,menu=mAr,button_text='SET',tooltip='Open new menu containing respective params')
anode_reaction_Rorg_energy=Parameter(name='anode_reaction_Rorg_energy',display_name=' ',\
              variable=str, menu=mAr, default='6.26e-20',show_par=False)
anode_reaction_filmres=Parameter(name='anode_reaction_filmres',display_name=' ',\
              variable=str, menu=mAr, default='0e-0',show_par=False)
Parameter(name='', menu=mAr, variable='label')			  

#Parameters for 'cathode Particle Properties'
Parameter(name='', menu=mCp, variable='label')
Parameter(name='cathode Particle Parameters', menu=mCp, variable='label')
Parameter(name='', menu=mCp, variable='label')
ANDspt_OPT_txt='CHR2: 2-concentration field Cahn-Hilliard-Reaction\n ACR: Allen-Cahn-Reaction particles\nCHR: Cahn-Hilliard-Reaction\n homog: homogeneously filling particles\n homog_sdn_homogeneously filling particles with size-dependent nucleation\nhomog2: 2-concentration field model with each layer treated as homogeneous.\n diffn: Particles with rxn at surface and internal solid diffusion'
ANDspt_OPT_pars = ('diffn','ACR','CHR2','CHR','homog','homog_sdn','homog2')
cathode_particle=Parameter(name='cathode_particle', display_name='Choose types of solid particles:', variable=tuple,menu=mCp, default=(ANDspt_OPT_pars),tooltip=ANDspt_OPT_txt)
cathode_particle_discretization=Parameter(name='cathode_particle_discretization',display_name=' ',variable=str, menu=mCp, default='5.0e-8',tooltip='Used for non-homogeneous particle type \n When using phase separating particles, must be less than characteristic interface width,\nlambda_b = (kappa/(rho_s*Omega_a))^(0.5)) [Bazant 2013]',show_par='False')
run_function=Parameter(name='  Set Parameters ', variable='function', default=func4,menu=mCp,button_text='SET',tooltip='Open new menu containing respective params')
cathode_particle_shape=Parameter(name='cathode_particle_shape',display_name=' ',variable=str, menu=mCp, default='sphere',tooltip='text_particle_shape',show_par='False')
cathode_particle_thickness=Parameter(name='cathode_particle_thickness',display_name=' ',variable=str, menu=mCp, default='20.0e-9',show_par='False')
#Parameters for 'cathode Material Properties'

Parameter(name='', menu=mCm, variable='label')
Parameter(name='cathode Material Parameters', menu=mCm, variable='label')
Parameter(name='', menu=mCm, variable='label')

cathode_material_Logpad=Parameter(name='cathode_material_Logpad',display_name=' ',  variable=str,\
              menu=mCm, default='false',show_par=False)
cathode_material_noise=Parameter(name='cathode_material_noise',display_name=' ', variable=str,\
              menu=mCm, default='false',show_par=False)
cathode_material_mag_noise=Parameter(name='cathode_material_mag_noise',display_name=' ', variable=str,\
              menu=mCm, default='1e-6',show_par=False)
cathode_material_num_noise=Parameter(name='cathode_material_num_noise',display_name=' ', variable=str,\
              menu=mCm, default='200',show_par=False)
cathode_material_oma=Parameter(name='cathode_material_oma',display_name=' ', variable=str,\
              menu=mCm, default='1.8560e-20',show_par=False)
AND_murf = ('LiMn2O4_ss2','LiFePO4','LiC6','LiMn2O4_ss','LiC6_coke_ss','LiC6_coke_ss2','LiC6_ss','LiC6_ss2','LiC6_2step_ss','Li_ss','NCA_ss1','NCA_ss2','testIS_ss','testRS_ss','ideal_sln','reg_sln','graphite_2param_homog','graphite_1param_homog','graphite_1param_homog_2','non_homog_rect_fixed_csurf','non_homog_round_wetting','general_non_homog','LiFePO4','LiC6_1param','testRS','testRS_ps')
cathode_material_muRfunc=Parameter(name='cathode_material_muRfunc', display_name='Choose muRfunctions :', variable=tuple,\
              menu=mCm, default=(AND_murf),tooltip='muRfunc defines the chemical potential of the reduced state\n(the electron conducting phase, usually solid) as a function of\nfilling fraction.')	
cathode_material_omb=Parameter(name='cathode_material_omb',display_name=' ',  variable=str,\
              menu=mCm, default='5.761532e-21',show_par=False)
cathode_material_omc=Parameter(name='cathode_material_omc',display_name=' ',  variable=str,\
              menu=mCm, default='8.23076e-20',show_par=False)	
cathode_material_grad_pen=Parameter(name='cathode_material_grad_pen',display_name=' ',  variable=str,\
              menu=mCm, default='5.0148e-10',show_par=False)
cathode_material_stress_coeff=Parameter(name='cathode_material_stress_coeff',display_name=' ',  variable=str,\
              menu=mCm, default='0.1916e9',show_par=False)	
cathode_material_EvdW=Parameter(name='cathode_material_EvdW',display_name=' ',  variable=str,\
              menu=mCm, default='0.0',show_par=False)
run_function=Parameter(name='  Set Parameters ', variable='function', default=func5,menu=mCm,button_text='SET',tooltip='Open new menu containing respective params') 		
cathode_material_rho_s=Parameter(name='cathode_material_rho_s',display_name=' ',  variable=str,\
              menu=mCm, default='1.42842e28',show_par=False)	
cathode_material_D=Parameter(name='cathode_material_D',display_name=' ', variable=str,\
              menu=mCm, default='1.0e-13',show_par=False)	

cathode_material_Dfunc=Parameter(name='cathode_material_Dfunc',display_name=' ',  variable=str,\
              menu=mCm, default='constant',show_par=False)
cathode_material_dgammadc=Parameter(name='cathode_material_dgammadc',display_name=' ', variable=str,\
              menu=mCm, default='0e-30',show_par=False)	
cathode_material_cwet=Parameter(name='cathode_material_cwet',display_name=' ',  variable=str,\
              menu=mCm, default='0.98',show_par=False)
#Parameters for 'cathode Reaction Properties'

Parameter(name='', menu=mCr, variable='label')
Parameter(name='cathode Reaction Parameters', menu=mCr, variable='label')
Parameter(name='', menu=mCr, variable='label')
ANDrxn_OPT = ('BV_mod01','BV','BV_gMod01','BV_raw','BV_mod02','Marcus','MHC')
ANDrxn_txt='#   BV: Butler-Volmer from Bazant Acc. Chem. Res. 2013 with gamma_TS = 1/(1-c)\n#   BV_gMod01: Like BV with gamma_TS = 1/(c*(1-c))\n#   BV_raw: BV with ecd = k0\n#   BV_mod01: As in Fuller, Doyle, Newman 1994 cathode\n#   BV_mod02: As in Fuller, Doyle, Newman 1994 cathode\n#   mCrcus: As in Bazant Acc. Chem. Res. 2013 \n#   MHC: See Zeng, Smith, Bai, Bazant, 2014'
cathode_reaction_type=Parameter(name='cathode_reaction_type', display_name='Choose Reaction type:', variable=tuple,\
              menu=mCr, default=(ANDrxn_OPT),tooltip=ANDrxn_txt)
cathode_reaction_ecd=Parameter(name='cathode_reaction_ecd',display_name=' ',\
              variable=str, menu=mCr, default='7.225',show_par=False)
cathode_reaction_charget=Parameter(name='cathode_reaction_charget',display_name=' ', \
              variable=str, menu=mCr, default='0.5',show_par=False)
run_function=Parameter(name='  Set Parameters ', variable='function', default=func6,menu=mCr,button_text='SET',tooltip='Open new menu containing respective params')
cathode_reaction_Rorg_energy=Parameter(name='cathode_reaction_Rorg_energy',display_name=' ',\
              variable=str, menu=mCr, default='6.26e-20',show_par=False)
cathode_reaction_filmres=Parameter(name='cathode_reaction_filmres',display_name=' ',\
              variable=str, menu=mCr, default='0e-0',show_par=False)
Parameter(name='', menu=mCr, variable='label')	

#####################################################
#simULATION PARAMETERS
#####################################################
#simulation Parameters

Parameter(name='', menu=mSp1, variable='label')
Parameter(name='Simulation Parameters', menu=mSp1, variable='label')
Parameter(name='', menu=mSp1, variable='label')

sim_profile=('CC','CV','CCsegments','CVsegments')
sim_pars_prof=Parameter(name='sim_pars_prof', display_name='simulation profile type CC or CV or combination both: ',variable=tuple, menu=mSp1, default=(sim_profile),tooltip='CV- constant voltage\n CC-constant current\n CCsegments- segments of constant current\n CVsegments- segments of constant voltage')
sim_pars_crate=Parameter(name='sim_pars_crate', display_name='Battery charge or discharge c-rate :',variable=str, menu=mSp1, default='0.66',tooltip='(number of capacities/hr) (positive for discharge, negative for charge)')
sim_pars_Vcutmax=Parameter(name='sim_pars_Vcutmax', display_name='Maximum Voltage cutoffs (V)::',variable=str, menu=mSp1, default='5',tooltip='Maximum voltage cutoff (V)')
sim_pars_Vcutmin=Parameter(name='sim_pars_Vcutmin', display_name='Minimum Voltage cutoffs (V):',variable=str, menu=mSp1, default='2',tooltip='Minimum Voltage cutoff (V)')
sim_pars_capfrac=Parameter(name='sim_pars_capfrac', display_name='Fraction of full (dis)charge to simulate (only used for CC):',variable=str, menu=mSp1, default='0.8',tooltip= 'Fraction of full (dis)charge to simulate (only used for CC)\n[E.g. to go from cathode filling fraction 0.1 to 0.5 with cathode limiting capacity \n and Crate>0, set cs0_c = 0.1, capFrac=0.4]')

sim_pars_applv=Parameter(name='sim_pars_applv', display_name='Battery applied voltage (only used for CV) (V):',variable=str, menu=mSp1, default='0.12',tooltip='Battery applied voltage (only used for CV), V')
sim_pars_segments=Parameter(name='sim_pars_segments', display_name='CC/CV segments defining profile for CCsegments or CVsegments profile types :',variable=str, menu=mSp1, default='[(0.3, 0.4),(-0.5, 0.1)]',tooltip='CC/CV segments defining profile for profileType = CCsegments or CVsegments \n Segments are input as a list of tuples, with the first entry of each row being the setpoint (CC or CV) \n and the second being the time (in minutes) for which that setpoint is held)')
#Parameter(name='(', menu=mSp1, variable='label')
sim_pars_rampt=Parameter(name='sim_pars_rampt', display_name='Ramp time (s):',variable=str, menu=mSp1, default='1e-3',tooltip='Ramp time to go from near-equilibrium initial condition to setpoint for CC/CV,\n as a fraction of total simulation time, non-dimensional \n OR Ramp time to linearly transition to each new segment for CCsegments/CVsegments')
#Parameter(name='""', menu=mSp1, variable='label')
sim_condir=('false','absolute directory path')
sim_pars_contdir=Parameter(name='sim_pars_contdir', display_name='Continuation directory for the simulation:',variable=tuple, menu=mSp1, default=(sim_condir),tooltip='Continuation directory. If false, begin a fresh simulation with the\nspecified input parameters here. Otherwise, this should be the \nabsolute path to the output directory of the simulation to continue. \n If continuing, keep the same input parameters as the simulation to \n continue, changing only `Sim Params` subheadings (other than Nvol \n and Npart options).\nOptions: false, absolute directory path')
sim_pars_fintime=Parameter(name='sim_pars_fintime', display_name='Final time (only used for CV) (s):',variable=str, menu=mSp1, default='1.2e3',tooltip='Final time (only used for CV), [s]')
sim_pars_numdisc=Parameter(name='sim_pars_numdisc', display_name='Number disc. in time:',variable=str, menu=mSp1, default='200',tooltip='Number disc. in time\n Note time stepping is adaptive, variable-order stepping, so this\naffects only the interpolated output values, not simulation \n accuracy. The output will have all simulation values at a linear\nspacing between initial and final times with tsteps total outputs.')
#Parameter(name='[ Note time stepping is adaptive, variable-order stepping, so this affects\n only the interpolated output values, not simulation accuracy. The output will have all\n simulation values at a linear spacing between initial and final times with tsteps total\n outputs.]', menu=mSp1, variable='label')
sim_pars_reltol=Parameter(name='sim_pars_reltol', display_name='Relative Tolerance:',variable=str, menu=mSp2, default='1.0e-6',tooltip='Relative Tolerance')
sim_pars_abstol=Parameter(name='sim_pars_abstol', display_name='Absolute Tolerance:',variable=str, menu=mSp2, default='1.0e-6',tooltip='Absolute Tolerance')
sim_pars_temp=Parameter(name='sim_pars_temp', display_name='Temperature (K):',variable=str, menu=mSp2, default='298',tooltip='Temp, K\nWARNING -- temperature dependence not fully implemented. Use 298 K.')
Parameter(name='[WARNING -- temperature dependence not fully implemented. Use 298 K.]', menu=mSp2, variable='label')
TF=('true','false')
sim_pars_randseed=Parameter(name='sim_pars_randseed', display_name='Random Seed: ',variable=tuple, menu=mSp2, default=(TF),tooltip='Set to true to give a random seed in the simulation (affects noise, particle size distribution). \n Set to true exactly reproducible results -- useful for testing.')
sim_pars_seed=Parameter(name='sim_pars_seed', display_name='Value of the random seed, must be an integer:',variable=str, menu=mSp2, default='0',tooltip='Value of the random seed, must be an integer')
sim_pars_Rseries=Parameter(name='sim_pars_Rseries', display_name='Series resistance, (Ohm m^2):',variable=str, menu=mSp2, default='0',tooltip='Series resistance, [Ohm m^2]')
nvoltxt='- Nvol_c must be >= 1 \n - If Nvol_c = 1 & Nvol_a = 0 & Nvol_s = 0, simulate a single volume with no \n separator and infinitely fast counter electrode \n - If Nvol_a = 0, simulate a Li foil electrode'
sim_pars_cat_numdisc=Parameter(name='sim_pars_cat_numdisc', display_name='Cathode numer disc. in x direction (volumes in electrodes):',variable=str, menu=mSp2, default='2',tooltip=nvoltxt)
sim_pars_an_numdisc=Parameter(name='sim_pars_an_numdisc', display_name='Anode numer disc. in x direction (volumes in electrodes):',variable=str, menu=mSp2, default='2',tooltip=nvoltxt)
sim_pars_sep_numdisc=Parameter(name='sim_pars_sep_numdisc', display_name='Seperator numer disc. in x direction (volumes in electrodes):',variable=str, menu=mSp2, default='2',tooltip=nvoltxt)
#Parameter(name=, menu=mSp2, variable='label')
sim_pars_cat_nPartVol=Parameter(name='sim_pars_cat_nPartVol', display_name='Number of particles per volume for cathode:',variable=str, menu=mSp2, default='1',tooltip='Number of particles per volume for cathode')
sim_pars_an_nPartVol=Parameter(name='sim_pars_an_nPartVol', display_name='Number of particles per volume for anode:',variable=str, menu=mSp2, default='1',tooltip='Number of particles per volume for anode')

Parameter(name='', menu=mSp2, variable='label')			  

#simulation Electrode Parameters.
Parameter(name='', menu=mSoe, variable='label')
Parameter(name='Electrode Parameters', menu=mSoe, variable='label')
Parameter(name='', menu=mSoe, variable='label')

sim_elect_rc_Lifoil=Parameter(name='sim_elect_rc_Lifoil', display_name='Rate constant of the Li foil electrode( A/m^2) "Used only if Nvol_a = 0":',variable=str, menu=mSoe, default='1e0',tooltip='Rate constant of the Li foil electrode, A/m^2 \nUsed only if Nvol_a = 0')
sim_elect_Lifoil_filmres=Parameter(name='sim_elect_Lifoil_filmres', display_name='Film resistance on the Li foil (Ohm m^2):',variable=str, menu=mSoe, default='0e-0',tooltip='Film resistance on the Li foil, Ohm m^2')

Parameter(name='', menu=mSoe, variable='label')	

#simulation Particles Distribution
Parameter(name='', menu=mSpd, variable='label')
Parameter(name='Particle Distribution Parameters', menu=mSpd, variable='label')
Parameter(name='', menu=mSpd, variable='label')

sim_PartDist_cat_mean=Parameter(name='sim_PartDist_cat_mean', display_name='cathode particle size distribution info Mean (m):',variable=str, menu=mSpd, default='1e-6',tooltip='cathode particle size distribution info Mean (m)\nelectrode particle size distribution info, m\n C3 -- size along [100] direction \n sphere or cylinder -- radius \n If using stddev = 0, set Npart = 1 for that electrode to avoid \n wasting computational effort.')
sim_PartDist_cat_std=Parameter(name='sim_PartDist_cat_std', display_name='cathode particle size distribution info Standard Deviation (m):',variable=str, menu=mSpd, default='0',tooltip='cathode particle size distribution info Standard Deviation (m)\nelectrode particle size distribution info, m\n C3 -- size along [100] direction \n sphere or cylinder -- radius \n If using stddev = 0, set Npart = 1 for that electrode to avoid \n wasting computational effort.')
sim_PartDist_an_mean=Parameter(name='sim_PartDist_an_mean', display_name='anode particle size distribution info Mean (m):',variable=str, menu=mSpd, default='18e-6',tooltip='anode particle size distribution info Mean (m)\nelectrode particle size distribution info, m\n C3 -- size along [100] direction \n sphere or cylinder -- radius \n If using stddev = 0, set Npart = 1 for that electrode to avoid \n wasting computational effort.')
sim_PartDist_an_std=Parameter(name='sim_PartDist_an_std', display_name='anode particle size distribution info Standard Deviation (m):',variable=str, menu=mSpd, default='0',tooltip='anode particle size distribution info Mean (m)\nelectrode particle size distribution info, m\n C3 -- size along [100] direction \n sphere or cylinder -- radius \n If using stddev = 0, set Npart = 1 for that electrode to avoid \n wasting computational effort.')

sim_PartDist_cat_fillfrac=Parameter(name='sim_PartDist_cat_fillfrac', display_name='Cathone initial filling fraction:',variable=str, menu=mSpd, default='0.2',tooltip='Initial cathode filling fractions\n(for disch, anode starts full, cathode starts empty)')
sim_PartDist_an_fillfrac=Parameter(name='sim_PartDist_an_fillfrac', display_name='Anode initial filling fraction:',variable=str, menu=mSpd, default='0.495',tooltip='Initial anode filling fractions\n(for disch, anode starts full, cathode starts empty)')

Parameter(name='', menu=mSpd, variable='label')	

#simulation Conductivity
Parameter(name='', menu=mSc, variable='label')
Parameter(name='Conductivity Parameters', menu=mSc, variable='label')
Parameter(name='', menu=mSc, variable='label')
TF=('false','true')
sim_cond_cat_bulk_par=Parameter(name='sim_cond_cat_bulk_par', display_name='simulate bulk cathode conductivity (Ohms Law)?:',variable=tuple, menu=mSc, default=(TF),tooltip='Simulate bulk cathode conductivity (Ohms Law)?')
sim_cond_cat_bulk=Parameter(name='sim_cond_cat_bulk', display_name='Dimensional Conductivity of cathode (S/m) \n "Specify if simulate bulk conductivty = true" ',variable=str, menu=mSc, default='1e-1',tooltip='Dimensional conductivity (used if simBulkCond = true), S/m')
sim_cond_an_bulk_par=Parameter(name='sim_cond_an_bulk_par', display_name='simulate bulk anode conductivity (Ohms Law)?:',variable=tuple, menu=mSc, default=(TF),tooltip='Simulate bulk anode conductivity (Ohms Law)?')
sim_cond_an_bulk=Parameter(name='sim_cond_an_bulk', display_name='Dimensional Conductivity of anode (S/m)\n "Specify if simulate bulk conductivty = true" ',variable=str, menu=mSc, default='1e-1',tooltip='Dimensional conductivity (used if simBulkCond = true), S/m')

sim_cond_cat_connect_par=Parameter(name='sim_cond_cat_connect_par', display_name='simulate cathode particle connectivity losses (Ohms Law)?:',variable=tuple, menu=mSc, default=(TF),tooltip='simulate cathode particle connectivity losses (Ohms Law)?')
sim_cond_an_connect_par=Parameter(name='sim_cond_an_connect_par', display_name='simulate anode particle connectivity losses (Ohms Law)?:',variable=tuple, menu=mSc, default=(TF),tooltip='simulate anode particle connectivity losses (Ohms Law)?')

sim_condpart_cat_mean=Parameter(name='sim_condpart_cat_mean', display_name='Mean conductance between cathode particles(1/Ohm):',variable=str, menu=mSc, default='1e-14',tooltip='Mean Conductance between cathode particles, S = 1/Ohm')
sim_condpart_cat_std=Parameter(name='sim_condpart_cat_std', display_name='Standard deviation in conductance between cathode particles(1/Ohm):',variable=str, menu=mSc, default='0',tooltip='Standard deviation in Conductance between cathode particles, S = 1/Ohm')
sim_condpart_an_mean=Parameter(name='sim_condpart_an_mean', display_name='Mean conductance between anode particles(1/Ohm):',variable=str, menu=mSc, default='1e-14',tooltip='Mean Conductance between anode particles, S = 1/Ohm')
sim_condpart_an_std=Parameter(name='sim_condpart_an_std', display_name='Standard deviation in conductance between anode particles(1/Ohm):',variable=str, menu=mSc, default='0',tooltip='Standard deviation in Conductance between anode particles, S = 1/Ohm')

Parameter(name='', menu=mSc, variable='label')	

#simulation Geometry
Parameter(name='', menu=mSg, variable='label')
Parameter(name='Geometry Parameters', menu=mSg, variable='label')
Parameter(name='', menu=mSg, variable='label')

sim_geo_cat_thickness=Parameter(name='sim_geo_cat_thickness', display_name='Cathode thickness (m):',variable=str, menu=mSg, default='200e-6',tooltip='Cathode thickness')
sim_geo_an_thickness=Parameter(name='sim_geo_an_thickness', display_name='Anode thickness (m):',variable=str, menu=mSg, default='243e-6',tooltip='Anode thickness')
sim_geo_sep_thickness=Parameter(name='sim_geo_sep_thickness', display_name='Seperator thickness (m):',variable=str, menu=mSg, default='50e-6',tooltip='Seperator thickness')

sim_geo_cat_volload=Parameter(name='sim_geo_cat_volload', display_name='Solid volume fraction of active cathode material:',variable=str, menu=mSg, default='0.849',tooltip='Volume loading percents of active cathode material (volume fraction of solid\nthat is active material)')
sim_geo_an_volload=Parameter(name='sim_geo_an_volload', display_name='Solid volume fraction of active anode material:',variable=str, menu=mSg, default='0.956',tooltip='Volume loading percents of active anode material (volume fraction of solid\nthat is active material)')

sim_geo_cat_por=Parameter(name='sim_geo_cat_por', display_name='Cathode material porsity or \n Liquid volume fraction in cathode region:',variable=str, menu=mSg, default='0.3',tooltip='Porosities (liquid volume fraction in cathode region)')
sim_geo_an_por=Parameter(name='sim_geo_an_por', display_name='Anode material porsity or \n Liquid volume fraction in anode region:',variable=str, menu=mSg, default='0.3',tooltip='Porosities (liquid volume fraction in anode region)')
sim_geo_sep_por=Parameter(name='sim_geo_sep_por', display_name='Liquid volume fraction in seperator region',variable=str, menu=mSg, default='0.4',tooltip='Porosities (liquid volume fraction in seperator region)')

sim_geo_cat_Brug=Parameter(name='sim_geo_cat_Brug', display_name='Bruggeman exponent for cathode:',variable=str, menu=mSg, default='-0.5',tooltip='Bruggeman exponent cathode(tortuosity = porosity^bruggExp)')
sim_geo_an_Brug=Parameter(name='sim_geo_an_Brug', display_name='Bruggeman exponent for anode:',variable=str, menu=mSg, default='-0.5',tooltip='Bruggeman exponent anode(tortuosity = porosity^bruggExp)')
sim_geo_sep_Brug=Parameter(name='sim_geo_sep_Brug', display_name='Bruggeman exponent for seperator:',variable=str, menu=mSg, default='-0.5',tooltip='Bruggeman exponent seperator(tortuosity = porosity^bruggExp)')

Parameter(name='', menu=mSc, variable='label')	

#simulation Electrolyte Properties
Parameter(name='', menu=mSel, variable='label')
Parameter(name='Electrolyte properties', menu=mSel, variable='label')
Parameter(name='', menu=mSel, variable='label')

sim_elprop_InConc=Parameter(name='sim_elprop_InConc', display_name='Initial electrolyte concentration (mol/m^3):',variable=str, menu=mSel, default='1000',tooltip='Initial electrolyte conc., mol/m^3')
sim_elprop_cation_chargeNum=Parameter(name='sim_elprop_cation_chargeNum', display_name='Cation charge number:',variable=str, menu=mSel, default='1',tooltip='Cation/anion charge number (e.g. 2, -1 for CaCl_2)')
sim_elprop_anion_chargeNum=Parameter(name='sim_elprop_anion_chargeNum', display_name='Anion charge number:',variable=str, menu=mSel, default='-1',tooltip='Cation/anion charge number (e.g. 2, -1 for CaCl_2)')
sim_elprop_disNup=Parameter(name='sim_elprop_disNup', display_name='Dissociation number for cation:',variable=str, menu=mSel, default='1',tooltip='Cation/anion dissociation number (e.g. 1, 2 for CaCl_2)')
sim_elprop_disNum=Parameter(name='sim_elprop_disNum', display_name='Dissociation number for anion:',variable=str, menu=mSel, default='1',tooltip='Cation/anion dissociation number (e.g. 1, 2 for CaCl_2)')
Parameter(name=' WARNING: Using SM model with BV reaction models for electrode particles \n assumes individual ion activities are given by their concentrations \n for the reaction rate exchange current density.', menu=mSel, variable='label')
EL_MOD=('SM','dilute')
sim_elprop_model=Parameter(name='sim_elprop_model', display_name='Electrolyte Model:',variable=tuple, menu=mSel, default=(EL_MOD),tooltip='Electrolyte model,\nOptions: dilute, SM\ndilute: assume dilute, binary electrolyte model; phi in\nelectrolyte is an "inner potential" or an "quasi-electrostatic \npotential"\nSM: use Stefan-Maxwell model for electrolyte transport; phi in\nelectrolyte is that measured with a Li metal reference electrode\nrelative to some fixed position in solution.\nWARNING: Using SM model with BV reaction models for electrode\nparticles assumes individual ion activities are given by their\nconcentrations for the reaction rate exchange current density.')
SM_PROP=('LiClO4_PC','valoen_bernardi','test1')
sim_elprop_SMpropset=Parameter(name='sim_elprop_SMpropset', display_name='Stefan-Maxwell (SM) model property set:',variable=tuple, menu=mSel, default=(SM_PROP),tooltip='Stefan-Maxwell property set,\nOptions:\ntest1: parameter set for testing\nLiClO4_PC: electrolyte/solvent used in Fuller, Doyle, Newman 1994\nconductivity taken from dualfoil5.2.f\nvaloen_bernardi: LiPF6 in carbonates as in Bernardi and Go 2011')
sim_elprop_electrans=Parameter(name='sim_elprop_electrans', display_name='Number of electrons transfered in the reaction, "1 for Li/Li+":',variable=str, menu=mSel, default='1',tooltip='Reference electrode (defining the electrolyte potential) information:\nnumber of electrons transfered in the reaction, 1 for Li/Li+')
sim_elprop_Stoicof=Parameter(name='sim_elprop_Stoicof', display_name='Stoichiometric coeffecient of cation, "-1 for Li/Li+":',variable=str, menu=mSel, default='-1',tooltip='Stoichiometric coefficient of cation, -1 for Li/Li+')
Parameter(name='Dilute solution properties (used only if electrolyte Model Type is "dilute")', menu=mSel, variable='label')
sim_elprop_catdiff=Parameter(name='sim_elprop_catdiff', display_name='Dilute Solution properties cation diffusion (m^2/s) :',variable=str, menu=mSel, default='2.2e-10',tooltip='Dilute solution properties (used only if elyteModelType = "dilute")\nCation/anion diff, m^2/s \ne.g. for LiPF6 in EC/DMC, Dp = 2.2e-10, Dm = 2.94e-10')
sim_elprop_andiff=Parameter(name='sim_elprop_andiff', display_name='Dilute Solution properties anion diffusion (m^2/s):',variable=str, menu=mSel, default='2.94e-10',tooltip='Dilute solution properties (used only if elyteModelType = "dilute")\nCation/anion diff, m^2/s \ne.g. for LiPF6 in EC/DMC, Dp = 2.2e-10, Dm = 2.94e-10')

Parameter(name='', menu=mSel, variable='label')

#'MPET simulation launch
Parameter(name='', display_name='', variable='label',\
              menu=mSlaunch)

AND_tf=('false','true')

sim_elect_cat_overpar=Parameter(name='sim__elect_cat_overpar', display_name='Overwrite cathode Data with specified input file:', variable=tuple,menu=mSlaunch, default=(AND_tf),tooltip='Overwrite cathode Data configuration file True/False')
Cathode_config_file=Parameter(name='Cathode_config_file', display_name='Specify cathode parameters input files for overwiting cathode parameters:', variable='file',menu=mSlaunch, default=(""),tooltip='Overwrite the cathode parameters configuration file: "params_c.cfg"')

sim_elect_an_overpar=Parameter(name='sim__elect_an_overpar', display_name='Overwrite anode Data with specified input file:', variable=tuple,menu=mSlaunch, default=(AND_tf),tooltip='Overwrite anode Data configuration file True/False')
Anode_config_file=Parameter(name='Anode_config_file', display_name='Specify anode parameters input files for overwiting anode parameters:', variable='file',menu=mSlaunch, default=(""),tooltip='Overwrite the anode parameters configuration file: "params_a.cfg"')

sim_elect_sim_overpar=Parameter(name='sim__elect_sim_overpar', display_name='Overwrite simulation setup Data with specified input file:', variable=tuple,menu=mSlaunch, default=(AND_tf),tooltip='Overwrite simulation Data configuration file True/False')
System_config_file=Parameter(name='System_config_file', display_name='Specify simulation config parameters input files for overwiting simulation parameters:', variable='file',menu=mSlaunch, default=(""),tooltip='Overwrite the simulation parameters configuration file: "params_system.cfg"')

run_function=Parameter(name='   SOLVE   ', variable='function', default=RunMPETgui,\
				menu=mSlaunch,button_text='Run')
Parameter(name='', display_name='', variable='label',\
              menu=mSlaunch)

p.set_credits("mpetUI created by Surya Mitra \n  Purdue University ")

if __name__=="__main__":
	if p.is_gui_mode()==False:
		p.add_command(RunMPETtext)
	else:
		pass
	p()

