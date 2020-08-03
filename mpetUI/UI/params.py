### function to find and replace the particular parameter in the parameter file.

def findstrNrep(fkey,keyword,value):
	if fkey=='sim':
		file='params_system.cfg'
	elif fkey=='anode':
		file='params_a.cfg'
	elif fkey=='cathode':
		file='params_c.cfg'
	else:
		print('Wrong Parameter Name, Check parameter command name')
	f=open(file,'r')
	lines=f.readlines()
	f.close()
	#print(lines)
	keyword1=keyword+' '
	keyword2=keyword+'='
	r=0
	for rstr in lines:
		sub_index=rstr.find(keyword)
		sub_index1=rstr.find(keyword1)
		sub_index2=rstr.find(keyword2)
		#print(sub_index==0,',',sub_index1==0,',',sub_index2==0)
		if ((sub_index==0 and sub_index1==0) or (sub_index==0 and sub_index2==0)) :
			rf=r
			# print ("Line number:", r)
			# print("The position of 'contains' word: ", sub_index)
		else:
			r=r+1
	lines[rf]=keyword + ' = '+str(value)+'\n'
	f=open(file,'w')
	f.writelines(lines)
	f.close()

mpet_paramsUI_dict = { 
'anode_material_D':'D', 
'anode_material_Dfunc':'Dfunc', 
'anode_material_EvdW':'EvdW', 
'anode_material_Logpad':'logPad', 
'anode_material_cwet':'cwet', 
'anode_material_dgammadc':'dgammadc', 
'anode_material_grad_pen':'kappa', 
'anode_material_mag_noise':'noise_prefac', 
'anode_material_muRfunc':'muRfunc', 
'anode_material_noise':'noise', 
'anode_material_num_noise':'numnoise', 
'anode_material_oma':'Omega_a', 
'anode_material_omb':'Omega_b', 
'anode_material_omc':'Omega_c', 
'anode_material_rho_s':'rho_s', 
'anode_material_stress_coeff':'B', 
'anode_particle':'type', 
'anode_particle_thickness':'thickness', 
'anode_particle_discretization':'discretization', 
'anode_particle_shape':'shape', 
'anode_reaction_Rorg_energy':'lambda', 
'anode_reaction_charget':'alpha', 
'anode_reaction_ecd':'k0', 
'anode_reaction_filmres':'Rfilm', 
'anode_reaction_type':'rxnType', 
'cathode_material_D':'D', 
'cathode_material_Dfunc':'Dfunc', 
'cathode_material_EvdW':'EvdW', 
'cathode_material_Logpad':'logPad', 
'cathode_material_cwet':'cwet', 
'cathode_material_dgammadc':'dgammadc', 
'cathode_material_grad_pen':'kappa', 
'cathode_material_mag_noise':'noise_prefac', 
'cathode_material_muRfunc':'muRfunc', 
'cathode_material_noise':'noise', 
'cathode_material_num_noise':'numnoise', 
'cathode_material_oma':'Omega_a', 
'cathode_material_omb':'Omega_b', 
'cathode_material_omc':'Omega_c', 
'cathode_material_rho_s':'rho_s', 
'cathode_material_stress_coeff':'B', 
'cathode_particle':'type', 
'cathode_particle_thickness':'thickness', 
'cathode_particle_discretization':'discretization', 
'cathode_particle_shape':'shape', 
'cathode_reaction_Rorg_energy':'lambda', 
'cathode_reaction_charget':'alpha', 
'cathode_reaction_ecd':'k0', 
'cathode_reaction_filmres':'Rfilm', 
'cathode_reaction_type':'rxnType', 
'sim_PartDist_an_fillfrac':'cs0_a', 
'sim_PartDist_an_mean':'mean_a', 
'sim_PartDist_an_std':'stddev_a', 
'sim_PartDist_cat_fillfrac':'cs0_c', 
'sim_PartDist_cat_mean':'mean_c', 
'sim_PartDist_cat_std':'stddev_c', 
'sim_cond_an_bulk':'sigma_s_a', 
'sim_cond_an_bulk_par':'simBulkCond_a', 
'sim_cond_an_connect_par':'simPartCond_a', 
'sim_cond_cat_bulk':'sigma_s_c', 
'sim_cond_cat_bulk_par':'simBulkCond_c', 
'sim_cond_cat_connect_par':'simPartCond_c', 
'sim_condpart_an_mean':'G_mean_a', 
'sim_condpart_an_std':'G_stddev_a', 
'sim_condpart_cat_mean':'G_mean_c', 
'sim_condpart_cat_std':'G_stddev_c', 
'sim_elect_Lifoil_filmres':'Rfilm_foil', 
'sim_elect_rc_Lifoil':'k0_foil', 
'sim_elprop_InConc':'c0', 
'sim_elprop_SMpropset':'SMset', 
'sim_elprop_Stoicof':'sp', 
'sim_elprop_andiff':'Dm', 
'sim_elprop_anion_chargeNum':'zm', 
'sim_elprop_catdiff':'Dp', 
'sim_elprop_cation_chargeNum':'zp', 
'sim_elprop_disNum':'num', 
'sim_elprop_disNup':'nup', 
'sim_elprop_electrans':'n', 
'sim_elprop_model':'elyteModelType', 
'sim_geo_an_Brug':'BruggExp_a', 
'sim_geo_an_por':'poros_a', 
'sim_geo_an_thickness':'L_a', 
'sim_geo_an_volload':'P_L_a', 
'sim_geo_cat_Brug':'BruggExp_c', 
'sim_geo_cat_por':'poros_c', 
'sim_geo_cat_thickness':'L_c', 
'sim_geo_cat_volload':'P_L_c', 
'sim_geo_sep_Brug':'BruggExp_s', 
'sim_geo_sep_por':'poros_s', 
'sim_geo_sep_thickness':'L_s', 
'sim_pars_Rseries':'Rser', 
'sim_pars_Vcutmax':'Vmax', 
'sim_pars_Vcutmin':'Vmin', 
'sim_pars_abstol':'absTol', 
'sim_pars_an_nPartVol':'Npart_a', 
'sim_pars_an_numdisc':'Nvol_a', 
'sim_pars_applv':'Vset', 
'sim_pars_capfrac':'capFrac', 
'sim_pars_cat_nPartVol':'Npart_c', 
'sim_pars_cat_numdisc':'Nvol_c', 
'sim_pars_contdir':'prevDir', 
'sim_pars_crate':'Crate', 
'sim_pars_fintime':'tend', 
'sim_pars_numdisc':'tsteps', 
'sim_pars_prof':'profileType', 
'sim_pars_rampt':'tramp', 
'sim_pars_randseed':'randomSeed', 
'sim_pars_reltol':'relTol', 
'sim_pars_seed':'seed', 
'sim_pars_segments':'segments', 
'sim_pars_sep_numdisc':'Nvol_s', 
'sim_pars_temp':'T', 
}
