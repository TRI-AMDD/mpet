import ConfigParser

class mpetIO():

    def getConfig(self, paramfile="params.cfg"):
        P = ConfigParser.RawConfigParser()
        P.optionxform = str
        P.read(paramfile)
        return P

    def getDictFromConfig(self, P):
        D = {}

        # Simulation Parameters
        D["profileType"] = P.get('Sim Params', 'profileType')
        D["Crate"] = P.getfloat('Sim Params', 'Crate')
        D["Vset"] = P.getfloat('Sim Params', 'Vset')
        D["capFrac"] = P.getfloat('Sim Params', 'capFrac')
        D["tend"] = P.getfloat('Sim Params', 'tend')
        D["tsteps"] = P.getfloat('Sim Params', 'tsteps')
        D["Tabs"] = P.getfloat('Sim Params', 'T')
        D["Nvol_ac"] = [P.getint('Sim Params', 'Nvol_a'),
                        P.getint('Sim Params', 'Nvol_c')]
        D["Nvol_s"] = P.getint('Sim Params', 'Nvol_s')
        D["Npart_ac"] = [P.getint('Sim Params', 'Npart_a'),
                         P.getint('Sim Params', 'Npart_c')]

        # Particle info
        D["psd_mean_ac"] = [P.getfloat('Particles', 'mean_a'),
                            P.getfloat('Particles', 'mean_c')]
        D["psd_stddev_ac"] = [P.getfloat('Particles', 'stddev_a'),
                              P.getfloat('Particles', 'stddev_c')]
        D["cs0_ac"] = [P.getfloat('Particles', 'cs0_a'),
                       P.getfloat('Particles', 'cs0_c')]
        D["solidType_ac"] = [P.get('Particles', 'solidType_a'),
                             P.get('Particles', 'solidType_c')]
        D["solidDisc_ac"] = [P.getfloat('Particles', 'solidDisc_a'),
                             P.getfloat('Particles', 'solidDisc_c')]
        D["solidShape_ac"] = [P.get('Particles', 'solidShape_a'),
                              P.get('Particles', 'solidShape_c')]
        D["partThick_ac"] = [P.getfloat('Particles', 'partThick_a'),
                             P.getfloat('Particles', 'partThick_c')]

        # Conductivity
        D["simBulkCond_ac"] = [P.getboolean('Conductivity', 'simBulkCond_a'),
                               P.getboolean('Conductivity', 'simBulkCond_c')]
        D["mcond_ac"] = [P.getfloat('Conductivity', 'mcond_a'),
                         P.getfloat('Conductivity', 'mcond_c')]
        D["simPartCond_ac"] = [P.getboolean('Conductivity', 'simPartCond_a'),
                               P.getboolean('Conductivity', 'simPartCond_c')]
        D["G_mean_ac"] = [P.getfloat('Conductivity', 'G_mean_a'),
                          P.getfloat('Conductivity', 'G_mean_c')]
        D["G_stddev_ac"] = [P.getfloat('Conductivity', 'G_stddev_a'),
                            P.getfloat('Conductivity', 'G_stddev_c')]
        D["simSurfCond_ac"] = [P.getboolean('Conductivity', 'simSurfCond_a'),
                               P.getboolean('Conductivity', 'simSurfCond_c')]
        D["scond_ac"] = [P.getfloat('Conductivity', 'scond_a'),
                         P.getfloat('Conductivity', 'scond_c')]

        # Materials
        D["Omga_ac"] = [P.getfloat('Materials', 'Omega_a_a'),
                        P.getfloat('Materials', 'Omega_a_c')]
        D["kappa_ac"] = [P.getfloat('Materials', 'kappa_a'),
                         P.getfloat('Materials', 'kappa_c')]
        D["B_ac"] = [P.getfloat('Materials', 'B_a'),
                     P.getfloat('Materials', 'B_c')]
        D["rhos_ac"] = [P.getfloat('Materials', 'rhos_a'),
                        P.getfloat('Materials', 'rhos_c')]
        D["Vstd_ac"] = [P.getfloat('Materials', 'Vstd_a'),
                        P.getfloat('Materials', 'Vstd_c')]
        D["Dsld_ac"] = [P.getfloat('Materials', 'Dsld_a'),
                        P.getfloat('Materials', 'Dsld_c')]
        D["dgammasdc_ac"] = [P.getfloat('Materials', 'dgammasdc_a'),
                             P.getfloat('Materials', 'dgammasdc_c')]
        D["cwet_ac"] = [P.getfloat('Materials', 'cwet_a'),
                        P.getfloat('Materials', 'cwet_c')]
        D["delPhiEqFit_ac"] = [P.getboolean('Materials', 'delPhiEqFit_a'),
                          P.getboolean('Materials', 'delPhiEqFit_c')]
        D["material_ac"] = [P.get('Materials', 'material_a'),
                            P.get('Materials', 'material_c')]

        # Reactions
        D["rxnType_ac"] = [P.get('Reactions', 'rxnType_a'),
                           P.get('Reactions', 'rxnType_c')]
        D["k0_ac"] = [P.getfloat('Reactions', 'k0_a'),
                      P.getfloat('Reactions', 'k0_c')]
        D["alpha_ac"] = [P.getfloat('Reactions', 'alpha_a'),
                         P.getfloat('Reactions', 'alpha_c')]
        D["lambda_ac"] = [P.getfloat('Reactions', 'lambda_a'),
                          P.getfloat('Reactions', 'lambda_c')]

        # Geometry
        D["L_ac"] = [P.getfloat('Geometry', 'L_a'),
                     P.getfloat('Geometry', 'L_c')]
        D["L_s"] = P.getfloat('Geometry', 'L_s')
        D["P_L_ac"] = [P.getfloat('Geometry', 'P_L_a'),
                       P.getfloat('Geometry', 'P_L_c')]
        D["poros_ac"] = [P.getfloat('Geometry', 'poros_a'),
                         P.getfloat('Geometry', 'poros_c')]

        # Electrolyte
        D["c0"] = P.getfloat('Electrolyte', 'c0')
        D["zp"] = P.getfloat('Electrolyte', 'zp')
        D["zm"] = P.getfloat('Electrolyte', 'zm')
        D["Dp"] = P.getfloat('Electrolyte', 'Dp')
        D["Dm"] = P.getfloat('Electrolyte', 'Dm')

        # Constants
        D["k"] = P.getfloat('Constants', 'k')
        D["Tref"] = P.getfloat('Constants', 'Tref')
        D["e"] = P.getfloat('Constants', 'e')
        D["N_A"] = P.getfloat('Constants', 'N_A')

        return D

    def writeConfigFile(self, P, filename="output_params.cfg"):
        fo = open(filename, "w")
        P.write(fo)
        fo.close()
        return
