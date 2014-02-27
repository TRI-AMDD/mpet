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
        D["Nvol_c"] = P.getint('Sim Params', 'Nvol_c')
        D["Nvol_s"] = P.getint('Sim Params', 'Nvol_s')
        D["Nvol_a"] = P.getint('Sim Params', 'Nvol_a')
        D["Npart_c"] = P.getint('Sim Params', 'Npart_c')
        D["Npart_a"] = P.getint('Sim Params', 'Npart_a')
        D["z"] = P.getfloat('Sim Params', 'z')

        # Particle info
        D["mean_c"] = P.getfloat('Particles', 'mean_c')
        D["stddev_c"] = P.getfloat('Particles', 'stddev_c')
        D["mean_a"] = P.getfloat('Particles', 'mean_a')
        D["stddev_a"] = P.getfloat('Particles', 'stddev_a')
        D["cs0_c"] = P.getfloat('Particles', 'cs0_c')
        D["cs0_a"] = P.getfloat('Particles', 'cs0_a')
        D["solidType_c"] = P.get('Particles', 'solidType_c')
        D["solidType_a"] = P.get('Particles', 'solidType_a')
        D["solidDisc_c"] = P.getfloat('Particles', 'solidDisc_c')
        D["solidDisc_a"] = P.getfloat('Particles', 'solidDisc_a')
        D["cwet_c"] = P.getfloat('Particles', 'cwet_c')
        D["cwet_a"] = P.getfloat('Particles', 'cwet_a')
        D["solidShape_c"] = P.get('Particles', 'solidShape_c')
        D["solidShape_a"] = P.get('Particles', 'solidShape_a')
        D["partThick_c"] = P.getfloat('Particles', 'partThick_c')
        D["partThick_a"] = P.getfloat('Particles', 'partThick_a')

        # Conductivity
        D["simBulkCond_c"] = P.getboolean('Conductivity', 'simBulkCond_c')
        D["simBulkCond_a"] = P.getboolean('Conductivity', 'simBulkCond_a')
        D["mcond_c"] = P.getfloat('Conductivity', 'mcond_c')
        D["mcond_a"] = P.getfloat('Conductivity', 'mcond_a')
        D["simSurfCond_c"] = P.getboolean('Conductivity', 'simSurCathCond_c')
        D["simSurfCond_a"] = P.getboolean('Conductivity', 'simSurCathCond_a')
        D["scond_c"] = P.getfloat('Conductivity', 'scond_c')
        D["scond_a"] = P.getfloat('Conductivity', 'scond_a')

        # Materials
        D["Omega_a_c"] = P.getfloat('Materials', 'Omega_a_c')
        D["Omega_a_a"] = P.getfloat('Materials', 'Omega_a_a')
        D["kappa_c"] = P.getfloat('Materials', 'kappa_c')
        D["kappa_a"] = P.getfloat('Materials', 'kappa_a')
        D["B_c"] = P.getfloat('Materials', 'B_c')
        D["B_a"] = P.getfloat('Materials', 'B_a')
        D["rhos_c"] = P.getfloat('Materials', 'rhos_c')
        D["rhos_a"] = P.getfloat('Materials', 'rhos_a')
        D["Vstd_c"] = P.getfloat('Materials', 'Vstd_c')
        D["Vstd_a"] = P.getfloat('Materials', 'Vstd_a')
        D["Dsld_c"] = P.getfloat('Materials', 'Dsld_c')
        D["Dsld_a"] = P.getfloat('Materials', 'Dsld_a')
        D["etaFit_c"] = P.getboolean('Materials', 'etaFit_c')
        D["etaFit_a"] = P.getboolean('Materials', 'etaFit_a')
        D["material_c"] = P.get('Materials', 'material_c')
        D["material_a"] = P.get('Materials', 'material_a')

        # Reactions
        D["rxnType_c"] = P.get('Reactions', 'rxnType_c')
        D["rxnType_a"] = P.get('Reactions', 'rxnType_a')
        D["k0_c"] = P.getfloat('Reactions', 'k0_c')
        D["k0_a"] = P.getfloat('Reactions', 'k0_a')
        D["alpha_c"] = P.getfloat('Reactions', 'alpha_c')
        D["alpha_a"] = P.getfloat('Reactions', 'alpha_a')
        D["lambda_c"] = P.getfloat('Reactions', 'lambda_c')
        D["lambda_a"] = P.getfloat('Reactions', 'lambda_a')

        # Geometry
        D["L_c"] = P.getfloat('Geometry', 'L_c')
        D["L_a"] = P.getfloat('Geometry', 'L_a')
        D["L_s"] = P.getfloat('Geometry', 'L_s')
        D["PL_c"] = P.getfloat('Geometry', 'Lp_c')
        D["PL_a"] = P.getfloat('Geometry', 'Lp_a')
        D["poros_c"] = P.getfloat('Geometry', 'poros_c')
        D["poros_a"] = P.getfloat('Geometry', 'poros_a')

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
