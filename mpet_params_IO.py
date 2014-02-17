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
        D["dim_crate"] = P.getfloat('Sim Params', 'dim_crate')
        D["dim_Vset"] = P.getfloat('Sim Params', 'dim_Vset')
        D["cs0"] = P.getfloat('Sim Params', 'cs0')
        D["ffend"] = P.getfloat('Sim Params', 'ffend')
        D["tend"] = P.getfloat('Sim Params', 'tend')
        D["mean"] = P.getfloat('Sim Params', 'mean')
        D["stddev"] = P.getfloat('Sim Params', 'stddev')
        D["Ntrode"] = P.getint('Sim Params', 'Ntrode')
        D["numpart"] = P.getint('Sim Params', 'numpart')
        D["tsteps"] = P.getfloat('Sim Params', 'tsteps')
        D["Tabs"] = P.getfloat('Sim Params', 'T')
        D["solidType"] = P.get('Sim Params', 'solidType')
        D["solidShape"] = P.get('Sim Params', 'solidShape')
        D["simBulkCathCond"] = P.getboolean('Sim Params',
                'simBulkCathCond')
        D["simSurfCathCond"] = P.getboolean('Sim Params',
                'simSurfCathCond')

        # Geometry
        D["Ltrode"] = P.getfloat('Geometry', 'Ltrode')
        D["Lsep"] = P.getfloat('Geometry', 'Lsep')
        D["Lp"] = P.getfloat('Geometry', 'Lp')
        D["poros"] = P.getfloat('Geometry', 'poros')

        # Electrolyte
        D["dim_c0"] = P.getfloat('Electrolyte Params', 'dim_c0')
        D["zp"] = P.getfloat('Electrolyte Params', 'zp')
        D["zm"] = P.getfloat('Electrolyte Params', 'zm')
        D["dim_Dp"] = P.getfloat('Electrolyte Params', 'dim_Dp')
        D["dim_Dm"] = P.getfloat('Electrolyte Params', 'dim_Dm')

        # Cathode Material Properties
        D["Omega_a"] = P.getfloat('Cathode Material Props', 'Omega_a')
        D["dim_kappa"] = P.getfloat('Cathode Material Props', 'dim_kappa')
        D["dim_b"] = P.getfloat('Cathode Material Props', 'dim_b')
        D["rhos"] = P.getfloat('Cathode Material Props', 'rhos')
        D["part_thick"] = P.getfloat('Cathode Material Props', 'part_thick')
        D["Vstd_c"] = P.getfloat('Cathode Material Props', 'Vstd')
        D["dim_mcond"] = P.getfloat('Cathode Material Props', 'dim_mcond')
        D["dim_scond"] = P.getfloat('Cathode Material Props', 'dim_scond')
        D["Dsld_c"] = P.getfloat('Cathode Material Props', 'Dsld')
        D["etaFit"] = P.getboolean('Cathode Material Props', 'etaFit')
        D["material_c"] = P.get('Cathode Material Props', 'material')

        # Cathode reaction
        D["rxnType_c"] = P.get('Cathode Reaction', 'rxnType')
        D["dim_k0"] = P.getfloat('Cathode Reaction', 'dim_k0')
        D["alpha"] = P.getfloat('Cathode Reaction', 'alpha')
        D["dim_lambda_c"] = P.getfloat('Cathode Reaction', 'dim_lambda_c')

        # ACR info
        D["cwet"] = P.getfloat('ACR info', 'cwet')
        D["solid_disc"] = P.getfloat('ACR info', 'solid_disc')

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
