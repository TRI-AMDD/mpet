import ConfigParser

class mpetIO():

    def readConfig(self, paramfile="params.cfg"):
        P = ConfigParser.RawConfigParser()
        P.optionxform = str
        P.read(paramfile)
        D = {}

        # Simulation Parameters
        D["profileType"] = P.get('Sim Params', 'profileType')
        D["dim_crate"] = P.getfloat('Sim Params', 'dim_crate')
        D["dim_Vset"] = P.getfloat('Sim Params', 'dim_Vset')
        D["cs0"] = P.getfloat('Sim Params', 'cs0')
        D["ffend"] = P.getfloat('Sim Params', 'ffend')
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
        D["Asep"] = P.getfloat('Geometry', 'Asep')
        D["Lp"] = P.getfloat('Geometry', 'Lp')
        D["poros"] = P.getfloat('Geometry', 'poros')

        # Electrolyte
        D["dim_c0"] = P.getfloat('Electrolyte Params', 'dim_c0')
        D["zp"] = P.getfloat('Electrolyte Params', 'zp')
        D["zm"] = P.getfloat('Electrolyte Params', 'zm')
        D["dim_Dp"] = P.getfloat('Electrolyte Params', 'dim_Dp')
        D["dim_Dm"] = P.getfloat('Electrolyte Params', 'dim_Dm')

        # Cathode Material Properties
        D["dim_a"] = P.getfloat('Cathode Material Props', 'dim_a')
        D["dim_kappa"] = P.getfloat('Cathode Material Props', 'dim_kappa')
        D["dim_b"] = P.getfloat('Cathode Material Props', 'dim_b')
        D["rhos"] = P.getfloat('Cathode Material Props', 'rhos')
        D["part_thick"] = P.getfloat('Cathode Material Props', 'part_thick')
        D["Vstd_c"] = P.getfloat('Cathode Material Props', 'Vstd')
        D["dim_mcond"] = P.getfloat('Cathode Material Props', 'dim_mcond')
        D["dim_scond"] = P.getfloat('Cathode Material Props', 'dim_scond')

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

#        # Simulation Parameters
#        profileType = P.get('Sim Params', 'profileType')
#        dim_crate = self.P.getfloat('Sim Params', 'dim_crate')
#        dim_Vset = self.P.getfloat('Sim Params', 'dim_Vset')
#        cs0 = self.P.getfloat('Sim Params', 'cs0')
#        ffend = P.getfloat('Sim Params', 'ffend')
#        mean = P.getfloat('Sim Params', 'mean')
#        stddev = P.getfloat('Sim Params', 'stddev')
#        Ntrode = P.getint('Sim Params', 'Ntrode')
#        numpart = P.getint('Sim Params', 'numpart')
#        tsteps = P.getfloat('Sim Params', 'tsteps')
#        Tabs = self.P.getfloat('Sim Params', 'T')
#        solidType = P.get('Sim Params', 'solidType')
#        solidShape = self.P.get('Sim Params', 'solidShape')
#        simBulkCathCond = self.P.getboolean('Sim Params',
#                'simBulkCathCond')
#        simSurfCathCond = self.P.getboolean('Sim Params',
#                'simSurfCathCond')
#
#        # Geometry
#        Ltrode = self.P.getfloat('Geometry', 'Ltrode')
#        Lsep = self.P.getfloat('Geometry', 'Lsep')
#        Asep = self.P.getfloat('Geometry', 'Asep')
#        Lp = self.P.getfloat('Geometry', 'Lp')
#        poros = self.P.getfloat('Geometry', 'poros')
#
#        # Electrolyte
#        c0 = self.P.getfloat('Electrolyte Params', 'c0')
#        zp = self.P.getfloat('Electrolyte Params', 'zp')
#        zm = self.P.getfloat('Electrolyte Params', 'zm')
#        Dp = self.P.getfloat('Electrolyte Params', 'Dp')
#        Dm = self.P.getfloat('Electrolyte Params', 'Dm')
#
#        # Cathode Material Properties
#        dim_a = self.P.getfloat('Cathode Material Props', 'dim_a')
#        dim_kappa = self.P.getfloat('Cathode Material Props', 'dim_kappa')
#        dim_b = self.P.getfloat('Cathode Material Props', 'dim_b')
#        rhos = self.P.getfloat('Cathode Material Props', 'rhos')
#        part_thick = self.P.getfloat('Cathode Material Props', 'part_thick')
#        Vstd_c = self.P.getfloat('Cathode Material Props', 'Vstd')
#        dim_mcond = self.P.getfloat('Cathode Material Props', 'dim_mcond')
#        dim_scond = self.P.getfloat('Cathode Material Props', 'dim_scond')
#
#        # Cathode reaction
#        rxnType_c = self.P.get('Cathode Reaction', 'rxnType')
#        dim_k0 = self.P.getfloat('Cathode Reaction', 'dim_k0')
#        alpha = self.P.getfloat('Cathode Reaction', 'alpha')
#        dim_lambda_c = self.P.getfloat('Cathode Reaction', 'dim_lambda_c')
#
#        # ACR info
#        cwet = self.P.getfloat('ACR info', 'cwet')
#        solid_disc = self.P.getfloat('ACR info', 'solid_disc')
#
#        # Constants
#        k = self.P.getfloat('Constants', 'k')
#        Tref = self.P.getfloat('Constants', 'Tref')
#        e = self.P.getfloat('Constants', 'e')
#        N_A = self.P.getfloat('Constants', 'N_A')

        return D, P

    def writeConfig(self, P, filename="output_params.cfg"):
        fo = open(filename, "w")
        P.write(fo)
        fo.close()
        return
