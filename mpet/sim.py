"""This module defines the actual simulation to be carried out.

In this, the model(s) are created and their variables are given initial conditions (if they are
differential variables) or initial guesses (if they never appear in equations in which they have
been differentiated in time).
"""
import os.path as osp

import sympy as sym
import daetools.pyDAE as dae
import numpy as np
import scipy.io as sio
import h5py

import mpet.mod_cell as mod_cell
import mpet.daeVariableTypes
import mpet.utils as utils

class SimMPET(dae.daeSimulation):
    def __init__(self, ndD_s=None, ndD_e=None, tScale=None):
        dae.daeSimulation.__init__(self)
        if (ndD_s is None) or (ndD_e is None):
            raise Exception("Need input parameter dictionaries")
        self.ndD_s = ndD_s
        self.ndD_e = ndD_e
        self.tScale = tScale
        ndD_s["currPrev"] = 0.
        ndD_s["phiPrev"] = 0.
        if ndD_s["prevDir"] != "false":
            # Get the data mat file from prevDir
            self.dataPrev = osp.join(ndD_s["prevDir"], "output_data")
            data, f_type = utils.open_data_file(self.dataPrev)
            ndD_s["currPrev"] = utils.get_dict_key(data, "current", f_type, final = True)
            ndD_s["phiPrev"] = utils.get_dict_key(data, "phi_applied", f_type, final = True)
            if f_type == "h5py":
                #close file if it is a h5py file
                data.close()

        #Set absolute tolerances for variableTypes
        mpet.daeVariableTypes.mole_frac_t.AbsoluteTolerance=ndD_s["absTol"]
        mpet.daeVariableTypes.conc_t.AbsoluteTolerance=ndD_s["absTol"]
        mpet.daeVariableTypes.elec_pot_t.AbsoluteTolerance=ndD_s["absTol"]

        # Define the model we're going to simulate
        self.m = mod_cell.ModCell("mpet", ndD_s=ndD_s, ndD_e=ndD_e)

    def SetUpParametersAndDomains(self):
        # Domains
        ndD = self.ndD_s
        if ndD["Nvol"]["s"] >= 1:
            self.m.DmnCell["s"].CreateArray(ndD["Nvol"]["s"])    
        for l in ndD["trodes"]:
            self.m.DmnCell[l].CreateArray(ndD["Nvol"][l])
            self.m.DmnPart[l].CreateArray(ndD["Npart"][l])
            for i in range(ndD["Nvol"][l]):
                for j in range(ndD["Npart"][l]):
                    self.m.particles[l][i, j].Dmn.CreateArray(
                        int(ndD["psd_num"][l][i,j]))

    def SetUpVariables(self):
        ndD_s = self.ndD_s
        ndD_e = self.ndD_e
        Nvol = ndD_s["Nvol"]
        Npart = ndD_s["Npart"]
        phi_cathode = ndD_s["phi_cathode"]
        if ndD_s["prevDir"] == "false":
            # Solids
            for l in ndD_s["trodes"]:
                cs0 = self.ndD_s['cs0'][l]
                # Guess initial filling fractions
                self.m.ffrac[l].SetInitialGuess(cs0)
                for i in range(Nvol[l]):
                    # Guess initial volumetric reaction rates
                    self.m.R_Vp[l].SetInitialGuess(i, 0.0)
                    self.m.R_no_deg_Vp[l].SetInitialGuess(i, 0.0)
                    # Guess initial value for the potential of the
                    # electrodes
                    if l == "a":  # anode
                        self.m.phi_bulk[l].SetInitialGuess(i, self.ndD_s["phiRef"]["a"])
                    else:  # cathode
                        self.m.phi_bulk[l].SetInitialGuess(i, phi_cathode)
                    for j in range(Npart[l]):
                        Nij = ndD_s["psd_num"][l][i,j]
                        part = self.m.particles[l][i,j]
                        # Guess initial value for the average solid
                        # concentrations and set initial value for
                        # solid concentrations
                        solidType = self.ndD_e[l]["indvPart"][i,j]["type"]
                        if solidType in ndD_s["1varTypes"]:
                            part.cbar.SetInitialGuess(cs0)
                            for k in range(Nij):
                                part.c.SetInitialCondition(k, cs0)
                            #Initialize degradation variables
                            part.Rxn_SEI.SetInitialGuess(0)
                            part.c_solv.SetInitialGuess(ndD_s["c0_solv"])
                            part.a_e_SEI.SetInitialGuess(1)
                            part.dcSEIbardt.SetInitialGuess(0)
                            part.L1.SetInitialCondition(ndD_e[l]["indvPart"][i,j]["L10"])
                            part.L2.SetInitialCondition(ndD_e[l]["indvPart"][i,j]["L20"])
                        elif solidType in ndD_s["2varTypes"]:
                            part.c1bar.SetInitialGuess(cs0)
                            part.c2bar.SetInitialGuess(cs0)
                            part.cbar.SetInitialGuess(cs0)
                            epsrnd = 0.0001
                            rnd1 = epsrnd*(np.random.rand(Nij) - 0.5)
                            rnd2 = epsrnd*(np.random.rand(Nij) - 0.5)
                            rnd1 -= np.mean(rnd1)
                            rnd2 -= np.mean(rnd2)
                            for k in range(Nij):
                                part.c1.SetInitialCondition(k, cs0+rnd1[k])
                                part.c2.SetInitialCondition(k, cs0+rnd2[k])

            phi_guess = 0
            #Cell potential initialization
            if ndD_s['tramp'] > 0:
                phi_guess=0
            elif ndD_s['profileType']=='CV':
                if "t" not in str(ndD_s['Vset']):
                    #if it is a float set value, we set initial guess to be set point
                    phi_guess = self.ndD_s['Vset']
                else:
                    t = sym.Symbol("t")
                    #lambdifies waveform so that we can run with numpy functions
                    f = sym.lambdify(t, self.ndD_s['Vset'], modules = "numpy") 
                    phi_guess = f(0)
            else:
                phi_guess = 0
            self.m.phi_applied.SetInitialGuess(phi_guess)
            self.m.phi_cell.SetInitialGuess(phi_guess)


            #Initialize the ghost points used for boundary conditions
            if not self.m.SVsim:
                self.m.c_lyteGP_L.SetInitialGuess(ndD_s["c0"])
                self.m.phi_lyteGP_L.SetInitialGuess(0)

            #Separator electrolyte initialization
            for i in range(Nvol["s"]):
                self.m.c_lyte["s"].SetInitialCondition(i, ndD_s['c0'])
                self.m.phi_lyte["s"].SetInitialGuess(i, 0)
            
            #Anode and cathode electrolyte initialization
            for l in ndD_s["trodes"]:
                for i in range(Nvol[l]):
                    self.m.c_lyte[l].SetInitialCondition(i, ndD_s['c0'])
                    self.m.phi_lyte[l].SetInitialGuess(i, 0)

                    #Set electrolyte concentration in each particle
                    for j in range(Npart[l]):
                        self.m.particles[l][i,j].c_lyte.SetInitialGuess(ndD_s["c0"])
    
            for j in range(Npart[l]): #initialize reaction plane potential = electrolyte potential
                part = self.m.particles[l][i,j]
                solidType = self.ndD_e[l]["indvPart"][i,j]["type"]
 
            #set last values
            self.m.last_current.AssignValue(0)
            self.m.last_phi_applied.AssignValue(phi_guess)

            #tracks which cycle counter we're on by maccor increments
            self.m.maccor_cycle_counter.AssignValue(1)
            #track the maccor step number we're on
            self.m.maccor_step_number.SetInitialGuess(1)

        else:
            dPrev = self.dataPrev
            data, f_type = utils.open_data_file(dPrev)
            for l in ndD_s["trodes"]:
                self.m.ffrac[l].SetInitialGuess(
                    utils.get_dict_key(data, "ffrac_" + l, f_type, final = True))
                for i in range(Nvol[l]):
                    self.m.R_Vp[l].SetInitialGuess(
                        i, np.asscalar(data["R_Vp_" + l][-1,i]))
                    self.m.phi_bulk[l].SetInitialGuess(
                        i, np.asscalar(data["phi_bulk_" + l][-1,i]))
                    for j in range(Npart[l]):
                        Nij = ndD_s["psd_num"][l][i,j]
                        part = self.m.particles[l][i,j]
                        solidType = self.ndD_e[l]["indvPart"][i, j]["type"]
                        partStr = "partTrode{l}vol{i}part{j}_".format(
                            l=l, i=i, j=j)
                        
                        #Set the inlet port variables for each particle
                        part.c_lyte.SetInitialGuess(np.asscalar(data["c_lyte_" + l][-1,i]))
                        part.phi_lyte.SetInitialGuess(np.asscalar(data["phi_lyte_" + l][-1,i]))
                        part.phi_m.SetInitialGuess(np.asscalar(data["phi_bulk_" + l][-1,i]))


                        if solidType in ndD_s["1varTypes"]:
                            part.cbar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "cbar", f_type, final = True))
                            for k in range(Nij):
                                part.c.SetInitialCondition(
                                    k, np.asscalar(data[partStr + "c"][-1,k]))
                        elif solidType in ndD_s["2varTypes"]:
                            part.c1bar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "c1bar", f_type, final = True))
                            part.c2bar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "c2bar", f_type, final = True))
                            part.cbar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "cbar", f_type, final = True))
                            for k in range(Nij):
                                part.c1.SetInitialCondition(
                                    k, np.asscalar(data[partStr + "c1"][-1,k]))
                                part.c2.SetInitialCondition(
                                    k, np.asscalar(data[partStr + "c2"][-1,k]))
            for i in range(Nvol["s"]):
                self.m.c_lyte["s"].SetInitialCondition(
                    i, np.asscalar(data["c_lyte_s"][-1,i]))
                self.m.phi_lyte["s"].SetInitialGuess(
                    i, np.asscalar(data["phi_lyte_s"][-1,i]))
            for l in ndD_s["trodes"]:
                for i in range(Nvol[l]):
                    self.m.c_lyte[l].SetInitialCondition(
                        i, np.asscalar(data["c_lyte_" + l][-1,i]))
                    self.m.phi_lyte[l].SetInitialGuess(
                        i, np.asscalar(data["phi_lyte_" + l][-1,i]))
            
            #Read in the ghost point values
            if not self.m.SVsim:
                self.m.c_lyteGP_L.SetInitialGuess(utils.get_dict_key(data, "c_lyteGP_L", f_type, final = True))
                self.m.phi_lyteGP_L.SetInitialGuess(utils.get_dict_key(data, "phi_lyteGP_L", f_type, final = True))
            
            # Guess the initial cell voltage
            self.m.phi_applied.SetInitialGuess(utils.get_dict_key(data, "phi_applied", f_type, final = True))
            self.m.phi_cell.SetInitialGuess(utils.get_dict_key(data, "phi_cell", f_type, final = True))
   
            #set last values
            self.m.last_current.AssignValue(utils.get_dict_key(data, "current", f_type, final = True))
            self.m.last_phi_applied.AssignValue(utils.get_dict_key(data, "phi_applied", f_type, final = True))

            #tracks which cycle number we're on, using the cycle numbers tracked by maccor files
            self.m.maccor_cycle_counter.AssignValue(utils.get_dict_key(data, "maccor_cycle_counter", f_type, final = True))
            #track the maccor step number we're on
            self.m.maccor_step_number.SetInitialGuess(utils.get_dict_key(data,"maccor_step_number", f_type, final = True))

            if f_type == "h5py":
                #close file if it is a h5py file
                data.close()

        self.m.time_counter.AssignValue(0) #used to determine new time cutoffs at each section
        #used to determine if the mpet simulation has finished
        self.m.cycle_number.AssignValue(1)
        #The simulation runs when the endCondition is 0
        self.m.endCondition.AssignValue(0)

    def Run(self):
        """
        Overload the simulation "Run" function so that the simulation
        terminates when the specified condition is satisfied.
        """
        time = 0.
        tScale = self.tScale
        for nextTime in self.ReportingTimes:
            self.Log.Message(
                "Integrating from {t0:.2f} to {t1:.2f} s ...".format(
                    t0=self.CurrentTime*tScale, t1=nextTime*tScale), 0)
            time = self.IntegrateUntilTime(
                #nextTime, dae.eDoNotStopAtDiscontinuity, True)
                nextTime, dae.eStopAtModelDiscontinuity, True)
            self.ReportData(self.CurrentTime)
            self.Log.SetProgress(int(100. * self.CurrentTime/self.TimeHorizon))

            #Break when an end condition has been met
            if self.m.endCondition.npyValues:
                description=mod_cell.endConditions[int(self.m.endCondition.npyValues)]
                self.Log.Message("Ending condition: " + description, 0)
                break
