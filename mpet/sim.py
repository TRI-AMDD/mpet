"""This module defines the actual simulation to be carried out.

In this, the model(s) are created and their variables are given initial conditions (if they are
differential variables) or initial guesses (if they never appear in equations in which they have
been differentiated in time).
"""
import os.path as osp

import sympy as sym
import daetools.pyDAE as dae
import numpy as np
import h5py

import mpet.mod_cell as mod_cell
import mpet.daeVariableTypes


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
            #self.dataPrev = h5py.File(osp.join(ndD_s["prevDir"], "output_data.hdf5"), 'r')
            self.dataPrev = osp.join(ndD_s["prevDir"], "output_data.hdf5")
            with h5py.File(self.dataPrev, 'r') as hf:
                ndD_s["currPrev"] = np.asscalar(hf["current"][-1])
                ndD_s["phiPrev"] = np.asscalar(hf["phi_applied"][-1])

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
         
            #set last values
            self.m.last_current.AssignValue(0)
            self.m.last_phi_applied.AssignValue(phi_guess)

            #tracks which cycle counter we're on by maccor increments
            self.m.maccor_cycle_counter.AssignValue(1)
            #track the maccor step number we're on
            self.m.maccor_step_number.SetInitialGuess(1)

        else:
            dPrev = self.dataPrev
            with h5py.File(dPrev, 'r') as hf:
                for l in ndD_s["trodes"]:
                    self.m.ffrac[l].SetInitialGuess(
                        np.asscalar(hf["ffrac_" + l][-1]))
                    for i in range(Nvol[l]):
                        self.m.R_Vp[l].SetInitialGuess(
                            i, np.asscalar(hf["R_Vp_" + l][-1,i]))
                        self.m.phi_bulk[l].SetInitialGuess(
                            i, np.asscalar(hf["phi_bulk_" + l][-1,i]))
                        for j in range(Npart[l]):
                            Nij = ndD_s["psd_num"][l][i,j]
                            part = self.m.particles[l][i,j]
                            solidType = self.ndD_e[l]["indvPart"][i, j]["type"]
                            partStr = "partTrode{l}vol{i}part{j}_".format(
                                l=l, i=i, j=j)
                            
                            #Set the inlet port variables for each particle
                            part.c_lyte.SetInitialGuess(np.asscalar(hf["c_lyte_" + l][-1,i]))
                            part.phi_lyte.SetInitialGuess(np.asscalar(hf["phi_lyte_" + l][-1,i]))
                            part.phi_m.SetInitialGuess(np.asscalar(hf["phi_bulk_" + l][-1,i]))


                            if solidType in ndD_s["1varTypes"]:
                                part.cbar.SetInitialGuess(
                                    np.asscalar(hf[partStr + "cbar"][-1]))
                                for k in range(Nij):
                                    part.c.SetInitialCondition(
                                        k, np.asscalar(hf[partStr + "c"][-1,k]))
                            elif solidType in ndD_s["2varTypes"]:
                                part.c1bar.SetInitialGuess(
                                    np.asscalar(hf[partStr + "c1bar"][-1]))
                                part.c2bar.SetInitialGuess(
                                    np.asscalar(hf[partStr + "c2bar"][-1]))
                                part.cbar.SetInitialGuess(
                                    np.asscalar(hf[partStr + "cbar"][-1]))
                                for k in range(Nij):
                                    part.c1.SetInitialCondition(
                                        k, np.asscalar(hf[partStr + "c1"][-1,k]))
                                    part.c2.SetInitialCondition(
                                        k, np.asscalar(hf[partStr + "c2"][-1,k]))
                for i in range(Nvol["s"]):
                    self.m.c_lyte["s"].SetInitialCondition(
                        i, np.asscalar(hf["c_lyte_s"][-1,i]))
                    self.m.phi_lyte["s"].SetInitialGuess(
                        i, np.asscalar(hf["phi_lyte_s"][-1,i]))
                for l in ndD_s["trodes"]:
                    for i in range(Nvol[l]):
                        self.m.c_lyte[l].SetInitialCondition(
                            i, np.asscalar(hf["c_lyte_" + l][-1,i]))
                        self.m.phi_lyte[l].SetInitialGuess(
                            i, np.asscalar(hf["phi_lyte_" + l][-1,i]))
                
                #Read in the ghost point values
                if not self.m.SVsim:
                    self.m.c_lyteGP_L.SetInitialGuess(np.asscalar(hf["c_lyteGP_L"][-1]))
                    self.m.phi_lyteGP_L.SetInitialGuess(np.asscalar(hf["phi_lyteGP_L"][-1]))
                
                # Guess the initial cell voltage
                self.m.phi_applied.SetInitialGuess(np.asscalar(hf["phi_applied"][-1]))
                self.m.phi_cell.SetInitialGuess(np.asscalar(hf["phi_cell"][-1]))

                #set last values
                self.m.last_current.AssignValue(np.asscalar(hf["current"][-1]))
                self.m.last_phi_applied.AssignValue(np.asscalar(hf["phi_applied"][-1]))

                #tracks which cycle number we're on, using the cycle numbers tracked by maccor files
                self.m.maccor_cycle_counter.AssignValue(np.asscalar(hf["maccor_cycle_counter"][-1]))
                #track the maccor step number we're on
                self.m.maccor_step_number.SetInitialGuess(np.asscalar(hf["maccor_step_number"][-1]))


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
                nextTime, dae.eDoNotStopAtDiscontinuity, True)
            self.ReportData(self.CurrentTime)
            self.Log.SetProgress(int(100. * self.CurrentTime/self.TimeHorizon))

            #Break when an end condition has been met
            if self.m.endCondition.npyValues:
                description=mod_cell.endConditions[int(self.m.endCondition.npyValues)]
                self.Log.Message("Ending condition: " + description, 0)
                break
