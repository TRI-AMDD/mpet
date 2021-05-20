"""This module defines the actual simulation to be carried out.

In this, the model(s) are created and their variables are given initial conditions (if they are
differential variables) or initial guesses (if they never appear in equations in which they have
been differentiated in time).
"""
import os.path as osp

import daetools.pyDAE as dae
import numpy as np
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
            data = utils.open_data_file(self.dataPrev)
            ndD_s["currPrev"] = utils.get_dict_key(data, "current", final=True)
            ndD_s["phiPrev"] = utils.get_dict_key(data, "phi_applied", final=True)

            # close file if it is a h5py file
            if isinstance(data, h5py._hl.files.File):
                data.close()

        # Set absolute tolerances for variableTypes
        mpet.daeVariableTypes.mole_frac_t.AbsoluteTolerance = ndD_s["absTol"]
        mpet.daeVariableTypes.conc_t.AbsoluteTolerance = ndD_s["absTol"]
        mpet.daeVariableTypes.elec_pot_t.AbsoluteTolerance = ndD_s["absTol"]

        # Define the model we're going to simulate
        self.m = mod_cell.ModCell("mpet", ndD_s=ndD_s, ndD_e=ndD_e)

    def SetUpParametersAndDomains(self):
        # Domains
        ndD = self.ndD_s
        if ndD["Nvol"]["s"] >= 1:
            self.m.DmnCell["s"].CreateArray(ndD["Nvol"]["s"])
        for tr in ndD["trodes"]:
            self.m.DmnCell[tr].CreateArray(ndD["Nvol"][tr])
            self.m.DmnPart[tr].CreateArray(ndD["Npart"][tr])
            for i in range(ndD["Nvol"][tr]):
                for j in range(ndD["Npart"][tr]):
                    self.m.particles[tr][i, j].Dmn.CreateArray(
                        int(ndD["psd_num"][tr][i,j]))

    def SetUpVariables(self):
        ndD_s = self.ndD_s
        Nvol = ndD_s["Nvol"]
        Npart = ndD_s["Npart"]
        phi_cathode = ndD_s["phi_cathode"]
        if ndD_s["prevDir"] == "false":
            # Set initial temperature
            # Solids
            for tr in ndD_s["trodes"]:
                cs0 = self.ndD_s['cs0'][tr]
                # Guess initial filling fractions
                self.m.ffrac[tr].SetInitialGuess(cs0)
                for i in range(Nvol[tr]):
                    # Guess initial volumetric reaction rates
                    self.m.R_Vp[tr].SetInitialGuess(i, 0.0)
                    # set initial temperature condition
                    self.m.T_lyte[tr].SetInitialCondition(i, ndD_s["T0"])
                    # Guess initial value for the potential of the
                    # electrodes
                    if tr == "a":  # anode
                        self.m.phi_bulk[tr].SetInitialGuess(i, self.ndD_s["phiRef"]["a"])
                    else:  # cathode
                        self.m.phi_bulk[tr].SetInitialGuess(i, phi_cathode)
                    for j in range(Npart[tr]):
                        Nij = ndD_s["psd_num"][tr][i,j]
                        part = self.m.particles[tr][i,j]
                        # Guess initial value for the average solid
                        # concentrations and set initial value for
                        # solid concentrations
                        solidType = self.ndD_e[tr]["indvPart"][i,j]["type"]
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

            # Cell potential initialization
            if ndD_s['tramp'] > 0:
                phi_guess = 0
            elif ndD_s['profileType'] == 'CV':
                phi_guess = self.ndD_s['Vset']
            elif ndD_s['profileType'] == 'CVsegments':
                phi_guess = self.ndD_s['segments'][0][0]
            else:
                phi_guess = 0
            self.m.phi_applied.SetInitialGuess(phi_guess)
            self.m.phi_cell.SetInitialGuess(phi_guess)

            # Initialize the ghost points used for boundary conditions
            if not self.m.SVsim:
                self.m.c_lyteGP_L.SetInitialGuess(ndD_s["c0"])
                self.m.phi_lyteGP_L.SetInitialGuess(0)

            # Separator electrolyte initialization
            for i in range(Nvol["s"]):
                self.m.c_lyte["s"].SetInitialCondition(i, ndD_s['c0'])
                self.m.T_lyte["s"].SetInitialCondition(i, ndD_s['T0'])
                self.m.phi_lyte["s"].SetInitialGuess(i, 0)

            # Anode and cathode electrolyte initialization
            for tr in ndD_s["trodes"]:
                for i in range(Nvol[tr]):
                    self.m.c_lyte[tr].SetInitialCondition(i, ndD_s['c0'])
                    self.m.T_lyte[tr].SetInitialCondition(i, ndD_s['T0'])
                    self.m.phi_lyte[tr].SetInitialGuess(i, 0)

                    # Set electrolyte concentration in each particle
                    for j in range(Npart[tr]):
                        self.m.particles[tr][i,j].c_lyte.SetInitialGuess(ndD_s["c0"])
                        self.m.particles[tr][i,j].T_lyte.SetInitialGuess(ndD_s["T0"])

        else:
            dPrev = self.dataPrev
            data = utils.open_data_file(dPrev)
            for tr in ndD_s["trodes"]:
                self.m.ffrac[tr].SetInitialGuess(
                    utils.get_dict_key(data, "ffrac_" + tr, final=True))
                for i in range(Nvol[tr]):
                    self.m.R_Vp[tr].SetInitialGuess(
                        i, data["R_Vp_" + tr][-1,i])
                    self.m.phi_bulk[tr].SetInitialGuess(
                        i, data["phi_bulk_" + tr][-1,i])
                    for j in range(Npart[tr]):
                        Nij = ndD_s["psd_num"][tr][i,j]
                        part = self.m.particles[tr][i,j]
                        solidType = self.ndD_e[tr]["indvPart"][i, j]["type"]
                        partStr = "partTrode{l}vol{i}part{j}_".format(
                            l=tr, i=i, j=j)

                        # Set the inlet port variables for each particle
                        part.c_lyte.SetInitialGuess(data["c_lyte_" + tr][-1,i])
                        part.T_lyte.SetInitialGuess(data["T_lyte_" + tr][-1,i])
                        part.phi_lyte.SetInitialGuess(data["phi_lyte_" + tr][-1,i])
                        part.phi_m.SetInitialGuess(data["phi_bulk_" + tr][-1,i])

                        if solidType in ndD_s["1varTypes"]:
                            part.cbar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "cbar", final=True))
                            for k in range(Nij):
                                part.c.SetInitialCondition(
                                    k, data[partStr + "c"][-1,k])
                        elif solidType in ndD_s["2varTypes"]:
                            part.c1bar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "c1bar", final=True))
                            part.c2bar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "c2bar", final=True))
                            part.cbar.SetInitialGuess(
                                utils.get_dict_key(data, partStr + "cbar", final=True))
                            for k in range(Nij):
                                part.c1.SetInitialCondition(
                                    k, data[partStr + "c1"][-1,k])
                                part.c2.SetInitialCondition(
                                    k, data[partStr + "c2"][-1,k])
            for i in range(Nvol["s"]):
                self.m.c_lyte["s"].SetInitialCondition(
                    i, data["c_lyte_s"][-1,i])
                self.m.T_lyte["s"].SetInitialCondition(
                    i, data["T_lyte_s"][-1,i])
                self.m.phi_lyte["s"].SetInitialGuess(
                    i, data["phi_lyte_s"][-1,i])
            for tr in ndD_s["trodes"]:
                for i in range(Nvol[tr]):
                    self.m.c_lyte[tr].SetInitialCondition(
                        i, data["c_lyte_" + tr][-1,i])
                    self.m.T_lyte[tr].SetInitialCondition(
                        i, data["T_lyte_" + tr][-1,i])
                    self.m.phi_lyte[tr].SetInitialGuess(
                        i, data["phi_lyte_" + tr][-1,i])

            # Read in the ghost point values
            if not self.m.SVsim:
                self.m.c_lyteGP_L.SetInitialGuess(
                    utils.get_dict_key(data, "c_lyteGP_L", final=True))
                self.m.phi_lyteGP_L.SetInitialGuess(
                    utils.get_dict_key(data, "phi_lyteGP_L", final=True))

            # Guess the initial cell voltage
            self.m.phi_applied.SetInitialGuess(utils.get_dict_key(data, "phi_applied", final=True))
            self.m.phi_cell.SetInitialGuess(utils.get_dict_key(data, "phi_cell", final=True))

            # close file if it is a h5py file
            if isinstance(data, h5py._hl.files.File):
                data.close()

        # The simulation runs when the endCondition is 0
        self.m.endCondition.AssignValue(0)

    def Run(self):
        """
        Overload the simulation "Run" function so that the simulation
        terminates when the specified condition is satisfied.
        """
        tScale = self.tScale
        for nextTime in self.ReportingTimes:
            self.Log.Message(
                "Integrating from {t0:.2f} to {t1:.2f} s ...".format(
                    t0=self.CurrentTime*tScale, t1=nextTime*tScale), 0)
            self.IntegrateUntilTime(nextTime, dae.eStopAtModelDiscontinuity, True)
            self.ReportData(self.CurrentTime)
            self.Log.SetProgress(int(100. * self.CurrentTime/self.TimeHorizon))

            # Break when an end condition has been met
            if self.m.endCondition.npyValues:
                description = mod_cell.endConditions[int(self.m.endCondition.npyValues)]
                self.Log.Message("Ending condition: " + description, 0)
                break
