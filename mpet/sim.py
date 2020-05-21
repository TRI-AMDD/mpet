"""This module defines the actual simulation to be carried out.

In this, the model(s) are created and their variables are given initial conditions (if they are
differential variables) or initial guesses (if they never appear in equations in which they have
been differentiated in time).
"""
import os.path as osp

import daetools.pyDAE as dae
import numpy as np
import scipy.io as sio

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
            self.dataPrev = sio.loadmat(
                osp.join(ndD_s["prevDir"], "output_data.mat"))
            ndD_s["currPrev"] = self.dataPrev["current"][0,-1]
            ndD_s["phiPrev"] = self.dataPrev["phi_applied"][0,-1]

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
                        self.m.phi_bulk[l].SetInitialGuess(i, 0.0)
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
            
            #Cell potential initialization
            ndDVref=self.ndD_s["phiRef"]["c"]-self.ndD_s["phiRef"]["a"]
            if ndD_s['tramp']>0:
                phi_guess=0
            elif ndD_s['profileType']=='CV':
                phi_guess = self.ndD_s['Vset']
            elif ndD_s['profileType']=='CVsegments':
                phi_guess = self.ndD_s['segments'][0][0]
            else:
                phi_guess = ndDVref
            self.m.phi_applied.SetInitialGuess(phi_guess)

            #Initialize the ghost points used for boundary conditions
            if not self.m.SVsim:
                self.m.c_lyteGP_L.SetInitialGuess(ndD_s["c0"])
                self.m.phi_lyteGP_L.SetInitialGuess(0)

            #Separator electrolyte initialization
            for i in range(Nvol["s"]):
                self.m.c_lyte["s"].SetInitialCondition(i, ndD_s['c0'])
                self.m.phi_lyte["s"].SetInitialGuess(i, .5*(self.ndD_s["phiRef"]["c"]
                                                           +self.ndD_s["phiRef"]["a"]))
            
            #Anode and cathode electrolyte initialization
            for l in ndD_s["trodes"]:
                for i in range(Nvol[l]):
                    self.m.c_lyte[l].SetInitialCondition(i, ndD_s['c0'])
                    self.m.phi_lyte[l].SetInitialGuess(i, self.ndD_s["phiRef"][l])

        else:
            dPrev = self.dataPrev
            for l in ndD_s["trodes"]:
                self.m.ffrac[l].SetInitialGuess(
                    dPrev["ffrac_" + l][0,-1])
                for i in range(Nvol[l]):
                    self.m.R_Vp[l].SetInitialGuess(
                        i, dPrev["R_Vp_" + l][-1,i])
                    self.m.phi_bulk[l].SetInitialGuess(
                        i, dPrev["phi_bulk_" + l][-1,i])
                    for j in range(Npart[l]):
                        Nij = ndD_s["psd_num"][l][i,j]
                        part = self.m.particles[l][i,j]
                        solidType = self.ndD_e[l]["indvPart"][i, j]["type"]
                        partStr = "partTrode{l}vol{i}part{j}_".format(
                            l=l, i=i, j=j)
                        if solidType in ndD_s["1varTypes"]:
                            part.cbar.SetInitialGuess(
                                dPrev[partStr + "cbar"][0,-1])
                            for k in range(Nij):
                                part.c.SetInitialCondition(
                                    k, dPrev[partStr + "c"][-1,k])
                        elif solidType in ndD_s["2varTypes"]:
                            part.c1bar.SetInitialGuess(
                                dPrev[partStr + "c1bar"][0,-1])
                            part.c2bar.SetInitialGuess(
                                dPrev[partStr + "c2bar"][0,-1])
                            part.cbar.SetInitialGuess(
                                dPrev[partStr + "cbar"][0,-1])
                            for k in range(Nij):
                                part.c1.SetInitialCondition(
                                    k, dPrev[partStr + "c1"][-1,k])
                                part.c2.SetInitialCondition(
                                    k, dPrev[partStr + "c2"][-1,k])
            for i in range(Nvol["s"]):
                self.m.c_lyte["s"].SetInitialCondition(
                    i, dPrev["c_lyte_s"][-1,i])
                self.m.phi_lyte["s"].SetInitialGuess(
                    i, dPrev["phi_lyte_s"][-1,i])
            for l in ndD_s["trodes"]:
                for i in range(Nvol[l]):
                    self.m.c_lyte[l].SetInitialCondition(
                        i, dPrev["c_lyte_" + l][-1,i])
                    self.m.phi_lyte[l].SetInitialGuess(
                        i, dPrev["phi_lyte_" + l][-1,i])
            # Guess the initial cell voltage
            self.m.phi_applied.SetInitialGuess(
                dPrev["phi_applied"][0,-1])

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
                nextTime, dae.eStopAtModelDiscontinuity, True)
            self.ReportData(self.CurrentTime)
            self.Log.SetProgress(int(100. * self.CurrentTime/self.TimeHorizon))

            #Break when an end condition has been met
            if self.m.endCondition.npyValues:
                description=mod_cell.endConditions[int(self.m.endCondition.npyValues)]
                self.Log.Message("Ending condition: " + description, 0)
                break