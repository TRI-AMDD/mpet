"""
Added module for cycling segments
"""

from daetools.pyDAE.variable_types import time_t
from pyUnits import s

import daetools.pyDAE as dae
import numpy as np

from mpet.daeVariableTypes import elec_pot_t


class CCCVCPcycle(dae.daeModel):

    def __init__(self, config, Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        self.config = config

        self.time_counter = dae.daeVariable(
            "time_counter", time_t, self, "restarts counter every time a new section is hit")
        self.last_current = dae.daeVariable(
            "last_current", dae.no_t, self,
            "tracks the current at the last step before a step is taken, used for ramp")
        self.last_phi_applied = dae.daeVariable(
            "last_phi_applied", elec_pot_t, self,
            "tracks the current at the last step before a step is taken, used for ramp")
        self.cycle_number = dae.daeVariable(
            "cycle_number", dae.no_t, self,
            "keeps track of which cycle number we are on in the mpet simulations")
        self.maccor_cycle_counter = dae.daeVariable(
            "maccor_cycle_counter", dae.no_t, self,
            "keeps track of which maccor cycle_number we are on")
        self.maccor_step_number = dae.daeVariable(
            "maccor_step_number", dae.no_t, self,
            "keeps track of which maccor step number we are on")

        # Get variables from the parent model
        self.current = Parent.current
        self.endCondition = Parent.endCondition
        self.phi_applied = Parent.phi_applied
        self.ffrac_limtrode = Parent.ffrac[config['limtrode']]

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)

        config = self.config

        limtrode = config["limtrode"]

        seg_array = np.array(config["segments"])
        constraints = seg_array[:,0]
        voltage_cutoffs = seg_array[:,1]
        capfrac_cutoffs = seg_array[:,2]
        crate_cutoffs = seg_array[:,3]
        time_cutoffs = seg_array[:,4]
        equation_type = seg_array[:,5]
        maccor_step_number = np.ones(seg_array.shape[0])

        if seg_array.shape[1] == 7:
            # if we also contain maccor step segments
            maccor_step_number = seg_array[:,6]

        ndDVref = config['c', 'phiRef']
        if 'a' in config['trodes']:
            ndDVref -= config['a', 'phiRef']

        # start state transition network
        self.stnCCCV = self.STN("CCCV")

        # start at a 0C state -1 for better initializing
        self.STATE("state_start")
        # if there is a ramp, add a ramp between this state and the last state
        if config["tramp"] > 0:
            self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
            eq = self.CreateEquation("Constraint_start")
            eq.Residual = self.current() + self.last_current() / \
                config["tramp"] * (dae.Time() - self.time_counter())/dae.Constant(1*s) - \
                self.last_current()
            self.ELSE()
            eq = self.CreateEquation("Constraint_start")
            eq.Residual = self.current()
            self.END_IF()
        else:
            eq = self.CreateEquation("Constraint_start")
            eq.Residual = self.current()

        # add new variable to assign maccor step number in equation
        self.maccor_step_number.AssignValue(1)

        # switch to next state unless cycle limit reached
        self.ON_CONDITION(self.cycle_number() <= dae.Constant(config['totalCycle']),
                          switchToStates=[('CCCV', 'state_0')],
                          setVariableValues=[(self.time_counter, dae.Time()),
                                             (self.last_current, self.current()),
                                             (self.last_phi_applied, self.phi_applied())])

        # loops through segments 1...N with indices 0...N-1
        # selects the correct charge/discharge type: 1 CCcharge, 2 CVcharge, 3 CCdisch, 4 CVdisch
        # if constraints is length of cycles, i is still the number in the nth cycle
        for i in range(0, len(constraints)):
            # creates new state

            self.STATE("state_" + str(i))
            new_state = "state_" + str(i+1)
            # if is CC charge, we set up equation and voltage cutoff
            # calculates time condition--if difference between current and prev start time is
            # larger than time cutof switch to next section
            # for our conditions for switching transition states, because we have multiple
            # different conditions. for switching transition states. if they are none the default
            # to switch is always false for that condition we use self.time_counter() < 0 as a
            # default false condition because daeTools does not accept False
            # if it is not None, we use a cutoff
            # if time is past the first cutoff, switch to nextstate
            time_cond = (self.time_counter() < dae.Constant(0*s)) if time_cutoffs[i] is None \
                else (dae.Time() - self.time_counter() >= dae.Constant(time_cutoffs[i]*s))

            # increment maccor cycle counter and switch to next state

            # add new variable to assign maccor step number in equation
            self.maccor_step_number.AssignValue(maccor_step_number[i])

            if equation_type[i] == 0:
                # if we are incrementing total_cycle
                if config["tramp"] > 0:
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() + self.last_current() / \
                        config["tramp"] * (dae.Time() - self.time_counter())/dae.Constant(1*s) - \
                        self.last_current()
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current()
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint" + str(i))
                    eq.Residual = self.current()

                # switch to next step immediately
                time_cond = (dae.Time() > self.time_counter())

                # here, we switch to next cycle if we hit the end of the previous cycle
                if i == len(constraints)-1:
                    # checks if the voltage, capacity fraction, or time segment conditions are
                    # broken
                    self.ON_CONDITION(time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.maccor_cycle_counter,
                                           self.maccor_cycle_counter() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                    # increases time_counter to increment to the beginning of the next segment

                else:
                    # checks if the voltage, capacity fraction, or time segment conditions are
                    # broken
                    self.ON_CONDITION(time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.maccor_cycle_counter,
                                           self.maccor_cycle_counter()+1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

            elif equation_type[i] == 1:

                # if not waveform input, set to constant value
                if config["tramp"] > 0:
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() - ((constraints[i] - self.last_current())
                                                    / config["tramp"]
                                                    * (dae.Time() - self.time_counter())
                                                    / dae.Constant(1*s) + self.last_current())
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() - constraints[i]
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() - constraints[i]

                # if hits voltage cutoff, switch to next state
                # if no voltage/capfrac cutoffs exist, automatically true, else is condition
                v_cond = (self.time_counter() < dae.Constant(0*s)) if voltage_cutoffs[i] is None \
                    else (-self.phi_applied() >= -voltage_cutoffs[i])
                # capacity fraction depends on cathode/anode. if charging, we cut off at
                # the capfrac
                cap_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if capfrac_cutoffs[i] is None \
                    else ((self.ffrac_limtrode() < 1 - capfrac_cutoffs[i]) if limtrode == "c"
                          else (self.ffrac_limtrode() > capfrac_cutoffs[i]))
                # for capacity condition, cathode is capped at 1-cap_frac, anode is at cap_Frac
                # if end state, then we send back to state 0 and also add one to cycle_number
                if i == len(constraints)-1:
                    # checks if the voltage, capacity fraction, or time segment conditions are
                    # broken
                    self.ON_CONDITION(v_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                    # increases time_counter to increment to the beginning of the next segment

                else:
                    # checks if the voltage, capacity fraction, or time segment conditions are
                    # broken
                    self.ON_CONDITION(v_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

            elif equation_type[i] == 2:

                if config["tramp"] > 0:
                    # if tramp, we use a ramp step to hit the value for better numerical
                    # stability
                    # if not waveform input, set to constant value
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.phi_applied() - \
                        ((constraints[i] - self.last_phi_applied())/config["tramp"] * (
                            dae.Time() - self.time_counter())/dae.Constant(1*s)
                            + self.last_phi_applied())
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.phi_applied() - constraints[i]
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.phi_applied() - constraints[i]

                # capacity fraction in battery is found by the filling fraction of the limiting
                # electrode
                # if is anode, will be capped at cap_frac, if is cathode needs to be capped at
                # 1-cap_frac
                # calc capacity and curr conditions (since Crate neg, we need to flip sign) to cut
                # off crate cutoff instead of voltage cutoff compared to CC
                crate_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if crate_cutoffs[i] is None \
                    else (-self.current() <= -crate_cutoffs[i])
                cap_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if capfrac_cutoffs[i] is None \
                    else ((self.ffrac_limtrode() <= 1 - capfrac_cutoffs[i]) if limtrode == "c"
                          else (self.ffrac_limtrode() >= capfrac_cutoffs[i]))
                # equation: constraining voltage
                # if past cap_frac, switch to next state
                # if end state, then we set endCondition to 3
                if i == len(constraints)-1:
                    # checks if crate, cap frac, or time segment conditions are broken
                    self.ON_CONDITION(crate_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                else:
                    # checks if crate, cap frac, or time segment conditions are broken
                    self.ON_CONDITION(crate_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

            elif equation_type[i] == 3:

                # constant power charge
                if config["tramp"] > 0:
                    # if tramp, we use a ramp step to hit the value for better numerical stability
                    # if not waveform input, set to constant value
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() * (self.phi_applied() + ndDVref) \
                        - ((constraints[i] - self.last_current()
                           * (self.last_phi_applied() + ndDVref))
                           / config["tramp"] * (dae.Time() - self.time_counter())
                           / dae.Constant(1 * s)
                           + self.last_current() * (self.last_phi_applied() + ndDVref))
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current()*(self.phi_applied() + ndDVref) - constraints[i]
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current()*(self.phi_applied() + ndDVref) - constraints[i]

                # if CC discharge, we set up capacity cutoff and voltage cutoff
                # needs to be minimized at capfrac for an anode and capped at 1-capfrac for a
                # cathode since discharging is delithiating anode and charging is lithiating anode
                cap_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if capfrac_cutoffs[i] is None \
                    else ((self.ffrac_limtrode() >= 1 - capfrac_cutoffs[i]) if limtrode == "c"
                          else (self.ffrac_limtrode() <= capfrac_cutoffs[i]))
                # voltage cutoff
                v_cond = (self.time_counter() < dae.Constant(0*s)) if voltage_cutoffs[i] is None \
                    else (-self.phi_applied() <= - voltage_cutoffs[i])
                crate_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if crate_cutoffs[i] is None \
                    else (-self.current() <= -crate_cutoffs[i])
                # if end state, then we set endCondition to 3
                if i == len(constraints)-1:
                    # if hits capacity fraction or voltage cutoff, switch to next state
                    self.ON_CONDITION(v_cond | cap_cond | crate_cond | time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                else:
                    # if hits capacity fraction or voltage cutoff, switch to next state
                    self.ON_CONDITION(v_cond | cap_cond | crate_cond | time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

            elif equation_type[i] == 4:

                # if not waveform input, set to constant value
                if config["tramp"] > 0:
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() - ((constraints[i] - self.last_current())
                                                    / config["tramp"]
                                                    * (dae.Time() - self.time_counter())
                                                    / dae.Constant(1*s) + self.last_current())
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() - constraints[i]
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current() - constraints[i]

                # if not waveform input, set to constant value
                # if CC discharge, we set up capacity cutoff and voltage cutoff
                # needs to be minimized at capfrac for an anode and capped at 1-capfrac for a
                # cathode since discharging is delithiating anode and charging is lithiating anode
                cap_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if capfrac_cutoffs[i] is None \
                    else ((self.ffrac_limtrode() > 1 - capfrac_cutoffs[i]) if limtrode == "c"
                          else (self.ffrac_limtrode() < capfrac_cutoffs[i]))
                # voltage cutoff
                v_cond = (self.time_counter() < dae.Constant(0*s)) if voltage_cutoffs[i] is None \
                    else (-self.phi_applied() <= -voltage_cutoffs[i])
                # if end state, then we set endCondition to 3
                if i == len(constraints)-1:
                    # if hits capacity fraction or voltage cutoff, switch to next state
                    self.ON_CONDITION(v_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                else:
                    # if hits capacity fraction or voltage cutoff, switch to next state
                    self.ON_CONDITION(v_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

            elif equation_type[i] == 5:

                # if CV discharge, we set up
                if config["tramp"] > 0:
                    # if tramp, we use a ramp step to hit the value for better numerical
                    # stability
                    # if not waveform input, set to constant value
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.phi_applied() - \
                        ((constraints[i] - self.last_phi_applied())/config["tramp"]
                         * (dae.Time() - self.time_counter())/dae.Constant(1*s)
                         + self.last_phi_applied())
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.phi_applied() - constraints[i]
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.phi_applied() - constraints[i]

                # conditions for cutting off: hits capacity fraction cutoff
                cap_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if capfrac_cutoffs[i] is None \
                    else ((self.ffrac_limtrode() >= 1 - capfrac_cutoffs[i]) if limtrode == "c"
                          else (self.ffrac_limtrode() <= capfrac_cutoffs[i]))
                # or hits crate limit
                crate_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if crate_cutoffs[i] is None \
                    else (self.current() <= crate_cutoffs[i])
                # if end state, then we set endCondition to 3
                if i == len(constraints)-1:
                    self.ON_CONDITION(crate_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                else:
                    self.ON_CONDITION(crate_cond | cap_cond | time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

            elif equation_type[i] == 6:

                if config["tramp"] > 0:
                    # if tramp, we use a ramp step to hit the value for better numerical stability
                    # if not waveform input, set to constant value
                    self.IF(dae.Time() < self.time_counter() + dae.Constant(config["tramp"]*s))
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current()*(self.phi_applied() + ndDVref) - \
                        ((constraints[i] - self.last_current()
                          * (self.last_phi_applied() + ndDVref))
                         / config["tramp"] * (dae.Time() - self.time_counter())
                         / dae.Constant(1*s) + self.last_current()
                         * (self.last_phi_applied() + ndDVref))
                    self.ELSE()
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current()*(self.phi_applied() + ndDVref) - constraints[i]
                    self.END_IF()
                else:
                    eq = self.CreateEquation("Constraint_" + str(i))
                    eq.Residual = self.current()*(self.phi_applied() + ndDVref) - constraints[i]

                # if CC discharge, we set up capacity cutoff and voltage cutoff
                # needs to be minimized at capfrac for an anode and capped at 1-capfrac for a
                # cathode
                # conditions for cutting off: hits capacity fraction cutoff
                cap_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if capfrac_cutoffs[i] is None \
                    else ((self.ffrac_limtrode() >= 1 - capfrac_cutoffs[i]) if limtrode == "c"
                          else (self.ffrac_limtrode() <= capfrac_cutoffs[i]))
                # or hits crate limit
                crate_cond = \
                    (self.time_counter() < dae.Constant(0*s)) if crate_cutoffs[i] is None \
                    else (self.current() <= crate_cutoffs[i])
                # voltage cutoff
                v_cond = (self.time_counter() < dae.Constant(0*s)) if voltage_cutoffs[i] is None \
                    else (-self.phi_applied() <= -voltage_cutoffs[i])
                # if end state, then we set endCondition to 3
                if i == len(constraints)-1:
                    # if hits capacity fraction or voltage cutoff, switch to next state
                    self.ON_CONDITION(v_cond | cap_cond | crate_cond | time_cond,
                                      switchToStates=[('CCCV', 'state_start')],
                                      setVariableValues=[
                                          (self.cycle_number, self.cycle_number() + 1),
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])
                else:
                    # if hits capacity fraction or voltage cutoff, switch to next state
                    self.ON_CONDITION(v_cond | cap_cond | crate_cond | time_cond,
                                      switchToStates=[('CCCV', new_state)],
                                      setVariableValues=[
                                          (self.time_counter, dae.Time()),
                                          (self.last_current, self.current()),
                                          (self.last_phi_applied, self.phi_applied())])

        self.END_STN()

        return
