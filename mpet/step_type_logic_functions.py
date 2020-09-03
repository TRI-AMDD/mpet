# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 18:43:07 2020

@authors: pasinger, dezhuang
"""

import numpy as np
import re
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr

from mpet.utils import *

#curr_step_process:
#[CC/CV set point, Vmin/Vmax, capfracmin/max, Cmin/Cmax, timemax, CC/CV type]
#where Vmin/max Cmin/Cmax depends on whether it is charge or discharge

#Cases for StepType


def StepTypeLogic(step, area, stepIndex, totalCycleCounter):
    """Processes normal step types (rest, charge, discharge and chg_func).
    Inputs: step-current step we are at. Outputs: curr_step_process-
    a 1*6 or 2*6 array of the steps to tack onto the step list coming from this
    process. Each row of the array is [CC/CV set point, Vmin/Vmax,
    capfracmin/max, Cmin/Cmax, timemax, CC/CV type]; next_step_index-the next
    step index of which step the step will go to after breaking.
    totalCycleCounter is the total number of cycles that we have by the maccor cycler"""
    StepType = step['StepType']
    switcher = {
            'Rest': case_Rest,
            'Charge': case_Charge,
            'Chg Func': case_ChgFunc,
            'Dischrge': case_Dischrge,
            'FastWave': case_Waveform
            }
    func = switcher.get(StepType, lambda: "invalid StepType")
    #return curr_step_process, which is #segments for the first 6 arrays, 7th array is maccor step number, 8th array
    #is plus for cycling index
    return func


def case_Rest(step, area, stepIndex, totalCycleCounter):
    """Processes rest steps as CC = 0 steps. Inputs and outputs same as
    SwitchTypeLogic."""
    #gets end entry time
    #+1 to stepIndex because it is started at 0
    curr_step_process = np.reshape(np.array([0, None, None, None, None, 1, stepIndex+1, totalCycleCounter]), (1, 8))
    #assign CC = 0 as set value and CC charge
    curr_step_process, next_step_index = process_ends(step["Ends"], curr_step_process, area, charge_type = 1) 
    return curr_step_process, next_step_index


def case_Charge(step, area, stepIndex, totalCycleCounter): # add Ends and Limits for current/voltage
    StepMode = step['StepMode']
    StepValue = step['StepValue']
    #negative C rates because charge
    curr_step_process = np.reshape(np.array([0, None, None, None, None, 0, stepIndex+1, totalCycleCounter]), (1, 8))
    #initialize guess
    if StepMode == 'Current':
       curr_step_process[0][5] = 1 #specify CC charge
       curr_step_process[0][0] = process_current(StepValue, 1, area) #set CC
       if step['Limits'] != None:
           if 'Voltage' in step['Limits'].keys():
                #CCCV step
                curr_step_process = np.vstack((curr_step_process, np.array([0, None, None, None, None, 0, stepIndex+1, totalCycleCounter])))
                #if first step CC, cutoff is Vlim
                curr_step_process[0][1] = float(step['Limits']['Voltage'])
                curr_step_process[1][5] = 2 #CV charge
                curr_step_process[1][0] = curr_step_process[0][1] #set Vset
           else:
                print("CC + CC step functionality not handled")
    elif StepMode == 'Voltage':
        #CV step
        curr_step_process[0][5] = 2 #specify CV charge
        curr_step_process[0][0] = float(StepValue) #set CV
        if step['Limits'] != None: 
            if 'Current' in step['Limits'].keys():
                #CVCC step
                curr_step_process = np.vstack((curr_step_process, np.array([0, None, None, None, None, 0, stepIndex+1, totalCycleCounter])))
                #if first step CV, cutofff is CC, negative because of charge
                curr_step_process[0][3] = process_current(step['Limits']['Current'], 1, area)
                curr_step_process[1][5] = 1 #CC charge
                curr_step_process[1][0] = curr_step_process[0][3] #set Cset
            else:
                print("CV + CV step functionality not handled")
 
    else:
        print('invalid StepMode in StepType=Charge')

    curr_step_process, next_step_index = process_ends(step['Ends'], curr_step_process, area, charge_type = 1)
    #processes step ends
    return curr_step_process, next_step_index


def case_ChgFunc(step, area, stepIndex, totalCycleCounter):
    #assigning temporarily
    StepValue = step['StepValue']
    y0 = float(StepValue.split('|')[0])
    text = StepValue.split('|')[1]
    #get symbols to exchange CRATE and STIME for
    symbol_map = {
                  'CRATE': '1',
                  'STIME': 't'
                    }
    to_symbols = re.compile('|'.join(re.escape(key) for key in symbol_map.keys())) 
    # run through the text looking for keys (regex) and replacing them with the values from the dict
    text = to_symbols.sub(lambda x: symbol_map[x.group()], text) 
    C = sym.Symbol('C')
    t = sym.Symbol('t')
    f = parse_expr(text) + y0
    chg_dichg = 0
    if f.subs(t, 0.001) >= 0: # if charge
        chg_dichg = 1
    else:
        chg_dichg = -1
    #assign CC = 0 as set value and CC charge
    curr_step_process = np.reshape(np.array([str(f), None, None, None, None, chg_dichg, stepIndex+1, totalCycleCounter]), (1, 8))
    #assign CC = 0 as set value and CC charge
    curr_step_process, next_step_index = process_ends(step["Ends"], curr_step_process, area, totalCycleCounter, charge_type = 1) 
    return curr_step_process, next_step_index
    # Do we intend to use this type of step? I thought we would use waveform instead.


def case_Dischrge(step, area, stepIndex, totalCycleCounter): # need to add in duration,EndType cases
    StepMode = step['StepMode']
    StepValue = step['StepValue']
    #we assume only the first end entry value has meaning
    curr_step_process = np.reshape(np.array([0, None, None, None, None, 0, stepIndex+1, totalCycleCounter]), (1, 8))

    if StepMode == 'Current':
        #CC step
        curr_step_process[0][5] = 4 #specify CC discharge
        curr_step_process[0][0] = process_current(StepValue, 0, area) #set CC
        if step['Limits'] != None:
            if 'Voltage' in step['Limits'].keys():
                #CCCV step
                curr_step_process = np.vstack((curr_step_process, np.array([0, None, None, None, None, 0, stepIndex+1, totalCycleCounter])))
                #if first step CC, cutoff is Vlim
                curr_step_process[0][1] = float(step['Limits']['Voltage'])
                curr_step_process[1][5] = 5 #CV discharge
                curr_step_process[1][0] = curr_step_process[0][1] #set Vset
            else:
                print("CC + CC step functionality not handled")
    elif StepMode == 'Voltage':
        #CV step
        curr_step_process[0][5] = 5 #specify CV charge
        curr_step_process[0][0] = float(StepValue) #set CV
        if step['Limits'] != None:
            if 'Current' in step['Limits'].keys():
                #CVCCstep
                curr_step_process = np.vstack((curr_step_process, np.array([0, None, None, None, None, 0, stepIndex+1, totalCycleCounter])))
                #if first step CV, cutofff is CC, negative because of charge
                curr_step_process[0][3] = process_current(step['Limits']['Current'], 0, area)
                curr_step_process[1][5] = 4 #CC charge
                curr_step_process[1][0] = curr_step_process[0][3] #set Cset
            else:
                print("CV + CV step functionality not handled")
    else:
        print('invalid StepMode in StepType=DisCharge')
   
    curr_step_process, next_step_index = process_ends(step['Ends'], curr_step_process, area, charge_type = 0)
    #processes step ends
    return curr_step_process, next_step_index



def case_Waveform(step, area, stepIndex, totalCycleCounter): # need to add in duration,EndType cases
    #read in waveform file
    wavefile = step['StepValue']
    data = np.loadtxt(wavefile + '.MWF', dtype={'names': ('charge_discharge', \
        'control_mode', 'control_value', 'time(s)', 'end_type', 'ineq', \
        'end_value', 'X', 'Y', 'Z'), 'formats': ('S1', 'S1', 'f4', 'f4', 'S1',\
        'S2', 'f4', 'S1', 'f4', 'S1')})
    #get waveform simulation process
    curr_step_process = process_waveform_segment(data, area, totalCycleCounter)
    next_step_index = stepIndex + 1
    #next step index is incrememented from current step
    return curr_step_process, next_step_index

