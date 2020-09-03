# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 14:52:36 2020

@authors: pasinger, dezhuang
"""

#%% import as needed
import numpy as np
import os
import sys
import configparser
import xmltodict

import mpet.io_utils as IO
import mpet.step_type_logic_functions as st

from mpet.utils import *

##########FUNCTIONS#############


def start_loop_process(step, loop_index, loop_counter, step_index):
    """Processes each begin do step of the JSON file. Starts a new set of steps to run
    and increments to next step, as well as initializes the new cycle counter.
    Inputs: step-current step we are at; loop_level is what level loop we are at
    loop_counter-an N*2 array keeping track of
    index we are at for each do loop as well as the total number of cycles to hit for
    loop_level: what number loop we are in
    each loop; step_index-in the JSON file, which step we are at.
    Updates loop_counter, step_list (new list of steps for simulator to run), and
    step_index"""
    #we assume next step for do step is just the next step index
    next_step_index = step_index + 1
    #loop_level incremenets by one from previous loop
    loop_level = 0 if len(loop_counter) == 0 else loop_counter[loop_index, 3] + 1
    #loop_level incremenets one because innermost loop increments by 1
    #int_ind is the internal index of the loop we just finished, it is set to loop_index by default (like if we need to add a new step)
    int_ind = loop_index #if no inner loops greater than it
    #we find the new loop_index since we've started looping through the beginning
    #if we are not in the last loop step
    if loop_counter.shape[0]-1 >= loop_index+1:
        #get all the loop steps that have a higher level than we are in
        mask1 = loop_counter[loop_index+1:,3] > loop_counter[loop_index,3]
        if np.any(mask1):
            #if not all loops full, then we need to find which loop we left off on
            #get the last index where they we are still in this inner loop cycle
            last_loop_index = get_end_of_loop(loop_counter, loop_index)
            #mask: if we have filled a certain loop step
            mask2 = (loop_counter[loop_index+1:last_loop_index,0] == loop_counter[loop_index+1:last_loop_index,1])

            if np.all(mask2):
                #if all loops are full, create a new one
                int_ind = loop_counter.shape[0]-1
            else:
                #if there are any internal loops for this code, get the first one for which [0] < [1] (that's not full)
                int_ind = np.argmax(~mask2) + loop_index
    #new loop index is incremented by 1 from internal index
    loop_index = int_ind + 1
    #prevents creating too many steps, creates new steps in loop_counter
    if loop_counter.shape[0] < loop_index + 1:
        #if we haven't created this loop before, append a loop counter step
        size_1 = loop_counter.shape[0] + 1
        #if any loop that pakcages it is not 0th array
        loop_counter = np.reshape(np.vstack((loop_counter, np.array([0, 0, next_step_index, loop_level]))), (size_1, 4))
    return loop_index, loop_counter, next_step_index


def finish_loop_process(step, loop_index, loop_counter):
    """Finishes loop by processing Loop step in a JSON file. Increments the loop
    step by +1 and sends back to first step if it is not the innnermost loop.
    If this loop is the innermost loop, then it increments all the way to the end
    value and runs a simulation for the number of innermost loop cycles.
    Inputs: step-current step; loop_level: what num loop by outermost to inner;
    loop_counter: tracker for loops"""
    #default cycle_number is set to 1 if not using battery cycler
    cycle_number = 1
    tot_cycles = step['Ends']['EndEntry']['Value'] #XXXhow to only get this for Loop Cnt?
    #tot_cycles = find_loop_limit(step)(step)
    loop_level = loop_counter[loop_index, 3]
    #sets total number of cycles in current loop
    loop_counter[loop_index, 1] = tot_cycles
    next_step_index = 0
    #if the loop level is the maximum of all the possible loop levels, it's in
    #the innermost loop
    if loop_counter[loop_index, 3] == np.amax(loop_counter[:, 3]):
        #if we are in the innermost loop, use internal battery cycler
        cycle_number = tot_cycles
        #set total number of cycles run to limit
        loop_counter[loop_index, 0] = loop_counter[loop_index, 1]
        #send to next loop
        loop_counter, loop_index, next_step_index = next_loop(loop_counter, loop_index, step)
    else:
        #loop counter for this loop is incremented by 1
        loop_counter[loop_index, 0] = loop_counter[loop_index, 0] + 1
        if loop_counter[loop_index, 0] == loop_counter[loop_index, 1]:
            #if we have reached the end of the loop
            loop_counter, loop_index, next_step_index = next_loop(loop_counter, loop_index, step)
        else:
            #if we haven't reached the counter limit, assign step to beginning of loop again
            next_step_index = int(loop_counter[loop_index, 2])
        #if we have any indexes that are in inner loops, we zero them out for the next counter
        #if our loop index is -1 however, we just leave it because it will be deleted
        if loop_index != -1: 
            ind_next_level = get_end_of_loop(loop_counter, loop_index) 
            #gets all internal loop indexes (level > current level until we hit next loop)
            loop_counter[loop_index+1:ind_next_level, 0] = 0 
    return cycle_number, loop_index, loop_counter, next_step_index


def process_basic_step(step, step_list, area, stepIndex):
    """Processes each step of the JSON file. Selects through step types
    and creates a state for each step.
    Inputs: step is current step, step_list is list of steps that need to be run.
    We should append the structure of current step to step_list
    Currently the structure of step list is a N*6 array.
    min/max are depending on charge/discharge
    [CC/CV set point, Vmin/Vmax, capfracmin/max, Cmin/Cmax, timemax, CC/CV type]
    CC/CV type: 1 is CC charge, 2 is CCdisch, 3 is CV charge, 4 is CV discharge
    """
    curr_step_properties, next_step_index = st.StepTypeLogic(step, area, stepIndex)(step, area, stepIndex)
    step_list.append(curr_step_properties)
    return step_list, next_step_index


def run_simulation(step_list, index_number, cycling_dict, cycle_number = 1):
    """Runs simulation for MPET, Liion, ets
    Inputs: step_list is the list of steps taken, simulation_type is the type of simulation:
    Currently MPET and LionSimba are supported. Cycle_number is the total number of cycles to be run,
    default is 1. returns empty step_list to restart"""
    if step_list != []:
        proc_step_list = [] #list of step segments
        
        #process step_list to make it readable by MPET
        for i in range(len(step_list)):
            if step_list[i].shape[0] == 1:
                proc_step_list.append(step_list[i].flatten().tolist())
                
            else:
                for j in range(len(step_list[i])):
                    proc_step_list.append(step_list[i][j].flatten().tolist())

        #runs simulation
        cycling_dict["step_" + str(index_number)] = {"segments": proc_step_list, "totalCycle": int(cycle_number)}
        #increments which set of simulations we are in
        index_number = index_number + 1
        #clears step list for next set of simulations
        step_list = []
    return step_list, index_number, cycling_dict

 
def get_cycling_dict(ndD_s, dD_s):
    """Inputs: ndD_s (uses the maccor file and processes it to return cycling_dict)
    which is a dictionary of steps to take in the cycler. Each step is a dictionary
    of format {'segments' (the step segments to take in this MPET run), and
    'totalCycle', the number of times to cycle this step}"""
    #process inputs
    maccor_file = ndD_s["profileType"]
    area = dD_s["active_area"]
    
    #%% Load the data and save variable of the TestSteps
    file_name = maccor_file
    if not os.path.isfile(maccor_file):
        file_name = ndD_s["paramfile_header"] + "/" + maccor_file
    
    with open(file_name) as f:
          data = xmltodict.parse(f.read(), process_namespaces=False, strip_whitespace=True)   

    f.close()

    Maccor = 'MaccorTestProcedure'
    proc = 'ProcSteps'
    test = 'TestStep'
    
    total_step_list = data[Maccor][proc][test]
    
    ##########
    
    total_cycle_counter = 0 #decided by placement of cycle counter
    
    #initializes step index
    step_index = 0 #what step index we are in in loop counter (first col of loop_counter)
    loop_counter = np.zeros((0, 4))
    #loop counter will be a N*4 array. first col is loop index, second col is total number of cycles in loop
    #third col will be step right after do step, fourth number will be what loop level it is in (-1 for not in loop, 0 for simple loop, 1 for nested loop in the loop above it, etc)
    loop_level = -1 # what level loop its on, no loop is -1, simple loop is 0, nested is 1...
    step_list = [] #keeps track of steps to run in the simulation
    index_number = 0 #we start with 0 simulations
    cycling_dict = {}
    
    #loops and processes the steps files
    
    #processes dict to add step_index as a key for easier state machine implementation
    max_step, total_step_list = get_dict_indexes(total_step_list)
    
    while step_index <= max_step-1:
        curr_step = total_step_list[step_index]
        #get the current step that we're on
        if curr_step["StepType"][:2] == "Do":
            #if we start tracking the loop, then we know that we need to run the previous process
            step_list, index_number, cycling_dict = run_simulation(step_list, index_number, cycling_dict) #mpet, lionsimba, or whatever
            #starts processing do loop
            loop_level, loop_counter, step_index = start_loop_process(curr_step, loop_level, loop_counter, step_index)
        elif curr_step["StepType"][:4] == "Loop": #finish tracking loop
            cycle_num, loop_level, loop_counter, step_index = finish_loop_process(curr_step, loop_level, loop_counter)
            step_list, index_number, cycling_dict = run_simulation(step_list, index_number, cycling_dict, cycle_number = cycle_num)#mpet, lionsimba, or whatever
        elif curr_step["StepType"] == "AdvCycle":
            total_cycle_counter += 1 #increments a step in total cycle counter
            step_index += 1
        elif curr_step["StepType"] == "End":
            step_list, index_number, cycling_dict = run_simulation(step_list, index_number, cycling_dict)
            step_index = 1e100 #steps should end
        else: #for a normal step, process as normal
            step_list, step_index = process_basic_step(curr_step, step_list, area, step_index)
    return cycling_dict
