#!/usr/bin/env python3

from matplotlib import pyplot as plt
from VKML.gui.tk import Parameter, Parser, Menu
from pylab import *
import tkinter as tk
from tkinter import messagebox as tkMessageBox
import re, os, sys
import os.path
import subprocess
from UI import vizhelp
def RunMPETplt():
	par=plt_pars()
	#print (par)
	if (save_pars()=='Yes' and conv_out_txt()=='No'):
		subprocess.call(["mpetplot.py","sim_output",par,"save"])
	elif(save_pars()=='No' and conv_out_txt()=='No'):
		subprocess.call(["mpetplot.py","sim_output",par])
	elif (save_pars()=='Yes' and conv_out_txt()=='Yes'):
		subprocess.call(["mpetplot.py","sim_output","text"])
		subprocess.call([ "mpetplot.py", "sim_output",par,"save"])
	elif (save_pars()=='No' and conv_out_txt()=='Yes'):
		subprocess.call([" mpetplot.py","sim_output","text"])
		

#Create parser
p=Parser(title='MPET Visualization',help_text=vizhelp.Vizhelp)

#Create menus
#welcome = Menu(title='Welcome', parser=p)
min_files = Menu(title='MPET output', parser=p)


#PARAMETERS FOR WELCOME MENU

#Parameters for plotting as specified in 

PLT_PAR=('v','vt','curr','elytec','elytecf','elytep','elytepf','elytei','elyteif','surf_c','surf_a','soc_c','soc_a','csld_c','csld_a','cbarLine_c','cbarLine_a','cbar_full','cbar_a','cbar_c','bulkp_c','bulkp_a')
pltpar_text='"v" or "vt"       -- voltage vs filling fraction or vs time\n \n"curr" 	-- current vs time\n\n"elytec{f}"       -- electrolyte concentration (movie) or final snapshot with "f"\n\n"elytep{f}"       -- electrolyte potential (movie) or final snapshot with "f"\n\n"elytei{f}"       -- electrolyte current density (movie) or final snapshot with "f"\n\n"surf_{c,a}"      -- solid surface concentrations\n\n"soc_{c,a}"       -- overall utilization / state of charge of electrode\n\n"csld_{c,a}"      -- solid concentrations of particles in electrode (movie; used with "solidType_{c,a}" not homog)\n\n"cbarLine_{c,a}"  -- average concentration in each particle of electrode\n\n"cbar_{full,c,a}" -- average solid concentrations as changing colors (movie)\n\n"bulkp_{c,a}"     -- macroscopic electrode solid phase potential (movie)'
Parameter(name='', menu=min_files, variable='label')
Parameter(name='MPET Output Plots', menu=min_files, variable='label')
Parameter(name='', menu=min_files, variable='label')
plt_pars=Parameter(name='plt_pars', display_name='Choose parameters for plotting:', variable=tuple,\
              menu=min_files, default=(PLT_PAR),tooltip=pltpar_text)


Parameter(name='("full", "c", "a" indicate full cell, cathode, and anode)', menu=min_files, variable='label')
Parameter(name='', menu=min_files, variable='label')
Parameter(name='', display_name='', variable='label',\
              menu=min_files)
TF=('No','Yes')
save_pars=Parameter(name='save_pars', display_name='Save the plot or movie:', variable=tuple,menu=min_files, default=(TF),tooltip='Save the plot to a pdf file')
conv_out_txt=Parameter(name='conv_out_txt', display_name='convert the output to plain text (csv) format:', variable=tuple,menu=min_files, default=(TF),tooltip='Convert the output to plain text and csv format')
Parameter(name='', display_name='', variable='label',\
              menu=min_files)
run_function=Parameter(name='SOLVE', variable='function', default=RunMPETplt,\
				menu=min_files,tooltip='Run visualization')
Parameter(name='', display_name='', variable='label',\
              menu=min_files)
			  
		  
if __name__=="__main__":
	if p.is_gui_mode()==False:
		p.add_command(RunMPETplt)
	else:
		pass
	p()
