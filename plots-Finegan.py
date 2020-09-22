import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def read_csv(file,skip=0):
    data=np.genfromtxt(file,delimiter=',',deletechars='#',names=True,skip_header=skip)
    return(data)


def read_csv_no_headers(file,skip=0):
    data=np.genfromtxt(file,delimiter=',',deletechars='#',skip_header=skip)
    return(data)


try:
    mpetDir=sys.argv[1]
except:
    mpetDir="sim_output"

#Run the simulation
import mpet.main as main
main.main("configs/params_system_Finegan20.cfg")

#Process the simulation results
import mpet.plot.outmat2txt as outmat2txt
outmat2txt.main(mpetDir)

#Load cell data
mpet=read_csv(os.path.join(mpetDir,"generalData.txt"))
anode=read_csv(os.path.join(mpetDir,"cbarAnodeData.txt"),skip=5)

La=101                  #Thickness of anode in um
Nvol_a=len(anode.dtype) #Number of anode volumes
dx=La/Nvol_a            #Physical width of each cell

#Load particle concentrations
cs=[]
depth=[]
for vol in range(Nvol_a):
    filename="sldAnodeVol{0:03n}Part000ConcData.txt".format(vol)
    part=read_csv_no_headers(os.path.join(mpetDir,filename),4)
    cs.append(part[:,-1])
    depth.append(La-dx*(.5+vol))
cs=np.array(cs)
depth=np.array(depth)


#Plot 1, particle surface concentration vs depth at different times
###############################################################################
fig,ax=plt.subplots(1)

for time in [150,300,800]:
    time_ind = int(np.interp(time, mpet['Time_[s]'], range(len(mpet))))
    ax.plot(depth,cs[:,time_ind],label="{0:.0f} s".format(mpet['Time_[s]'][time_ind]))
    ax.set_xlabel("Depth (um)")
    ax.set_ylabel("Surface concentration")
    ax.set_ylim(0,1)

ax.legend()
plt.show()


#Plot #2, graphite SOC vs time ad different depths
###############################################################################
fig,ax=plt.subplots(1)

for x_pos in [1,30,60,90]:
    vol=int((La-x_pos)/dx)
    ax.plot(mpet['Time_[s]'],anode[str(vol)+"/0"],label="{pos:.0f} um".format(pos=x_pos))

#Plot the mean value
# ax.plot(mpet['Time_[s]'],mpet['Filling_fraction_of_anode'],label="mean")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Intercalation fraction")
ax.set_ylim(0,1)

ax.legend()
plt.show()