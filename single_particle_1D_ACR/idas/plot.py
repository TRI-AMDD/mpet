import sys
import os
import time

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

def plot(plot_type, save_flag):
    infile = "acr_sp_0.001C.mat"
    infile = "acr_sp.mat"
    data = sio.loadmat(infile)
#    print data.keys()
    # Plot voltage profile
    if plot_type == "v":
        Vstd = 3.422        # Standard potential, V
        k = 1.381e-23       # Boltzmann constant
        T = 298             # Temp, K
        e = 1.602e-19       # Charge of proton, C
        fig = plt.figure()
#        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        ax = fig.add_subplot(111)
        times = data['Ray.phi_times'][0]
        voltage = Vstd - (k*T/e)*data['Ray.phi'][0]
        ax.plot(times, voltage)
        plt.show()
    # Plot concentration profiles
    elif plot_type == "c":
        if save_flag:
            outputdir = "{name}_image_output".format(name=infile[0:-4])
            os.mkdir(outputdir)
        else:
            # "interactive mode on"
            plt.ion()
        numtimes = len(data['Ray.c_times'][0])
        fig = plt.figure()
#        ax = fig.add_axes([0.08, 0.08, 0.85, 0.85])
        ax = fig.add_subplot(111)
        ax.set_ylim((0, 1))
        ax.set_xlabel('Particle Length [nm]')
        ax.set_ylabel('Concentration [fraction]')
        datay = data['Ray.c'][0]
        Lx = data['Ray.Lx'][0]
        numy = len(datay)
        datax = np.linspace(0, Lx, numy)
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datay)
#        for i in range(5):
        for i in range(numtimes):
#            print "img {num}/{tot}".format(num=i, tot=numtimes)
            datay = data['Ray.c'][i]
            line1.set_ydata(datay)
            if save_flag:
                filename = os.path.join(outputdir,
                        'c{num:06d}.png'.format(num=i))
                fig.savefig(filename)
            else:
                fig.canvas.draw()
                time.sleep(1e-3)
#            fig.clear()
    else:
        raise Exception("Unexpected plot type argument."
                + "Try 'v' or 'c' or 'd'")
    return

if __name__ == "__main__":
    if len(sys.argv) > 1:
        plot_type = sys.argv[1]
    else:
        plot_type = "v"
    save_flag = False
    plot(plot_type, save_flag)
