import sys
import os
import time

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

def plot(infile, plot_type, save_flag):
#    infile = "acr_sp_0.0C_2x1_40.mat"
#    infile = "acr_sp_0.0C.mat"
#    infile = "acr_sp_1C.mat"
#    infile = "acr_sp.mat"
    data = sio.loadmat(infile)
#    print 1./data['mpet.currset']
#    print data.keys()
#    print data['mpet.poros_sep']
#    print data['mpet.poros_trode']
#    print data['mpet.currset']
#    print data['mpet.tp']
#    print data['mpet.dim_Damb']
#    print data['mpet.c_lyte_sep']
#    print data['mpet.c_lyte_trode']
#    print data['mpet.phi_applied']
#    print data['mpet.phi_lyte_sep']
#    print data['mpet.phi_lyte_trode']
#    print data['mpet.solid_vol0_part0']
#    print data['mpet.solid_vol9_part0']
#    print data['mpet.j_plus']
#    print data['mpet.currset']
#    print np.sum(data['mpet.j_plus'][0])
#    print np.sum(data['mpet.j_plus'][1])
#    print data['mpet.particle_numVols']
    # Pick out some useful parameters
    Vstd = float(data['mpet.Vstd'][0][0])        # Standard potential, V
    k = float(data['mpet.k'][0][0])       # Boltzmann constant
    T = float(data['mpet.Tabs'][0][0])             # Temp, K
    e = float(data['mpet.e'][0][0])       # Charge of proton, C
    # Extract the reported simulation times
    times = data['mpet.phi_applied_times'][0]
    numtimes = len(times)
    Ntrode = int(data['mpet.NumTrode'][0][0])
    Nsep = int(data['mpet.NumSep'][0][0])
    numpart = int(data['mpet.NumPart'][0][0])

    # Plot voltage profile
    if plot_type == "v":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        voltage = Vstd - (k*T/e)*data['mpet.phi_applied'][0]
        ffvec = data['mpet.ffrac_cathode'][0]
        ax.plot(ffvec, voltage)
        ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        ax.set_ylabel("Voltage [V]")
        plt.show()

    # Plot electrolyte concentration
    elif plot_type == "elytec":
        if save_flag:
            outputdir = "{name}_image_output".format(name=infile[0:-4])
            os.mkdir(outputdir)
        else:
            # "interactive mode on"
            plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ymin = 0
        ymax = 1.5
        ax.set_xlabel('Battery Position [um]')
        ax.set_ylabel('Concentration of electrolyte [nondim]')
        sep = 'mpet.c_lyte_sep'
        trode = 'mpet.c_lyte_trode'
        datay_sep = data[sep][0]
        datay_trode = data[trode][0]
        Lsep = data['mpet.Lsep'][0] * 1e6
        Ltrode = data['mpet.Ltrode'][0] * 1e6
        datay = np.hstack((datay_sep, datay_trode))
        numy = len(datay)
        xmin = 0
        xmax = Lsep + Ltrode
        datax = np.linspace(xmin, xmax, numy)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay)
        ax.axvline(x=Lsep, ymin=ymin, ymax=ymax, linestyle='--', color='g')
        for i in range(numtimes):
            datay_sep = data[sep][i]
            datay_trode = data[trode][i]
            datay = np.hstack((datay_sep, datay_trode))
            line1.set_ydata(datay)
            if save_flag:
                filename = os.path.join(outputdir,
                        'elyte_c{num:06d}.png'.format(num=i))
                fig.savefig(filename)
            else:
                fig.canvas.draw()
                time.sleep(1e-3)

    # Plot electrolyte electrostatic potential
    elif plot_type == "elytep":
        if save_flag:
            outputdir = "{name}_image_output".format(name=infile[0:-4])
            os.mkdir(outputdir)
        else:
            # "interactive mode on"
            plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ymin = -5
        ymax = 5
        ax.set_xlabel('Battery Position [um]')
        ax.set_ylabel('Potential of electrolyte [nondim]')
        sep = 'mpet.phi_lyte_sep'
        trode = 'mpet.phi_lyte_trode'
        datay_sep = data[sep][0]
        datay_trode = data[trode][0]
        Lsep = data['mpet.Lsep'][0] * 1e6
        Ltrode = data['mpet.Ltrode'][0] * 1e6
        datay = np.hstack((datay_sep, datay_trode))
        numy = len(datay)
        xmin = 0
        xmax = Lsep + Ltrode
        datax = np.linspace(xmin, xmax, numy)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay)
        ax.axvline(x=Lsep, ymin=ymin, ymax=ymax, linestyle='--', color='g')
        for i in range(numtimes):
            datay_sep = data[sep][i]
            datay_trode = data[trode][i]
            datay = np.hstack((datay_sep, datay_trode))
            line1.set_ydata(datay)
            if save_flag:
                filename = os.path.join(outputdir,
                        'elyte_p{num:06d}.png'.format(num=i))
                fig.savefig(filename)
            else:
                fig.canvas.draw()
                time.sleep(1e-3)

    # Plot all solid concentrations
    elif plot_type == "csld":
        if save_flag:
            outputdir = "{name}_image_output".format(name=infile[0:-4])
            os.mkdir(outputdir)
        else:
            # "interactive mode on"
            plt.ion()
        fig, ax = plt.subplots(numpart, Ntrode, squeeze=False,
                sharey=True)
        sol = np.empty((numpart, Ntrode), dtype=object)
        lens = np.zeros((numpart, Ntrode))
        lines = np.empty((numpart, Ntrode), dtype=object)
        for i in range(numpart):
            for j in range(Ntrode):
                sol[i, j] = "mpet.solid_vol{j}_part{i}".format(i=i, j=j)
                lens[i, j] = data["mpet.psd_lengths"][0][j, i]
                # Remove axis ticks
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
                ax[i, j].yaxis.set_major_locator(plt.NullLocator())
                datay = data[sol[i, j]][0]
                numy = len(datay)
                datax = np.linspace(0, lens[i, j], numy)
                ax[i, j].set_ylim((0, 1))
                ax[i, j].set_xlim((0, lens[i, j]))
                line, = ax[i, j].plot(datax, datay)
                lines[i, j] = line
        for tval in range(numtimes):
            for i in range(numpart):
                for j in range(Ntrode):
                    datay = data[sol[i, j]][tval]
                    lines[i, j].set_ydata(datay)
            if save_flag:
                filename = os.path.join(outputdir,
                        'c{num:06d}.png'.format(num=tval))
                fig.savefig(filename)
            else:
                fig.canvas.draw()
                time.sleep(1e-3)

    # Plot average solid concentrations
    elif plot_type == "cbar":
        if save_flag:
            outputdir = "{name}_image_output".format(name=infile[0:-4])
            os.mkdir(outputdir)
        else:
            # "interactive mode on"
            plt.ion()
        # Get particle sizes (and max size) (length-based)
        psd_len = data['mpet.psd_lengths'][0]
        len_max = np.max(psd_len)
        len_min = np.min(psd_len)
        size_min = 0.10
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title("% = 00.0")
        tmin = np.min(times)
        tmax = np.max(times)
#        ax.patch.set_facecolor('gray')
        ax.patch.set_facecolor('white')
        ax.set_aspect('equal', 'box')
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
#        ax.invert_yaxis()
        rects = np.empty((Ntrode, numpart), dtype=object)
        cbar_mat = data['mpet.cbar_sld'][0]
        for (i,j),c in np.ndenumerate(cbar_mat):
            p_len = psd_len[i, j]
            size_frac = (p_len - len_min)/(len_max - len_min)
            size = ((size_frac) * (1-size_min) +
                    size_min)
            color = 'green'
            rects[i, j] = plt.Rectangle([i - size / 2, j - size / 2],
                    size, size, facecolor=color, edgecolor=color)
            ax.add_patch(rects[i, j])
        ax.autoscale_view()
        for tval in range(1,numtimes):
            t_current = times[tval]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ax.set_title("% = {perc:2.1f}".format(perc=tfrac))
            cbar_mat = data['mpet.cbar_sld'][tval]
            for (i,j),c in np.ndenumerate(cbar_mat):
                if c < 0.4:
                    color = 'green'
                elif c < 0.6:
                    color = 'yellow'
                else: # c > 0.6
                    color = 'red'
                rects[i, j].set_facecolor(color)
                rects[i, j].set_edgecolor(color)
            if save_flag:
                filename = os.path.join(outputdir,
                        'cbar{num:06d}.png'.format(num=tval))
                fig.savefig(filename)
            else:
                fig.canvas.draw()
                time.sleep(5e-2)

    # Plot cathode potential
    elif plot_type == "cathp":
        if save_flag:
            outputdir = "{name}_image_output".format(name=infile[0:-4])
            os.mkdir(outputdir)
        else:
            # "interactive mode on"
            plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ymin = -5
        ymax = 5
        ax.set_xlabel('Position in cathode [um]')
        ax.set_ylabel('Potential of cathode [nondim]')
        cathp = 'mpet.phi_cath'
        Ltrode = data['mpet.Ltrode'][0] * 1e6
        datay = data[cathp][0]
        numy = len(datay)
        xmin = 0
        xmax = Ltrode
        datax = np.linspace(xmin, xmax, numy)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay)
        for i in range(numtimes):
            datay = data[cathp][i]
            line1.set_ydata(datay)
            if save_flag:
                filename = os.path.join(outputdir,
                        'elyte_p{num:06d}.png'.format(num=i))
                fig.savefig(filename)
            else:
                fig.canvas.draw()
                time.sleep(1e-3)

    else:
        raise Exception("Unexpected plot type argument."
                + "Try 'v', 'elytec', 'elytep', 'cbar', 'csld'")
    return

if __name__ == "__main__":
    if len(sys.argv) > 2:
        plot_type = sys.argv[2]
    else:
        plot_type = "v"
    if len(sys.argv) < 2:
        raise Exception("Need input data file name")
    infile = sys.argv[1]
    if not os.path.isfile(os.path.join(os.getcwd(), infile)):
        raise Exception("Input file doesn't exist")
    save_flag = False
    plot(infile, plot_type, save_flag)
