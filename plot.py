import sys
import os
import time

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

def plot(infile, plot_type, save_flag):
    pfx = 'mpet.'
    data = sio.loadmat(infile)
    # Pick out some useful parameters
    Vstd = float(data[pfx + 'Vstd'][0][0])        # Standard potential, V
    k = float(data[pfx + 'k'][0][0])       # Boltzmann constant
    T = float(data[pfx + 'Tabs'][0][0])             # Temp, K
    e = float(data[pfx + 'e'][0][0])       # Charge of proton, C
    # Discretization info
    Ntrode = int(data[pfx + 'NumTrode'][0][0])
    Nsep = int(data[pfx + 'NumSep'][0][0])
    numpart = int(data[pfx + 'NumPart'][0][0])
    # Extract the reported simulation times
    times = data[pfx + 'phi_applied_times'][0]
    numtimes = len(times)
    # Simulation type
    sim_ACR = data[pfx + "type_ACR"][0][0]
    sim_homog = data[pfx + "type_homog"][0][0]
    shape_sphere = data[pfx + "shape_sphere"][0][0]
    shape_C3 = data[pfx + "shape_C3"][0][0]

    # Print relevant simulation info
#    if plot_type == "params":
    if True:
        if sim_ACR:
            print "Type: ACR"
        elif sim_homog:
            print "Type: homog"
        else:
            print "Type: ?"
        if shape_sphere:
            print "Shape: sphere"
        elif shape_C3:
            print "Shape: C3"
        else:
            print "Shape: ?"
        print "C_rate:", data[pfx + "C_rate"][0][0]
        print "psd_mean [nm]:", data[pfx + "psd_mean"][0][0]*1e9
        print "psd_stddev [nm]:", data[pfx + "psd_stddev"][0][0]*1e9
        print "Nsep:", Nsep
        print "Ntrode:", Ntrode
        print "Npart:", numpart
        print "dim_Dp [m^2/s]:", data[pfx + "dim_Dp"][0][0]
        print "dim_Dm [m^2/s]:", data[pfx + "dim_Dm"][0][0]
        print "dim_Damb [m^2/s]:", data[pfx + "dim_Damb"][0][0]
        print "alpha:", data[pfx + "alpha"][0][0]
        print "dim_k0 [A/m^2]:", data[pfx + "dim_k0"][0][0]
        mcond_bool = data[pfx + "cath_bulk_cond"][0][0]
        if mcond_bool:
            print ("cathode conductivity loss: Yes -- " +
                    "dim_mcond [S/m]: " +
                    str(data[pfx + "dim_mcond"][0][0]))
        else:
            print "cathode conductivity loss: No"

    # Plot voltage profile
    if plot_type == "v":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        voltage = Vstd - (k*T/e)*data[pfx + 'phi_applied'][0]
        ffvec = data[pfx + 'ffrac_cathode'][0]
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
        sep = pfx + 'c_lyte_sep'
        trode = pfx + 'c_lyte_trode'
        datay_sep = data[sep][0]
        datay_trode = data[trode][0]
        Lsep = data[pfx + 'Lsep'][0] * 1e6
        Ltrode = data[pfx + 'Ltrode'][0] * 1e6
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
        sep = pfx + 'phi_lyte_sep'
        trode = pfx + 'phi_lyte_trode'
        datay_sep = data[sep][0]
        datay_trode = data[trode][0]
        Lsep = data[pfx + 'Lsep'][0] * 1e6
        Ltrode = data[pfx + 'Ltrode'][0] * 1e6
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
        psd_len = data[pfx + 'psd_lengths'][0]
        len_max = np.max(psd_len)
        len_min = np.min(psd_len)
        size_min = 0.10
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
#        title = ax.set_title("% = 00.0")
#        tmin = np.min(times)
#        tmax = np.max(times)
#        ax.patch.set_facecolor('gray')
        ax.patch.set_facecolor('white')
        ax.set_aspect('equal', 'box')
#        ax.set_aspect('equal')
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
#        yleft = fig.text(0.3, 0.5, "Separator",
#                rotation=90,
#                verticalalignment="center")
#        yright = fig.text(0.7, 0.5, "Current Collector",
#                rotation=90,
#                verticalalignment="center")
#        ax.invert_yaxis()
        rects = np.empty((Ntrode, numpart), dtype=object)
        cbar_mat = data[pfx + 'cbar_sld'][0]
        converter = matplotlib.colors.ColorConverter()
        rgba_green = converter.to_rgba('green')
        rgba_yellow = converter.to_rgba('yellow')
        rgba_red = converter.to_rgba('red')
        for (i,j),c in np.ndenumerate(cbar_mat):
            p_len = psd_len[i, j]
            if len_max == len_min:
                size_frac = 0.4
            else:
                size_frac = (p_len - len_min)/(len_max - len_min)
            size = ((size_frac) * (1-size_min) +
                    size_min)
            color = 'green'
            rects[i, j] = plt.Rectangle([i - size / 2, j - size / 2],
                    size, size, facecolor=color, edgecolor=color)
            ax.add_patch(rects[i, j])
        ax.autoscale_view()
        for tval in range(1,numtimes):
#            t_current = times[tval]
#            tfrac = (t_current - tmin)/(tmax - tmin) * 100
#            ax.set_title("% = {perc:2.1f}".format(perc=tfrac))
#            title.set_text("% = {perc:2.1f}".format(perc=tfrac))
#            ax.draw_artist(title)
            cbar_mat = data[pfx + 'cbar_sld'][tval]
            for (i,j),c in np.ndenumerate(cbar_mat):
                if c < 0.4:
                    color = 'green'
                    color_code = rgba_green
                elif c < 0.6:
                    color = 'yellow'
                    color_code = rgba_yellow
                else: # c > 0.6
                    color = 'red'
                    color_code = rgba_red
                fc = rects[i, j].get_facecolor()
                if fc != color_code:
                    rects[i, j].set_facecolor(color)
                    rects[i, j].set_edgecolor(color)
                    ax.draw_artist(rects[i, j])
            if save_flag:
                filename = os.path.join(outputdir,
                        'cbar{num:06d}.png'.format(num=tval))
                fig.savefig(filename)
            else:
                fig.canvas.blit(ax.bbox)
                time.sleep(1e-3)

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
        ymin = -1
        ymax = 10
        ax.set_xlabel('Position in cathode [um]')
        ax.set_ylabel('Potential of cathode [nondim]')
        cathp = pfx + 'phi_cath'
        Ltrode = data[pfx + 'Ltrode'][0] * 1e6
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
        raise Exception("Unexpected plot type argument. "
                + "Try 'v', 'elytec', 'elytep', 'cbar', 'csld', "
                + "'cathp'.")
    return

if __name__ == "__main__":
    # Get input file from script parameters
    if len(sys.argv) < 2:
        raise Exception("Need input data file name")
    infile = sys.argv[1]
    if not os.path.isfile(os.path.join(os.getcwd(), infile)):
        raise Exception("Input file doesn't exist")
    # Get plot type from script parameters
    if len(sys.argv) > 2:
        plot_type = sys.argv[2]
    else:
        plot_type = "v"
    # Save the plot instead of showing on screen?
    # Get from script parameters
    save_flag = False
    if len(sys.argv) > 3:
        if sys.argv[3] == "save":
            save_flag = True
    plot(infile, plot_type, save_flag)
