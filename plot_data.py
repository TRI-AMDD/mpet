import sys
import os
import time

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import matplotlib.collections as mcollect

def show_data(infile, plot_type, save_flag):
    pfx = 'mpet.'
    ttl_fmt = "% = {perc:2.1f}"
    data = sio.loadmat(infile)
    # Pick out some useful parameters
    Vstd = data[pfx + 'Vstd_c'][0][0]  # Standard potential, V
    k = data[pfx + 'k'][0][0]        # Boltzmann constant
    Tref = data[pfx + 'Tref'][0][0]     # Temp, K
    e = data[pfx + 'e'][0][0]        # Charge of proton, C
    td = data[pfx + 'td'][0][0]     # diffusive time
    # Discretization info
    Ntrode = int(data[pfx + 'NumTrode'][0][0])
    Nsep = int(data[pfx + 'NumSep'][0][0])
    numpart = int(data[pfx + 'NumPart'][0][0])
    # Extract the reported simulation times
    times = data[pfx + 'phi_applied_times'][0]
    numtimes = len(times)
    tmin = np.min(times)
    tmax = np.max(times)
    # Simulation type
    sim_ACR = data[pfx + "type_ACR"][0][0]
    sim_homog = data[pfx + "type_homog"][0][0]
    shape_sphere = data[pfx + "shape_sphere"][0][0]
    shape_C3 = data[pfx + "shape_C3"][0][0]

    # Print relevant simulation info
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
#    print "lambda:", data[pfx + "lambda_c"][0][0]
    print "dim_k0 [A/m^2]:", data[pfx + "dim_k0"][0][0]
    mcond_bool = data[pfx + "cath_bulk_cond"][0][0]
    if mcond_bool:
        print ("cathode bulk conductivity loss: Yes -- " +
                "dim_mcond [S/m]: " +
                str(data[pfx + "dim_mcond"][0][0]))
    else:
        print "cathode bulk conductivity loss: No"
    scond_bool = data[pfx + "cath_surf_cond"][0][0]
    if scond_bool:
        print ("cathode surface conductivity loss: Yes -- " +
                "dim_scond [S]: " +
                str(data[pfx + "dim_scond"][0][0]))
    else:
        print "cathode surface conductivity loss: No"

    # Plot voltage profile
    if plot_type == "v":
        fig, ax = plt.subplots()
        voltage = Vstd - (k*Tref/e)*data[pfx + 'phi_applied'][0]
        ffvec = data[pfx + 'ffrac_cathode'][0]
        ax.plot(ffvec, voltage)
        xmin = np.min(ffvec)
        xmax = np.max(ffvec)
        ax.axhline(y=Vstd, xmin=xmin, xmax=xmax, linestyle='--', color='g')
        ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        ax.set_ylabel("Voltage [V]")
        ax.set_ylim((Vstd - 0.3, Vstd + 0.4))
        if save_flag:
            fig.savefig("mpet_v.png")
        plt.show()
        return

    # Plot current profile
    if plot_type == "curr":
        fig, ax = plt.subplots()
        current = data[pfx + "current"][0] * 3600/td
        ffvec = data[pfx + 'ffrac_cathode'][0]
#        ax.plot(ffvec, current)
        ax.plot(times, current)
        xmin = np.min(ffvec)
        xmax = np.max(ffvec)
        ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        ax.set_ylabel("Current [C-rate]")
#        ax.set_ylim((Vstd - 0.3, Vstd + 0.4))
        if save_flag:
            fig.savefig("mpet_current.png")
        plt.show()
        return

    # Plot electrolyte concentration
    elif plot_type == "elytec":
        matplotlib.animation.Animation._blit_draw = _blit_draw
        fig, ax = plt.subplots()
        ymin = 0
        ymax = 1.5
        ax.set_xlabel('Battery Position [um]')
        ax.set_ylabel('Concentration of electrolyte [nondim]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
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
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            datay_sep = data[sep][tind]
            datay_trode = data[trode][tind]
            datay = np.hstack((datay_sep, datay_trode))
            line1.set_ydata(datay)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return line1, ttl

    # Plot electrolyte electrostatic potential
    elif plot_type == "elytep":
        matplotlib.animation.Animation._blit_draw = _blit_draw
        fig, ax = plt.subplots()
        ymin = -5
        ymax = 5
        ax.set_xlabel('Battery Position [um]')
        ax.set_ylabel('Potential of electrolyte [nondim]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
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
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            datay_sep = data[sep][tind]
            datay_trode = data[trode][tind]
            datay = np.hstack((datay_sep, datay_trode))
            line1.set_ydata(datay)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return line1, ttl

    # Plot all solid concentrations or potentials
    elif (plot_type == "csld") or (plot_type == "phisld"):
        fig, ax = plt.subplots(numpart, Ntrode, squeeze=False,
                sharey=True)
        sol = np.empty((numpart, Ntrode), dtype=object)
        lens = np.zeros((numpart, Ntrode))
        lines = np.empty((numpart, Ntrode), dtype=object)
        if plot_type == "csld":
            str_base = "mpet.solid_c_vol{j}_part{i}"
            ylim = (0, 1.01)
        else: # plot_type == "phisld"
            str_base = "mpet.solid_p_vol{j}_part{i}"
            ylim = (-10, 20)
        for i in range(numpart):
            for j in range(Ntrode):
                sol[i, j] = str_base.format(i=i, j=j)
                lens[i, j] = data["mpet.psd_lengths"][0][j, i]
                # Remove axis ticks
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
#                ax[i, j].yaxis.set_major_locator(plt.NullLocator())
                datay = data[sol[i, j]][0]
                numy = len(datay)
                datax = np.linspace(0, lens[i, j], numy)
                ax[i, j].set_ylim(ylim)
                ax[i, j].set_xlim((0, lens[i, j]))
                line, = ax[i, j].plot(datax, datay)
                lines[i, j] = line
        def init():
            for i in range(numpart):
                for j in range(Ntrode):
                    datax = np.zeros(data[sol[i, j]][0].shape)
                    lines[i, j].set_ydata(np.ma.array(datax, mask=True))
            return tuple(lines.reshape(-1))
        def animate(tind):
            for i in range(numpart):
                for j in range(Ntrode):
                    datay = data[sol[i, j]][tind]
                    lines[i, j].set_ydata(datay)
            return tuple(lines.reshape(-1))

    # Plot average solid concentrations
    elif plot_type == "cbar":
        # Set up colors.
        # Define if you want smooth or discrete color changes
        # Option: "smooth" or "discrete"
        color_changes = "discrete"
        # Discrete color changes:
        if color_changes == "discrete":
            to_yellow = 0.4
            to_red = 0.6
            # Make a discrete colormap that goes from green to yellow
            # to red instantaneously
            cdict = {
                    "red" : [(0.0, 0.0, 0.0),
                             (to_yellow, 0.0, 1.0),
                             (1.0, 1.0, 1.0)],
                    "green" : [(0.0, 0.502, 0.502),
                               (to_yellow, 0.502, 1.0),
                               (to_red, 1.0, 0.0),
                               (1.0, 0.0, 0.0)],
                    "blue" : [(0.0, 0.0, 0.0),
                              (1.0, 0.0, 0.0)]
                    }
            cmap = matplotlib.colors.LinearSegmentedColormap(
                    "discrete", cdict)
        # Smooth colormap changes:
        if color_changes == "smooth":
#            cmap = matplotlib.cm.RdYlGn_r # A default green-yellow-red map
            # generated with colormap.org
            cmaps = np.load("colormaps_custom.npz")
            cmap_data = cmaps["GnYlRd_3"]
            cmap = matplotlib.colors.ListedColormap(cmap_data/255.)

        # Implement hack to be able to animate title
        matplotlib.animation.Animation._blit_draw = _blit_draw
        # Get particle sizes (and max size) (length-based)
        psd_len = data[pfx + 'psd_lengths'][0]
        len_max = np.max(psd_len)
        len_min = np.min(psd_len)
        size_frac_min = 0.10
        fig, ax = plt.subplots()
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        ax.patch.set_facecolor('white')
        # Don't stretch axes to fit figure -- keep 1:1 x:y ratio.
        ax.set_aspect('equal', 'box')
        # Don't show axis ticks
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_xlim(0, 1.)
        ax.set_ylim(0, float(numpart)/Ntrode)
        # Label parts of the figure
        ylft = ax.text(-0.07, 0.5, "Separator",
                transform=ax.transAxes, rotation=90,
                verticalalignment="center",
                horizontalalignment="center")
        yrht = ax.text(1.09, 0.5, "Current Collector",
                transform=ax.transAxes, rotation=90,
                verticalalignment="center",
                horizontalalignment="center")
        xbtm = ax.text(.50, -0.05, "Electrode Depth -->",
                transform=ax.transAxes, rotation=0,
                verticalalignment="center",
                horizontalalignment="center")
        # Geometric parameters for placing the rectangles on the axes
        spacing = 1.0 / Ntrode
        size_fracs = 0.4*np.ones((Ntrode, numpart))
        if len_max != len_min:
            size_fracs = (psd_len - len_min)/(len_max - len_min)
        sizes = (size_fracs * (1 - size_frac_min) + size_frac_min) / Ntrode
        # Create rectangle "patches" to add to figure axes.
        rects = np.empty((Ntrode, numpart), dtype=object)
        color = 'green' # value is irrelevant -- it will be animated
        for (i, j), c in np.ndenumerate(sizes):
            size = sizes[i, j]
            center = np.array([spacing*(i + 0.5), spacing*(j + 0.5)])
            bottom_left = center - size / 2
            rects[i, j] = plt.Rectangle(bottom_left,
                    size, size, color=color)
        # Create a group of rectange "patches" from the rects array
        collection = mcollect.PatchCollection(rects.reshape(-1),
                animated=True)
        # Put them on the axes
        ax.add_collection(collection)
        # Have a "background" image of rectanges representing the
        # initial state of the system.
        def init():
            cbar_mat = data[pfx + 'cbar_sld'][0]
            colors = cmap(cbar_mat.reshape(-1))
            collection.set_color(colors)
            ttl.set_text('')
            return collection, ttl
        def animate(tind):
            cbar_mat = data[pfx + 'cbar_sld'][tind]
            colors = cmap(cbar_mat.reshape(-1))
            collection.set_color(colors)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return collection, ttl

    # Plot cathode potential
    elif plot_type == "cathp":
        matplotlib.animation.Animation._blit_draw = _blit_draw
        fig, ax = plt.subplots()
        ymin = -1
        ymax = 10
        ax.set_xlabel('Position in cathode [um]')
        ax.set_ylabel('Potential of cathode [nondim]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
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
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            datay = data[cathp][tind]
            line1.set_ydata(datay)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return line1, ttl

    else:
        raise Exception("Unexpected plot type argument. " +
                "Try 'v', 'elytec', 'elytep', 'cbar', 'csld', " +
                "'cathp'.")

    ani = manim.FuncAnimation(fig, animate, frames=numtimes,
            interval=50, blit=True, repeat=False, init_func=init)
    if save_flag:
        ani.save("mpet_{type}.mp4".format(type=plot_type), fps=30)
#                extra_args=['-vcodec', 'libx264'])
    plt.show()

    return

# This is a block of code which messes with some matplotlib internals
# to allow for animation of a title. See
# http://stackoverflow.com/questions/17558096/animated-title-in-matplotlib
def _blit_draw(self, artists, bg_cache):
    # Handles blitted drawing, which renders only the artists given instead
    # of the entire figure.
    updated_ax = []
    for a in artists:
        # If we haven't cached the background for this axes object, do
        # so now. This might not always be reliable, but it's an attempt
        # to automate the process.
        if a.axes not in bg_cache:
            # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            # change here
            bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
        a.axes.draw_artist(a)
        updated_ax.append(a.axes)

    # After rendering all the needed artists, blit each axes individually.
    for ax in set(updated_ax):
        # and here
        # ax.figure.canvas.blit(ax.bbox)
        ax.figure.canvas.blit(ax.figure.bbox)

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
    show_data(infile, plot_type, save_flag)
