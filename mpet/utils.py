def set_plot_defaults(mpl):
    axtickfsize = 18
    labelfsize = 20
    legfsize = labelfsize - 2
    txtfsize = labelfsize - 2
    lwidth = 3
    markersize = 10
    markeredgewidth = 0.1
    mpl.rcParams['xtick.labelsize'] = axtickfsize
    mpl.rcParams['ytick.labelsize'] = axtickfsize
    mpl.rcParams['axes.labelsize'] = labelfsize
    mpl.rcParams['font.size'] = txtfsize
    mpl.rcParams['legend.fontsize'] = legfsize
    mpl.rcParams['lines.linewidth'] = lwidth
    mpl.rcParams['lines.markersize'] = markersize
    mpl.rcParams['lines.markeredgewidth'] = markeredgewidth
#    mpl.rcParams['text.usetex'] = True
