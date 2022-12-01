Analyze the results with built-in tools
=======================================


Analyze output with ``mpetplot.py``. Pass the output data directory, then use the optional plotting arguments. The options for ``mpetplot.py`` are:
   #. ``-pt`` for plotting types
   #. ``-t`` for saving output to text format
   #. ``-s`` for options to save the plot
   #. ``-c`` for color_map options that are used with plot type ``cbar_{full,c,a}``
   #. ``-st`` to specify the smooth colormap used with plot type ``cbar_{full,c,a}``

1.  Analyze output with plots using ``mpetplot.py``. Pass output data directory, then use ``-pt [plottype]`` with one (or more) of the plot types listed below. Default is ``v``.
    - e.g., voltage plot: ``mpetplot.py sim_output -pt v``
    - other options (``full``, ``c``, ``a`` indicate full cell, cathode, and anode):
   #. ``v`` or ``vt`` -- voltage vs filling fraction or vs time
   #. ``curr`` -- current vs time
   #. ``elytec{f}`` -- electrolyte concentration (movie) or final snapshot with f
   #. ``elytep{f}`` -- electrolyte potential (movie) or final snapshot with f
   #. ``elytei{f}`` -- electrolyte current density (movie) or final snapshot with f
   #. ``surf_{c,a}`` -- solid surface concentrations
   #. ``soc_{c,a}`` -- overall utilization / state of charge of electrode
   #. ``csld_{c,a}`` -- solid concentrations of particles in electrode (movie; used with solidType_{c,a} not homog)
   #. ``cbarLine_{c,a}`` -- average concentration in each particle of electrode
   #. ``bulkp_{c,a}`` -- macroscopic electrode solid phase potential (movie)
   #. ``cbar_{full,c,a}`` -- average solid concentrations as changing colors (movie)
    - There are two options for the color map type that is used: ``smooth`` or ``discrete``. This can be set with the ``-c`` option, e.g., ``mpetplot.py sim_output -pt cbar_full -c discrete``. The default value is ``discrete``.
    - When using the ``smooth`` color map option, the colors are selected from colormao_custom.npz, which includes three options (``GnYlRd_1``, ``GnYlRd_2``, and ``GnYlRd_3``) that can be selected with the ``st`` option, e.g., ``mpetplot.py sim_output -pt cbar_full -c discrete -st GnYlRd_1``. The default value is ``GnYlRd_3``.

2.  Alternatively, convert the output to plain text (csv) format using : ``mpetplot.py sim_output text`` (or replace sim_output with any subfolder in the history folder).
Then analyze using whatever tools you prefer. If you want to save output to a movie (or figure), add save as an extra argument to ``mpetplot.py``: ``mpetplot.py sim_output cbar save``.

Movie output requires that you have ``ffmpeg`` or ``mencoder`` (part of ``MPlayer``) installed.
