Analyze the results with built-in tools
=============================================


Analyze output with mpetplot.py (pass output data directory, then plot-type as arguments)

 * e.g., voltage plot: ``mpetplot.py sim_output v``
 * other options (full, c, a indicate full cell, cathode, and anode):

   #. ``v`` or ``vt`` -- voltage vs filling fraction or vs time
   #. ``curr`` -- current vs time
   #. ``elytec{f}`` -- electrolyte concentration (movie) or final snapshot with f
   #. ``elytep{f}`` -- electrolyte potential (movie) or final snapshot with f
   #. ``elytei{f}`` -- electrolyte current density (movie) or final snapshot with f
   #. ``surf_{c,a}`` -- solid surface concentrations
   #. ``soc_{c,a}`` -- overall utilization / state of charge of electrode
   #. ``csld_{c,a}`` -- solid concentrations of particles in electrode (movie; used with solidType_{c,a} not homog)
   #. ``cbarLine_{c,a}`` -- average concentration in each particle of electrode
   #. ``cbar_{full,c,a}`` -- average solid concentrations as changing colors (movie)
   #. ``bulkp_{c,a}`` -- macroscopic electrode solid phase potential (movie)




Alternatively, convert the output to plain text (csv) format using : ``mpetplot.py sim_output text`` (or replace sim_output with any subfolder in the history folder).
Then analyze using whatever tools you prefer. If you want to save output to a movie (or figure), add save as an extra argument to ``mpetplot.py``: ``mpetplot.py sim_output cbar save``.

Movie output requires that you have ``ffmpeg`` or ``mencoder`` (part of ``MPlayer``) installed.