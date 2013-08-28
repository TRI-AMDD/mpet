ACR w/ surface wetting and coherency strain
[t,cpcs,ffvec,vvec,disc,part] = acr_mpet_swcs_rev2(dim_crate,part,cpcs0,disc,ffend)

OUTPUTS
t = time
cpcs = results vector
ffvec = filling fraction
vvec = voltage
disc = discretizations
part = particle sizes

INPUTS
dim_crate = dimensional C-rate
part = particle sizes (enter 0 if you dont have one)
cpcs0 = start point (enter 0 if you dont have one, or use last time from previous output)
disc = discretizations (required for use with part vector, enter 0 if not applicable)
ffend = ending filling fraction of discretization (1 for full discharge, 0 for full charge)

Solid solution model with variable particle sizes
[t,cpcs,disc,psd,ffvec,vvec] = pm_2dvarsz_rev1(disc,psd)

OUTPUTS
t = time
cpcs = results vector
disc = discretizations
psd = particle size vector
ffvec = filling fraction
vvec = voltage

INPUTS
disc = discretizations vector (0 if not applicable)
psd = particle sizes (0 if not applicable)

The animate scripts are straightforward.  Just make sure you check the 'fig' and 'output' variables.  Fig denotes the plot variables and output is a boolean (0 or 1) denoting whether or not to record the video.
