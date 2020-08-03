params_a='''[Particles]
type = diffn
discretization = 50.0e-8
shape = sphere
thickness = 20.0e-9

[Material]
muRfunc = LiC6_coke_ss2
logPad = false
noise = false
noise_prefac = 1e-6
numnoise = 200
Omega_a = 1.8560e-20
Omega_b = 5.761532e-21
Omega_c = 8.23076e-20
kappa = 5.0148e-10
B = 0.1916e9
EvdW = 0.0
rho_s = 1.58981e28
D = 5e-13
Dfunc = constant
dgammadc = 0e-30
cwet = 0.98

[Reactions]
rxnType = BV_mod02
k0 = 8.3038
alpha = 0.5
lambda = 6.26e-20
Rfilm = 0e-0
'''
params_c='''[Particles]
type = diffn
discretization = 5.0e-8
shape = sphere
thickness = 20.0e-9

[Material]
muRfunc = LiMn2O4_ss2
logPad = false
noise = false
noise_prefac = 1e-6
numnoise = 200
Omega_a = 1.8560e-20
Omega_b = 5.761532e-21
Omega_c = 8.23076e-20
kappa = 5.0148e-10
B = 0.1916e9
EvdW = 0.0
rho_s = 1.42842e28
D = 1.0e-13
Dfunc = constant
dgammadc = 0e-30
cwet = 0.98

[Reactions]
rxnType = BV_mod01
k0 = 7.225
alpha = 0.5
lambda = 6.26e-20
Rfilm = 0e-0
'''
params_system='''[Sim Params]
# Constant voltage or current or segments of one of them
profileType = CC
# Battery (dis)charge c-rate (only used for CC), number of capacities / hr
Crate = 0.72
# Voltage cutoffs, V
Vmax = 5
Vmin = 2
# Fraction of full (dis)charge to simulate (only used for CC)
capFrac = 0.8
Vset = 0.12
segments = [(0.3, 0.4),(-0.5, 0.1)]
tramp = 1e-3
prevDir = false
tend = 1.2e3
tsteps = 200
relTol = 1.0e-6
absTol = 1.0e-6
T = 298
randomSeed = true
seed = 0
Rser = 0
# Cathode, anode, and separator numer disc. in x direction (volumes in electrodes)
Nvol_c = 2
Nvol_s = 2
Nvol_a = 2
# Number of particles per volume for cathode and anode
Npart_c = 1
Npart_a = 1

[Electrodes]
# The name of the parameter file describing the cathode particles
cathode = params_c.cfg
# The name of the parameter file describing the anode particles
anode = params_a.cfg
k0_foil = 1e0
Rfilm_foil = 0e-0

[Particles]
mean_c = 1e-6
stddev_c = 0
mean_a = 18e-6
stddev_a = 0
# Initial electrode filling fractions
# (for disch, anode starts full, cathode starts empty)
cs0_c = 0.2
cs0_a = 0.495

[Conductivity]
simBulkCond_c = false
simBulkCond_a = false
sigma_s_c = 1e-1
sigma_s_a = 1e-1
simPartCond_c = false
simPartCond_a = false
G_mean_c = 1e-14
G_stddev_c = 0
G_mean_a = 1e-14
G_stddev_a = 0

[Geometry]
# Thicknesses, m
L_c = 200e-6
L_a = 243e-6
L_s = 50e-6
# Volume loading percents of active material (volume fraction of solid
# that is active material)
P_L_c = 0.849
P_L_a = 0.956
# Porosities (liquid volume fraction in each region)
poros_c = 0.3
poros_a = 0.3
poros_s = 0.4
# Bruggeman exponent (tortuosity = porosity^bruggExp)
BruggExp_c = -0.5
BruggExp_a = -0.5
BruggExp_s = -0.5

[Electrolyte]
# Initial electrolyte conc., mol/m^3
c0 = 1000
# Cation/anion charge number (e.g. 2, -1 for CaCl_2)
zp = 1
zm = -1
# Cation/anion dissociation number (e.g. 1, 2 for CaCl_2)
nup = 1
num = 1
# Electrolyte model
elyteModelType = SM
SMset = LiClO4_PC
n = 1
sp = -1
Dp = 2.2e-10
Dm = 2.94e-10
'''

fs=open('params_system.cfg','w+')
fs.write(params_system)
fs.close()

fa=open('params_a.cfg','w+')
fa.write(params_a)
fa.close()

fc=open('params_c.cfg','w+')
fc.write(params_c)
fc.close()
	