%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
k = 1.381e-23;      % Boltzmann constant
T = 298;            % Temp, K
e = 1.602e-19;      % Charge of proton, C
Na = 6.02e23;       % Avogadro's number
F = e*Na;           % Faraday's number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DIMENSIONAL VALUES HERE

% Discharge settings
dim_crate = 0.1;                  % C-rate (electrode capacity per hour)
dim_io = 0.1;                      % Exchange current density, A/m^2 (0.1 for H2/Pt)

% Electrode properties
Lx = 50e-6;                         % electrode thickness, m
Lsep = 25e-6;                        % separator thickness, m
Asep = 1e-4;                        % area of separator, m^2
Lp = 0.69;                          % Volume loading percent active material
poros = 0.4;                        % Porosity
c0 = 1000;                          % Initial electrolyte conc., mol/m^3
zp = 1;                             % Cation charge number
zm = 1;                             % Anion charge number
Dp = 2.2e-10;                       % Cation diff, m^2/s, LiPF6 in EC/DMC
Dm = 2.94e-10;                      % Anion diff, m^2/s, LiPF6 in EC/DMC

% Particle size distribution
mean = 160e-9;                      % Average particle size, m
stddev = 20e-9;                     % Standard deviation, m

% Material properties
dim_a = 1.8560e-20;                 % Regular solution parameter, J
dim_kappa = 5.0148e-10;             % Gradient penalty, J/m
dim_b = 0.1916e9;                   % Stress, Pa
% dim_b = 0; dim_kappa = 0;

rhos = 1.3793e28;                   % site density, 1/m^3
cwet = 0.98;                        % Dimensionless wetted conc.
wet_thick = 2e-9;                   % Thickness of wetting on surf.
Vstd = 3.422;                       % Standard potential, V
alpha = 0.5;                        % Charge transfer coefficient

% Discretization settings
Nx = 5;                            % Number disc. in x direction
solid_disc = 1e-9;                 % Discretization size of solid, m (MUST BE LESS THAN LAMBDA)
numpart = 4;                       % Particles per volume
tsteps = 200;                      % Number disc. in time
ffend = 0.95;                       % Final filling fraction

save('./Data/params.mat')
