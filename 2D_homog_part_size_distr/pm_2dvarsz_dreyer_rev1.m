function [t,cs,disc,psd,ffvec,vvec] = pm_2dvarsz_dreyer_rev1(disc,psd)

% This script simulates a 2D electrode with variable size particles.  The
% particles are all homogeneous and use the regular solution model (ONLY).
% The purpose of this script is to simulate the case of a simple constant
% non-monotonic OCP for all particles, with particle size effects.

% The user enters a C-rate and particle size distribution.  The area:volume
% ratio can be user defined but is assumed to be either plate particles or
% spherical particles.

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
dim_crate = .01;                    % C-rate (electrode capacity per hour)
dim_io = 1;                         % Exchange current density, A/m^2 (0.1 for H2/Pt)

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
Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm);   % Ambipolar diffusivity
tp = zp*Dp / (zp*Dp + zm*Dm);       % Cation transference number

% Particle size distribution
mean = 160e-9;                      % Average particle size, m
stddev = 20e-9;                     % Standard deviation, m

% Material properties
dim_a = 1.8560e-20;                 % Regular solution parameter, J
% dim_kappa = 5.0148e-10;             % Gradient penalty, J/m
% dim_b = 0.1916e9;                   % Stress, Pa
% dim_b = 0; dim_kappa = 0;

rhos = 1.3793e28;                   % site density, 1/m^3
csmax = rhos/Na;                    % maximum concentration, mol/m^3
% cwet = 0.98;                        % Dimensionless wetted conc.
% wet_thick = 2e-9;                   % Thickness of wetting on surf.
Vstd = 3.422;                       % Standard potential, V
alpha = 0.5;                        % Charge transfer coefficient

% Discretization settings
Nx = 20;                            % Number disc. in x direction
Ny = 10;                            % Number disc. in y direction
tsteps = 200;                       % Number disc. in time
ffend = .4;                          % Final filling fraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we take care of our particle size distributions
num_particles = Nx * Ny;
if (max(size(psd))>1)
    part_size = psd;
    ap = 3.475 ./ part_size;
else
    part_size = abs(normrnd(mean,stddev,num_particles,1));
    ap = 3.475 ./ part_size;     % Area:volume ratio for plate particles
end
psd = part_size;

% Now we calculate the dimensionless quantities used in the simulation
td = Lx^2 / Damb;       % Diffusive time
currset = dim_crate * (td/3600);

if currset ~= 0
    tr = linspace(0,1/abs(currset),tsteps);
else    
    tr = linspace(0,30,100);
end
io = (ap .* dim_io .* td) ./ (F .* csmax);

% Need some noise
% noise = 0.0000001*randn(max(size(tr)),Nx*Ny);
noise = zeros(max(size(tr)),Nx*Ny);

% Material properties
% kappa = dim_kappa ./ (k*T*rhos*part_size.^2);
a = dim_a / (k*T);
% b = dim_b / (k*T*rhos);

% Set the discretization values
if ~isa(disc,'struct')
    sf = Lsep/Lx;
    ssx = ceil(sf*Nx);
    ss = ssx*Ny;
    disc = struct('ss',ss,'steps',Nx*Ny,...
                    'len',2*(ss+Nx*Ny)+Nx*Ny+1, ...
                    'sol',2*(ss+Nx*Ny)+1,'Nx',Nx,'Ny',Ny,'sf',sf);
else
    Nx = disc.Nx;
    Ny = disc.Ny;
end

cs0 = 0.01;                 
phi_init = calcmu(cs0,a);
% Assemble it all
csinit = zeros(Nx*Ny+1,1);
csinit(1:end-1) = cs0;
csinit(end) = phi_init;

% Before we can call the solver, we need a Mass matrix
M = genMass(Nx,Ny);

% Prepare to call the solver
% options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none');
options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none','Events',@events);
disp('Calling ode15s solver...')
[t,cs]=ode15s(@calcRHS,tr,csinit,options,io,currset,a,alpha,disc,tr,ffend,noise);

% Now we analyze the results before returning
disp('Done.')                
disp('Calculating the voltage and filling fraction vectors...')                
                
% First we calculate the voltage                 
vvec = Vstd - (k*T/e)*cs(:,end);

% Now the filling fraction vector - we only care about the active parts of
% the particle.  That is, we ignore the surface wetting as it does not move
ffvec = zeros(max(size(t)),1);
for i=1:max(size(t))
    ffvec(i) = sum(cs(i,1:end-1))/(Nx*Ny);
end
disp('Finished.')

return;

function val = calcRHS(t,cs,io,currset,a,alpha,disc,tr,ffend,noise)

% Initialize output
val = zeros(max(size(cs)),1);             
             
% REACTION RATE OF PARTICLES
muvec = calcmu(cs(1:end-1),a);
ecd = io.*sqrt(exp(muvec)).*(1-cs(1:end-1));
eta = muvec-cs(end);
val(1:end-1) = ecd.*(exp(-alpha.*eta)-exp((1-alpha).*eta));
val(1:end-1) = val(1:end-1) + interp1q(tr,noise,t)';

% CURRENT CONDITION
val(end) = currset;
val = real(val);

return;

function mu = calcmu(cs,a)

% This function calculates the chemical potential
mu = log(cs./(1-cs))+a.*(1-2.*cs);

return;

function M = genMass(Nx,Ny)
    
% Initialize
M = sparse(Nx*Ny+1,Nx*Ny+1);
M(1:end-1,1:end-1) = speye(Nx*Ny);
M(end,1:end-1) = 1./(Nx*Ny);

return;

function [value, isterminal, direction] = events(t,cs,io,currset,a,alpha,disc,tr,ffend,noise)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)      

% Calculate the filling fraction 
ffvec = sum(cs(1:end-1))/(disc.Nx*disc.Ny);
value = ffvec - ffend;
isterminal = 1;
direction = 0;
                        
return;