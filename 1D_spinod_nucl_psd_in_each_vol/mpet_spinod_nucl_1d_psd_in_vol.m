function [t,cpcs,csmat,disc,psd,ffvec,vvec] = mpet_spinod_nucl_1d_psd_in_vol(psd,disc)

% This script simulates a 1D electrode with variable size particles.  The
% particles are all homogeneous and use the regular solution model (ONLY).
% The purpose of this script is to simulate the case of a simple constant
% non-monotonic OCP for all particles, with particle size effects.

% Each volume contains a number of particles taken from a distribution of
% particle sizes.  A function for area and volume are required.

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
dim_crate = 1;                    % C-rate (electrode capacity per hour)
dim_io = .1;                         % Exchange current density, A/m^2 (0.1 for H2/Pt)

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
Nx = 25;                            % Number disc. in x direction
numpart = 50;                       % Particles per volume
tsteps = 200;                       % Number disc. in time
ffend = .95;                          % Final filling fraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we take care of our particle size distributions
if max(size(psd))==1
    totalpart = Nx * numpart;
    psd = abs(normrnd(mean,stddev,totalpart,1));
    pareavec = (4*pi).*psd.^2;
    pvolvec = (4/3).*pi.*psd.^3;
else
    totalpart = disc.steps*disc.numpart;
    pareavec = (4*pi).*psd.^2;
    pvolvec = (4/3).*pi.*psd.^3;
end

% Now we calculate the dimensionless quantities used in the simulation
td = Lx^2 / Damb;       % Diffusive time
nDp = Dp / Damb;
nDm = Dm / Damb;
currset = dim_crate * (td/3600);

if currset ~= 0
    tr = linspace(0,1/abs(currset),tsteps);
else    
    tr = linspace(0,30,100);
end
io = ((pareavec./pvolvec) .* dim_io .* td) ./ (F .* csmax);
epsbeta = ((1-poros)*Lp*csmax) / c0; % Vs/V*csmax/c0 = poros*beta

% Need some noise
% noise = 0.0000001*randn(max(size(tr)),Nx*Ny);
noise = zeros(max(size(tr)),totalpart);
a = dim_a / (k*T);

% Set the discretization values
if ~isa(disc,'struct')
    sf = Lsep/Lx;
    ssx = ceil(sf*Nx);
    ss = ssx;
    disc = struct('ss',ss,'steps',Nx,'numpart',numpart,...
                    'len',2*(ss+Nx)+totalpart+1, ...
                    'sol',2*(ss+Nx)+1,'Nx',Nx,'sf',sf);                      
else
    Nx = disc.Nx;
end

% Calculate the composition bounds within which the particles are
% phase separated.
% cslow is the spinodal point at the lower composition
% (other is 0.5*(1+sqrt(1-2/a)) )
cslow = 0.5*(1-sqrt(1-2/a));
% csup is the binodal, the upper comp where mu_homogeneous = 0
csup = fsolve(@(cs) log(cs/(1-cs))+a*(1-2*cs), 0.995);
% disp('cslow')
% disp(cslow)
% disp('csup')
% disp(csup)
% if cslow > csup
%     error('cs bound calcs wrong')
% end
% cstestvec = 0.001:0.001:0.999;
% plot(cstestvec, calcmu(cstestvec,a,cslow,csup))
% return
% Create global information for nucleation events
global nucl_vec t_last
t_last = 0;
nucl_vec = ones(totalpart,1); % Everything starts homogeneous
nu_nucl = 1e-10*td; % dimensionless (dimensional = 1e-10 1/sec)
mu_crit = log(cslow/(1-cslow))+a*(1-2*cslow);

cs0 = 0.01;                 
%phi_init = calcmu(cs0,a,cslow,csup,nu_nucl,mu_crit,t_last);
phi_init = log(cs0/(1-cs0))+a.*(1-2.*cs0);
cinit = 1;
% Assemble it all
cpcsinit = zeros(disc.len,1);
cpcsinit(1:disc.ss+disc.steps) = cinit;
cpcsinit(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = phi_init;
cpcsinit(disc.sol:end-1) = cs0;
cpcsinit(end) = phi_init;

% Before we can call the solver, we need a Mass matrix
M = genMass(disc,poros,Nx,epsbeta,tp,pvolvec);

% Porosity vector
porosvec = ones(disc.ss+disc.steps+1,1);
porosvec(disc.ss+1:end) = poros;
porosvec = porosvec.^(3/2);     % Bruggeman

% Prepare to call the solver
% options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none');
options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none',...
    'RelTol',5e-2,'AbsTol',1e-3','Events',@events);
disp('Calling ode15s solver...')
[t,cpcs]=ode15s(@calcRHS,tr,cpcsinit,options,io,currset,a,alpha,porosvec,numpart,...
                 Nx,disc,tp,zp,zm,nDp,nDm,tr,epsbeta,cslow,csup,nu_nucl,mu_crit,...
                 ffend,noise);

% Now we analyze the results before returning
disp('Done.')                
disp('Calculating the voltage and filling fraction vectors...')                
                
% First we calculate the voltage                 
vvec = Vstd - (k*T/e)*cpcs(:,end);

% Now the filling fraction vector - we only care about the active parts of
% the particle.  That is, we ignore the surface wetting as it does not move
ffvec = zeros(max(size(t)),1);
for i=1:max(size(t))
    ffvec(i) = sum(pvolvec.*cpcs(i,disc.sol:end-1)')/sum(pvolvec);
end

% Create cs matrix
csmat = zeros(max(size(t)),Nx,numpart);
for i=1:max(size(t))
    for j=1:Nx
        for k=0:numpart-1
            csmat(i,j,k+1) = cpcs(i,disc.sol+(j-1)*numpart+k);
        end
    end
end

disp('Finished.')

return;

function val = calcRHS(t,cpcs,io,currset,a,alpha,porosvec,numpart,...
                 Nx,disc,tp,zp,zm,nDp,nDm,tr,epsbeta,cslow,csup,nu_nucl,mu_crit,...
                 ffend,noise)

% Initialize output
val = zeros(max(size(cpcs)),1);             
             
% Pull out the concentrations first
cvec = cpcs(1:disc.ss+disc.steps);
phivec = cpcs(disc.ss+disc.steps+1:2*(disc.ss+disc.steps));
phi0 = cpcs(end);
csvec = cpcs(disc.sol:end-1);

% MASS CONSERVATION - ELECTROLYTE DIFFUSION
ctmp = zeros(disc.ss+disc.steps+2,1);
ctmp(2:end-1) = cvec;
% Boundary conditions
ctmp(1) = ctmp(2) + currset*epsbeta*(1-tp)/Nx;
ctmp(end) = ctmp(end-1);
% Porosity effects
cflux = -porosvec.*diff(ctmp).*Nx;
val(1:disc.ss+disc.steps) = -diff(cflux).*Nx;

% CHARGE CONSERVATION - DIVERGENCE OF CURRENT DENSITY
phitmp = zeros(disc.ss+disc.steps+2,1);
phitmp(2:end-1) = phivec;
% Boundary conditions
phitmp(1) = phi0;
phitmp(end) = phitmp(end-1);
% Current density
cavg = (ctmp(1:end-1)+ctmp(2:end))/2;
currdens = -((zp*nDp-zm*nDm).*diff(ctmp).*Nx) - ...
                ((zp*nDp+zm*nDm).*cavg.*diff(phitmp).*Nx);
val(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = -diff(porosvec.*currdens).*Nx;                                                            

% REACTION RATE OF PARTICLES
rxncmat = repmat(cvec(disc.ss+1:end),1,numpart);
rxnphimat = repmat(phivec(disc.ss+1:end),1,numpart);
muvec = calcmu(csvec,a,cslow,csup,nu_nucl,mu_crit,t);
ecd = io.*sqrt(reshape(rxncmat',numpart*disc.steps,1)).*sqrt(exp(muvec)).*(1-csvec);
eta = muvec-reshape(rxnphimat',numpart*disc.steps,1);
val(disc.sol:end-1) = ecd.*(exp(-alpha.*eta)-exp((1-alpha).*eta));
val(disc.sol:end-1) = val(disc.sol:end-1) + interp1q(tr,noise,t)';

% CURRENT CONDITION
val(end) = currset;
val = real(val);

return;

function mu = calcmu(cs,a,cslow,csup,nu_nucl,mu_crit,t)
% This function calculates the chemical potential. If we're
% passed the low-concentration spinodal point but below the
% high-concentration binodal, we're phase separated, with mu = 0.
% Alternatively, if a particle has already nucleated and still
% below csup, mu = 0.

global nucl_vec t_last
% Establish where the particles are
gtr_than_cslow = (cslow < cs); % 1 if cs > cslow
% "if we've passed cslow and not nucleated, nucleate that
% particle"
% it has _not_ nucleated if nucl_vec is 1
nucleated_spinod = bitand(gtr_than_cslow, nucl_vec);
nucl_vec = nucl_vec - nucleated_spinod;
% Figure out the probability that particles nucleated in the last
% function call (delta t will go to zero unless this is one in
% which a real time step was taken)
delta_t = t-t_last;
mu_homog = log(cs./(1-cs))+a.*(1-2.*cs);
nucl_rate = nu_nucl*exp(-(mu_crit-mu_homog)); % -> 0 for low cs
% 0 percent chance to nucleate if already have nuclated (nucl_vec=0)
% 1 * Poisson distribution probability of nucleating otherwise
prob_vec = nucl_vec.*(1 - exp(-delta_t*nucl_rate));
% Generate some uniform random numbers between 0 and 1 to compare to.
rand_vec = rand(max(size(cs)),1);
% If the random number is lower than the probability, the
% particle nucleated.
nucleated_rand = (rand_vec < prob_vec);
nucl_vec = nucl_vec - nucleated_rand;
% If we've gone beyond csup, we're homogeneous again
gtr_than_csup = (csup < cs); % 1 if cs > csup
mumask = nucl_vec + gtr_than_csup;
mu = mu_homog.*mumask;

return;

function M = genMass(disc,poros,Nx,epsbeta,tp,pvolvec)

% Initialize
M = sparse(disc.len,disc.len);

% Electrolyte terms
M(1:disc.ss,1:disc.ss) = speye(disc.ss);
M(disc.ss+1:disc.ss+disc.steps,disc.ss+1:disc.ss+disc.steps) = poros*speye(disc.steps);

% Mass conservation between electrolyte and solid particles
numpart = max(size(pvolvec))/Nx;
for i=1:Nx
    for j=0:numpart-1
        M(disc.ss+i, disc.sol+(i-1)*numpart+j) = ...
                epsbeta*(1-tp)*pvolvec((i-1)*numpart+j+1,1) ...
                / sum(pvolvec((i-1)*numpart+1:i*numpart));
    end
end

% Potential terms
for i=1:Nx
    for j=0:numpart-1
        M(2*disc.ss+disc.steps+i, disc.sol+(i-1)*numpart+j) = ...
                epsbeta*pvolvec((i-1)*numpart+j+1,1) ...
                / sum(pvolvec((i-1)*numpart+1:i*numpart));
    end
end

% Solid particles
M(disc.sol:end-1,disc.sol:end-1) = speye(Nx*numpart,Nx*numpart);

% Current conservation
for i=1:Nx
    for j=0:numpart-1
        M(end,disc.sol+(i-1)*numpart+j) = ...
                pvolvec((i-1)*numpart+j+1,1) ...
                / sum(pvolvec((i-1)*numpart+1:i*numpart))/Nx;
    end
end

return;

function [value, isterminal, direction] = events(t,cpcs,io,currset,a,alpha,porosvec,numpart,...
                 Nx,disc,tp,zp,zm,nDp,nDm,tr,epsbeta,cslow,csup,nu_nucl,mu_crit,...
                 ffend,noise)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)      

% Calculate the filling fraction 
ffvec = sum(cpcs(disc.sol:end-1))/(Nx*numpart);
value = ffvec - ffend;
isterminal = 1;
direction = 0;
                        
return;
