function [t,cpcs,ffvec,vvec,disc,part] = acr_mpet_swcs_rev2(dim_crate,part,cpcs0,disc,ffend)

% NAME: Allen-Cahn surface wetted coherent strain particles in a porous electrode
% AUTHOR: TR Ferguson, Bazant Group, MIT
% DATE: 4/30/2013

% This script simulates a porous electrode at constant current.  The
% particles are assumed to be 1D, with surfaces (parallel to direction of
% intercalation) wetted.  The particle sizes are set by a normal
% distribution.

% OUTPUTS
% t = time vector 
% cpcs = conc., potential, and solid conc. vector
% ffvec = filling fraction vector
% vvec = voltage vector
% disc = discretization settings (struct)
% part = particle sizes and discretizations (struct)

% INPUTS
% dim_crate = dimensional C-rate
% part = particle sizes and discretizations (if from another simulation)
%           this allows you to run the same setup at a different C-rate               

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
% dim_crate = 0.1;                  % C-rate (electrode capacity per hour)
dim_io = 0.1;                      % Exchange current density, A/m^2 (0.1 for H2/Pt)

% Pick starting voltage based on charge/discharge
if dim_crate > 0
    init_voltage = 3.45;
else
    init_voltage = 3.35;                 
end

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
dim_kappa = 5.0148e-10;             % Gradient penalty, J/m
dim_b = 0.1916e9;                   % Stress, Pa
% dim_b = 0; dim_kappa = 0;

rhos = 1.3793e28;                   % site density, 1/m^3
csmax = rhos/Na;                    % maximum concentration, mol/m^3
cwet = 0.98;                        % Dimensionless wetted conc.
wet_thick = 2e-9;                   % Thickness of wetting on surf.
Vstd = 3.422;                       % Standard potential, V
alpha = 0.5;                        % Charge transfer coefficient

% Discretization settings
Nx = 8;                            % Number disc. in x direction
Ny = 2;                            % Number disc. in y direction
solid_disc = 1e-9;                 % Discretization size of solid, m (MUST BE LESS THAN LAMBDA)
tsteps = 200;                      % Number disc. in time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we take care of our particle size distributions
num_particles = Nx * Ny;
wet_steps = ceil(wet_thick / solid_disc);
if ~isa(part,'struct')
    part_size = abs(normrnd(mean,stddev,num_particles,1));
    part_steps = [0; ceil(part_size/solid_disc)];
    ap = 3.475 ./ part_size;     % Area:volume ratio for plate particles
    part = struct('steps',part_steps,'sizes',part_size, ...
                'wetsteps',wet_steps,'cwet',cwet);
elseif (isa(part,'struct') && (max(size(part.sizes)) == Nx * Ny))
    part_size = part.sizes;
    part_steps = part.steps;
    ap = 3.475 ./ part_size;    
    wet_steps = part.wetsteps;
    cwet = part.cwet;
else
    error('Error retrieving particle sizes. Be sure your input values match the discretization settings') 
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
io = (ap .* dim_io .* td) ./ (F .* csmax);
beta = ((1-poros)*Lp*csmax) / (poros*c0);

% Material properties
kappa = dim_kappa ./ (k*T*rhos*part_size.^2);
a = dim_a / (k*T);
b = dim_b / (k*T*rhos);


if (isa(disc,'struct'))
    % disc is a struct, and present, do nothing
    Nx = disc.Nx;
    Ny = disc.Ny;
    sf = disc.sf;
    ssx = ceil(sf*Nx);
else
    sf = Lsep/Lx;
    ssx = ceil(sf*Nx);
    ss = ssx*Ny;
    disc = struct('ss',ss,'steps',Nx*Ny,'ssteps',sum(part_steps),...
                'len',2*(ss+Nx*Ny)+sum(part_steps)+1, ...
                'sol',2*(ss+Nx*Ny)+1,'Nx',Nx,'Ny',Ny,'sf',sf);
end


if max(size(cpcs0)) == 1
    % First we need to find the equilibrium solid concentration profiles for
    % each particle - use a parfor loop to loop over the particles
    opt = optimset('Display','Off');
    outarr = cell(Nx*Ny,1);       % cell array for output
    cs0 = 0.01;                 % Guess an initial filling fraction
    phi_init = -(init_voltage-Vstd)*(e/(k*T));
    cinit = 1;
    disp('Generating initial solid concentration profiles...')
    parfor i=1:Nx*Ny
        % Use fsolve to get solid concentration profiles
        csinit = cs0*ones(part_steps(i+1),1);
        outarr{i} = csinit;
%         outarr{i} = fsolve(@calc_dcs_dt,csinit,opt,phi_init,cinit,io(i),kappa(i),...
%                                     a,b,part_steps(i+1),alpha,cwet,wet_steps);
    end

    % Now we move our solid concentrations vectors
    initcsvec = real(outarr{1});
    for i=2:Nx*Ny
       initcsvec = [initcsvec; real(outarr{i})];
    end
    disp('Done.')

    % Create the initial vector for the time stepper
    disp('Generating initial vector...')
    
    % Assemble it all
    cpcsinit = zeros(disc.len,1);
    cpcsinit(1:disc.ss+disc.steps) = cinit;
    cpcsinit(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = phi_init;
    cpcsinit(disc.sol:end-1) = initcsvec;
    cpcsinit(end) = phi_init;
else
    cpcsinit = cpcs0;
end

% Create a grid for the porosity
pxg = ones(Ny,ssx+Nx+1);
pyg = ones(Ny+1,ssx+Nx);
pxg(:,ssx+1:end) = poros;
pyg(:,ssx+1:end) = poros;

% Bruggeman relation
pxg = pxg.^(3/2);
pyg = pyg.^(3/2);

% Before we can call the solver, we need a Mass matrix
M = genMass(disc,poros,part_steps,Nx,Ny,beta,tp);

% Prepare to call the solver
% options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none');
options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none','Events',@events);
disp('Calling ode15s solver...')
[t,cpcs]=ode15s(@calcRHS,tr,cpcsinit,options,io,currset,kappa,a,b,alpha, ...
                    cwet,wet_steps,part_steps,Nx,Ny,ssx,disc,tp,zp,zm,nDp,nDm,pxg,pyg,tr,beta,ffend);

% Now we analyze the results before returning
disp('Done.')                
disp('Calculating the voltage and filling fraction vectors...')                
                
% First we calculate the voltage                 
vvec = Vstd - (k*T/e)*cpcs(:,end);

% Now the filling fraction vector - we only care about the active parts of
% the particle.  That is, we ignore the surface wetting as it does not move
ffvec = zeros(max(size(t)),1);
for i=1:max(size(t))
    for j=1:Nx*Ny
        ind1 = sum(part_steps(1:j))+1;
        ind2 = sum(part_steps(1:j+1));
        ffvec(i) = ffvec(i) + ...
                sum(cpcs(i,disc.sol+ind1-1:disc.sol+ind2-1))/(Nx*Ny*part_steps(j+1));
    end
end
disp('Finished.')

return;

function val = calcRHS(t,cpcs,io,currset,kappa,a,b,alpha,cwet,wet_steps,...
                 part_steps,Nx,Ny,ssx,disc,tp,zp,zm,nDp,nDm,pxg,pyg,tr,beta,ffend)

% Initialize output
val = zeros(max(size(cpcs)),1);             
             
% Pull out the concentrations first
cvec = cpcs(1:disc.ss+disc.steps);
phivec = cpcs(disc.ss+disc.steps+1:2*(disc.ss+disc.steps));
phi0 = cpcs(end);

% MASS CONSERVATION - ELECTROLYTE DIFFUSION
cxtmp = zeros(Ny,ssx+Nx+2);
cytmp = zeros(Ny+2,ssx+Nx);
cxtmp(:,2:end-1) = reshape(cvec,Ny,ssx+Nx);
cytmp(2:end-1,:) = reshape(cvec,Ny,ssx+Nx);
% No flux conditions
cxtmp(:,end) = cxtmp(:,end-1);
cytmp(1,:) = cytmp(2,:);
cytmp(end,:) = cytmp(end-1,:);
% Flux into separator
cxtmp(:,1) = cxtmp(:,2) + currset*beta*(1-tp)/Nx;
% Get fluxes, multiply by porosity
cxflux = -pxg.*diff(cxtmp,1,2).*Nx;
cyflux = -pyg.*diff(cytmp,1,1).*Ny;
% Divergence of the fluxes
val(1:disc.ss+disc.steps) = reshape(-diff(cxflux,1,2)*Nx - ...
                                diff(cyflux,1,1)*Ny,disc.ss+disc.steps,1);
                           
% CHARGE CONSERVATION - DIVERGENCE OF CURRENT DENSITY
phixtmp = zeros(Ny,ssx+Nx+2);
phiytmp = zeros(Ny+2,ssx+Nx);
phixtmp(:,2:end-1) = reshape(phivec,Ny,ssx+Nx);
phiytmp(2:end-1,:) = reshape(phivec,Ny,ssx+Nx);
% No flux conditions
phixtmp(:,end) = phixtmp(:,end-1);
phiytmp(1,:) = phiytmp(2,:);
phiytmp(end,:) = phiytmp(end-1,:);
% Potential BC
phixtmp(:,1) = phi0;
% Average c values for current density on boundaries
cavgx = (cxtmp(:,1:end-1)+cxtmp(:,2:end))/2;
cavgy = (cytmp(1:end-1,:)+cytmp(2:end,:))/2;
cdx = -((zp*nDp-zp*nDm).*diff(cxtmp,1,2).*Nx)- ...
            ((zp*nDp+zm*nDm).*cavgx.*(diff(phixtmp,1,2).*Nx));
cdx = cdx.*pxg;
cdy = -((zp*nDp-zp*nDm).*diff(cytmp,1,1).*Ny)- ...
            ((zp*nDp+zm*nDm).*cavgy.*(diff(phiytmp,1,1).*Ny));
cdy = cdy.*pyg;
val(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = ...
        reshape(-diff(cdx,1,2).*Nx-diff(cdy,1,1).*Ny,disc.ss+disc.steps,1);

% Calculate the reaction rate of the particles
outarr = cell(Nx*Ny,1);
for i=1:Nx*Ny
    ind1 = sum(part_steps(1:i))+1;
    ind2 = sum(part_steps(1:i+1));
    % Get the solid concentration and potential values
    csi = cpcs(disc.sol+ind1-1:disc.sol+ind2-1);
    ci = cpcs(disc.ss+i);
    phii = cpcs(2*disc.ss+disc.steps+i);
    outarr{i} = calc_dcs_dt(csi,phii,ci,io(i),kappa(i),...
                                a,b,part_steps(i+1),alpha,cwet,wet_steps);
end

% Move to output vector
val(disc.sol:end-1) = real(cell2mat(outarr));

% Finally the current condition
val(end) = currset;

return;

function dcsdt = calc_dcs_dt(cs,phi,ccath,io,kappa,a,b,ssteps,alpha,cwet,wet_steps)

% This function returns the reaction rate for each channel in the particle
% - first we need the curvature and the average filling, X
cstmp = zeros(max(size(cs))+2,1);
cstmp(2:end-1) = cs;
cstmp(1) = cwet;
cstmp(end) = cwet;
dxs = 1/(ssteps);
curv = diff(diff(cstmp))./(dxs^2);

meanfill = (sum(cs)) / (ssteps);
mu = log(cs./(1-cs)) + a.*(1-2.*cs) - kappa.*curv + b.*(cs-meanfill);
act = exp(mu);

% Standard finite volume transition state
gamma_ts = (1./(1-cs));
% New transition state with regular solution parameter
% gamma_ts = (1./(1-cs)).*exp(a.*(1-2.*cs));

ecd = (io.*ccath^(1-alpha).*act.^(alpha))./gamma_ts;
eta = mu-phi;
dcsdt = ecd.*(exp(-alpha.*eta)-exp((1-alpha).*eta));

return;

function M = genMass(disc,poros,part_steps,Nx,Ny,beta,tp)
    
M = sparse(disc.len,disc.len);
% Electrolyte terms
M(1:disc.ss,1:disc.ss) = speye(disc.ss);
M(disc.ss+1:disc.ss+disc.steps,disc.ss+1:disc.ss+disc.steps) = poros*speye(disc.steps);

% Mass conservation (solid particles)
for i=1:Nx*Ny
    % Loop through particles
    ind1 = sum(part_steps(1:i))+1;
    ind2 = sum(part_steps(1:i+1));
    M(disc.ss+i, disc.sol+ind1:disc.sol+ind2) = (beta*(1-tp))/part_steps(i+1);
end

% Potential terms
for i=1:Nx*Ny
    % Loop through particles
    ind1 = sum(part_steps(1:i))+1;
    ind2 = sum(part_steps(1:i+1));
    M(2*disc.ss+disc.steps+i, disc.sol+ind1:disc.sol+ind2) = beta/part_steps(i+1);
end

% Solid particles
M(disc.sol:end-1,disc.sol:end-1) = speye(sum(part_steps));

% Potential terms
for i=1:Nx*Ny
    % Loop through particles
    ind1 = sum(part_steps(1:i))+1;
    ind2 = sum(part_steps(1:i+1));
    M(end, disc.sol+ind1:disc.sol+ind2) = 1/(part_steps(i+1)*Nx*Ny);
end

return;

function [value, isterminal, direction] = events(t,cpcs,io,currset,kappa,a,b,alpha,cwet,wet_steps,...
                 part_steps,Nx,Ny,ssx,disc,tp,zp,zm,nDp,nDm,pxg,pyg,tr,beta,ffend)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)      

% Calculate the filling fraction
ffvec = 0;
for j=1:Nx*Ny
        ind1 = sum(part_steps(1:j))+1;
        ind2 = sum(part_steps(1:j+1));
        ffvec = ffvec + ...
                sum(cpcs(disc.sol+ind1-1:disc.sol+ind2-1))/(Nx*Ny*part_steps(j+1));
end
value = ffvec - ffend;
isterminal = 1;
direction = 0;

                        
return;
