function [t,cs,ffvec,vvec] = acr_sp_swse_rev1(currset,partsize,csinit,ffend)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
k = 1.381e-23;      % Boltzmann constant
T = 298;            % Temp, K
e = 1.602e-19;      % Charge of proton, C
Na = 6.02e23;       % Avogadro's number
F = e*Na;           % Faraday's number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DIMENSIONAL VALUES HERE
init_voltage = 3.5;                 % Initial cell voltage, V

% Particle size
part_size = partsize * 10^-9;                  % Average particle size, m

% Material properties
dim_a = 1.8560e-20;                 % Regular solution parameter, J
dim_kappa = 5.0148e-10;             % Gradient penalty, J/m
dim_b = 0.1916e9;                   % Stress, Pa
rhos = 1.3793e28;                   % site density, 1/m^3
csmax = rhos/Na;                    % maximum concentration, mol/m^3
cwet = 0.98;                        % Dimensionless wetted conc.
wet_thick = 2e-9;                   % Thickness of wetting on surf.
Vstd = 3.422;                       % Standard potential, V
alpha = 0.5;                        % Charge transfer coefficient

% Discretization settings
solid_disc = 1e-9;                 % Discretization size of solid, m (MUST BE LESS THAN LAMBDA)
tsteps = 200;                      % Number disc. in time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we take care of our particle size distributions
part_steps = ceil(part_size/solid_disc);

% Now we calculate the dimensionless quantities used in the simulation
if (currset ~= 0)
    tr = linspace(0,1/abs(currset),tsteps);
else
    tr = linspace(0,300,tsteps);
end
noise = .001*randn(max(size(tr)),part_steps);

% Material properties
kappa = dim_kappa ./ (k*T*rhos*part_size.^2);
a = dim_a / (k*T);
b = dim_b / (k*T*rhos);

% First we need to find the equilibrium solid concentration profiles for
% each particle - use a parfor loop to loop over the particles
opt = optimset('Display','Off');

if (csinit == 0)
    cs0base = 0.01;                 % Guess an initial filling fraction
    phi_init = -(init_voltage-Vstd)*(e/(k*T));
    cinit = 1;
    cs0init = cs0base*ones(part_steps,1);
    cs0 = fsolve(@calc_dcs_dt,cs0init,opt,phi_init,cinit,1,kappa,...
                                a,b,part_steps,alpha,cwet);

    % Assemble it all
    csinit = zeros(part_steps+1,1);
    csinit(1:part_steps) = real(cs0);
    csinit(end) = phi_init;
end
    
% Before we can call the solver, we need a Mass matrix
M = genMass(part_steps);

% Prepare to call the solver
% options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none');
options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none','Events',@events);
[t,cs]=ode15s(@calcRHS,tr,csinit,options,kappa,a,b,part_steps,alpha,cwet,currset,noise,tr,ffend);              
                
% First we calculate the voltage                 
vvec = Vstd - (k*T/e)*cs(:,end);

% Now the filling fraction vector - we only care about the active parts of
% the particle.  That is, we ignore the surface wetting as it does not move
ffvec = zeros(max(size(t)),1);
for i=1:max(size(t))
    ffvec(i) = ffvec(i) + ...
            sum(cs(i,1:part_steps)/part_steps);
end

return;

function val = calcRHS(t,cs,kappa,a,b,part_steps,alpha,cwet,currset,noise,tr,ffend)

val = zeros(part_steps+1,1);
val(1:part_steps) = calc_dcs_dt(cs(1:part_steps),cs(end),1,1, ...
                        kappa,a,b,part_steps,alpha,cwet) + interp1q(tr,noise,t)';
val(end) = real(currset);

return;

function dcsdt = calc_dcs_dt(cs,phi,ccath,io,kappa,a,b,ssteps,alpha,cwet)

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
ecd = io.*ccath^(1-alpha).*act.^(alpha).*(1-cs);
eta = mu-phi;
dcsdt = ecd.*(exp(-alpha.*eta)-exp((1-alpha).*eta));

return;

function M = genMass(part_steps)
    
M = sparse(part_steps+1,part_steps+1);
M(1:part_steps,1:part_steps) = speye(part_steps);
M(end,1:part_steps) = 1/part_steps;

return;


function [value, isterminal, direction] = events(t,cs,kappa,a,b, ...
                        part_steps,alpha,cwet,currset,noise,tr,ffend)
                        
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)      

% Calculate the filling fraction
ffvec = sum(cs(1:end-1))/part_steps;
value = ffvec - ffend;
isterminal = 1;
direction = 0;

                        
return;

