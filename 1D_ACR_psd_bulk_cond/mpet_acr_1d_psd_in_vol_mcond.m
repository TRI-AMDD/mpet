function [t,cpcs,csmat,ffvec,vvec,disc,part] = mpet_acr_1d_psd_in_vol_mcond_rev2b(dim_crate,part,disc,cpcs0)

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
% dim_crate = 1;                  % C-rate (electrode capacity per hour)
dim_io = .1;                      % Exchange current density, A/m^2 (0.1 for H2/Pt)

% Pick starting voltage based on charge/discharge
if dim_crate > 0
    init_voltage = 3.45;
else
    init_voltage = 3.35;                 
end

% Electrode properties
Lx = 24e-6;                         % electrode thickness, m
Lsep = 300e-6;                        % separator thickness, m
% Asep = 1e-4;                        % area of separator, m^2
Lp = 0.69;                          % Volume loading percent active material
poros = 0.44;                        % Porosity
c0 = 1000;                          % Initial electrolyte conc., mol/m^3
zp = 1;                             % Cation charge number
zm = 1;                             % Anion charge number
Dp = 2.2e-10;                       % Cation diff, m^2/s, LiPF6 in EC/DMC
Dm = 2.94e-10;                      % Anion diff, m^2/s, LiPF6 in EC/DMC
Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm);   % Ambipolar diffusivity
tp = zp*Dp / (zp*Dp + zm*Dm);       % Cation transference number
mcond = .001;                      % Dimensional electronic conductivity
ASRcont = 20e-4;  %24                     % Dimensional area specific contact resistance, ohm*m^2


% Particle size distribution
mean = 208e-9;%152e-9;                      % Average particle size, m
stddev = 73e-9;%57e-9;                     % Standard deviation, m

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
alpha = 0.18;                        % Charge transfer coefficient

% Discretization settings
Nx = 8;                            % Number disc. in x direction
solid_disc = 2e-9;                 % Discretization size of solid, m (MUST BE LESS THAN LAMBDA)
numpart = 5;                       % Particles per volume
tsteps = 200;                      % Number disc. in time
ffend = 0.98;                       % Final filling fraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we take care of our particle size distributions
totalpart = Nx*numpart;
wet_steps = ceil(wet_thick / solid_disc);
if ~isa(part,'struct')
    part_size = abs(normrnd(mean,stddev,totalpart,1));
    part_steps = [0; ceil(part_size/solid_disc)];
    pareavec = (4*pi).*part_size.^2;
    pvolvec = (4/3).*pi.*part_size.^3;
%    ap = 3.475 ./ part_size;     % Area:volume ratio for plate particles
    part = struct('steps',part_steps,'sizes',part_size, ...
                'wetsteps',wet_steps,'cwet',cwet, ...
                'areas', pareavec, 'vols', pvolvec);
elseif (isa(part,'struct') && (max(size(part.sizes)) == totalpart))
    part_size = part.sizes;
    part_steps = part.steps;
%    ap = 3.475 ./ part_size;    
    pareavec = part.areas;
    pvolvec = part.vols;
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
ASRcont = ASRcont * Lx * csmax * F * (1-eps) * Lp / td * e / (k * T);

% scaling of conductivity
mcond =  mcond *  (td * k * Na * T) / (Lx^2 * F^2 * c0);

if currset ~= 0
    tr = linspace(0,1/abs(currset),tsteps);
else    
    tr = linspace(0,30,100);
end
io = ((pareavec./pvolvec) .* dim_io .* td) ./ (F .* csmax);
epsbeta = ((1-poros)*Lp*csmax) / c0; % Vs/V*csmax/c0 = poros*beta

% Material properties
kappa = dim_kappa ./ (k*T*rhos*part_size.^2);
a = dim_a / (k*T);
b = dim_b / (k*T*rhos);


if (isa(disc,'struct'))
    % disc is a struct, and present, do nothing
    Nx = disc.Nx;
%    Ny = disc.Ny;
    numpart = disc.numpart
    sf = disc.sf;
    ssx = ceil(sf*Nx);
else
    sf = Lsep/Lx; % separator fraction
    ssx = ceil(sf*Nx);
    ss = ssx; % separator steps?
    disc = struct('ss',ss,'steps',Nx,'ssteps',sum(part_steps),...
                'len',2*(ss+Nx)+sum(part_steps)+Nx+1, ...
                'sol',2*(ss+Nx)+1,'Nx',Nx,'numpart',numpart,'sf',sf);
end


if max(size(cpcs0)) == 1
    % First we need to find the equilibrium solid concentration profiles for
    % each particle - use a parfor loop to loop over the particles
    opt = optimset('Display','Off');
    outarr = cell(totalpart,1);       % cell array for output
    % ME: set initial condition depending on current direction
    if dim_crate > 0
        cs0 = 0.01;                 % Guess an initial filling fraction
    else
        cs0 = 0.98;
        ffend = 1 - ffend;
    end
    
    phi_init = -(init_voltage-Vstd)*(e/(k*T));
    cinit = 1;
    disp('Generating initial solid concentration profiles...')
    parfor i=1:totalpart
        % Use fsolve to get solid concentration profiles
        csinit = cs0*ones(part_steps(i+1),1);
        outarr{i} = csinit;
%         outarr{i} = fsolve(@calc_dcs_dt,csinit,opt,phi_init,cinit,io(i),kappa(i),...
%                                     a,b,part_steps(i+1),alpha,cwet,wet_steps);
    end

    % Now we move our solid concentrations vectors, appending to
    % the full vector, initcsvec
    initcsvec = real(outarr{1});
    for i=2:totalpart
       initcsvec = [initcsvec; real(outarr{i})];
    end
    disp('Done.')

    % Create the initial vector for the time stepper
    disp('Generating initial vector...')
    
    % Assemble it all
    cpcsinit = zeros(disc.len,1);
    cpcsinit(1:disc.ss+disc.steps) = cinit;
    cpcsinit(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = phi_init;
    cpcsinit(disc.sol:disc.sol+disc.ssteps-1) = initcsvec;
    cpcsinit(disc.sol+disc.ssteps:end-1) = 0;
    cpcsinit(end) = phi_init;
else
    cpcsinit = cpcs0;
end

% Create a grid for the porosity
%pxg = ones(Ny,ssx+Nx+1);
%pyg = ones(Ny+1,ssx+Nx);
%pxg(:,ssx+1:end) = poros;
%pyg(:,ssx+1:end) = poros;
porosvec = ones(disc.ss+disc.steps+1, 1);
porosvec(disc.ss+1:end) = poros;
porosvec = porosvec.^(3/2);     % Bruggeman

%% Bruggeman relation
%pxg = pxg.^(3/2);
%pyg = pyg.^(3/2);

% Before we can call the solver, we need a Mass matrix
M = genMass(disc,poros,part_steps,Nx,epsbeta,tp,pvolvec);

% Prepare to call the solver
% options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none');
options=odeset('Mass',M,'MassSingular','yes','MStateDependence','none','Events',@events);
disp('Calling ode15s solver...')
[t,cpcs]=ode15s(@calcRHS,tr,cpcsinit,options,io,currset,kappa,a,b,alpha, ...
                    cwet,wet_steps,part_steps,Nx,disc,tp,zp,zm,nDp,nDm,mcond,porosvec,pvolvec,tr,epsbeta,ffend);

% Now we analyze the results before returning
disp('Done.')                
disp('Calculating the voltage and filling fraction vectors...')                
                
% First we calculate the voltage                 
vvec = Vstd - (k*T/e)*cpcs(:,end) - (k*T/e) * currset * ASRcont;

% Converting scaled time to actual time
t = t.*td;

% Now the filling fraction vector - we only care about the active parts of
% the particle.  That is, we ignore the surface wetting as it does not move
ffvec = zeros(max(size(t)),1);
for i=1:max(size(t))
    for j=1:Nx
        for k=0:numpart-1
            ind1 = sum(part_steps(1:(j-1)*numpart+1+k))+1;
            ind2 = sum(part_steps(1:(j-1)*numpart+1+k+1));
%                    sum(cpcs(i,disc.sol+ind1-1:disc.sol+ind2-1))
            ffvec(i) = ffvec(i) + ...
                    sum(cpcs(i,disc.sol+ind1:disc.sol+ind2)) ...
                    / (Nx*part_steps((j-1)*numpart+1+k+1)) ...
                    * pvolvec((j-1)*numpart+k+1,1) ...
                    / sum(pvolvec((j-1)*numpart+1:j*numpart));
        end
%        ind1 = sum(part_steps(1:j))+1;
%        ind2 = sum(part_steps(1:j+1));
%        ffvec(i) = ffvec(i) + ...
%                sum(cpcs(i,disc.sol+ind1-1:disc.sol+ind2-1))/(Nx*part_steps(j+1));
    end
end

% Create cs matrix
csmat = zeros(max(size(t)), Nx, numpart);
for i=1:max(size(t))
    for j=1:Nx
        for k=0:numpart-1
            ind1 = sum(part_steps(1:(j-1)*numpart+1+k))+1;
            ind2 = sum(part_steps(1:(j-1)*numpart+1+k+1));
            csmat(i,j,k+1) = ...
                    sum(cpcs(i,disc.sol+ind1-1:disc.sol+ind2-1)) ...
                    / (part_steps((j-1)*numpart+1+k+1));
        end
    end
end

disp('Finished.')

return;

function val = calcRHS(t,cpcs,io,currset,kappa,a,b,alpha,cwet,wet_steps,...
                 part_steps,Nx,disc,tp,zp,zm,nDp,nDm,mcond,porosvec,pvolvec,tr,epsbeta,ffend)

% Initialize output
val = zeros(max(size(cpcs)),1);             
             
% Pull out the concentrations first
cvec = cpcs(1:disc.ss+disc.steps);
phivec = cpcs(disc.ss+disc.steps+1:2*(disc.ss+disc.steps));
phi0 = cpcs(end);

% MASS CONSERVATION - ELECTROLYTE DIFFUSION
%cxtmp = zeros(Ny,ssx+Nx+2);
%cytmp = zeros(Ny+2,ssx+Nx);
%cxtmp(:,2:end-1) = reshape(cvec,Ny,ssx+Nx);
%cytmp(2:end-1,:) = reshape(cvec,Ny,ssx+Nx);
ctmp = zeros(disc.ss+disc.steps+2,1);
ctmp(2:end-1) = cvec;
% No flux conditions
ctmp(end) = ctmp(end-1);
%cxtmp(:,end) = cxtmp(:,end-1);
%cytmp(1,:) = cytmp(2,:);
%cytmp(end,:) = cytmp(end-1,:);
% Flux into separator
%cxtmp(:,1) = cxtmp(:,2) + currset*beta*(1-tp)/Nx;
ctmp(1) = ctmp(2) + currset*epsbeta*(1-tp)/Nx;
% Get fluxes, multiply by porosity
cflux = -porosvec.*diff(ctmp).*Nx;
%cxflux = -pxg.*diff(cxtmp,1,2).*Nx;
%cyflux = -pyg.*diff(cytmp,1,1).*Ny;
% Divergence of the fluxes
val(1:disc.ss+disc.steps) = -diff(cflux).*Nx;
%val(1:disc.ss+disc.steps) = reshape(-diff(cxflux,1,2)*Nx - ...
%                                diff(cyflux,1,1)*Ny,disc.ss+disc.steps,1);
                           
% CHARGE CONSERVATION - DIVERGENCE OF CURRENT DENSITY
%phixtmp = zeros(Ny,ssx+Nx+2);
%phiytmp = zeros(Ny+2,ssx+Nx);
phitmp = zeros(disc.ss+disc.steps+2,1);
%phixtmp(:,2:end-1) = reshape(phivec,Ny,ssx+Nx);
%phiytmp(2:end-1,:) = reshape(phivec,Ny,ssx+Nx);
phitmp(2:end-1) = phivec;
% No flux conditions
phitmp(end) = phitmp(end-1);
%phixtmp(:,end) = phixtmp(:,end-1);
%phiytmp(1,:) = phiytmp(2,:);
%phiytmp(end,:) = phiytmp(end-1,:);
% Potential BC
phitmp(1) = phi0;
%phixtmp(:,1) = phi0;
% Average c values for current density on boundaries
cavg = (ctmp(1:end-1)+ctmp(2:end))/2;
%cavgx = (cxtmp(:,1:end-1)+cxtmp(:,2:end))/2;
%cavgy = (cytmp(1:end-1,:)+cytmp(2:end,:))/2;
currdens = -((zp*nDp-zp*nDm).*diff(ctmp).*Nx)- ...
            ((zp*nDp+zm*nDm).*cavg.*diff(phitmp).*Nx);
%cdx = -((zp*nDp-zp*nDm).*diff(cxtmp,1,2).*Nx)- ...
%            ((zp*nDp+zm*nDm).*cavgx.*(diff(phixtmp,1,2).*Nx));
%cdx = cdx.*pxg;
%cdy = -((zp*nDp-zp*nDm).*diff(cytmp,1,1).*Ny)- ...
%            ((zp*nDp+zm*nDm).*cavgy.*(diff(phiytmp,1,1).*Ny));
%cdy = cdy.*pyg;
val(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = ...
        -diff(porosvec.*currdens).*Nx;
%val(disc.ss+disc.steps+1:2*(disc.ss+disc.steps)) = ...
%        reshape(-diff(cdx,1,2).*Nx-diff(cdy,1,1).*Ny,disc.ss+disc.steps,1);

% Calculate the reaction rate of the particles
totalpart = max(size(io));
numpart = totalpart/Nx;
outarr = cell(totalpart,1);
for i=1:Nx
    % Get the electrolyte concentration and potential values
    ci = cpcs(disc.ss+i);
    phii = cpcs(2*disc.ss+disc.steps+i);
    phimi = cpcs(disc.sol+disc.ssteps+i-1);   % ME: check if indices are correct -> should be ok now...
    % Loop through particles in each volume
    for j=0:numpart-1
        ind1 = sum(part_steps(1:(i-1)*numpart+1+j))+1;
        ind2 = sum(part_steps(1:(i-1)*numpart+1+j+1));
        csi = cpcs(disc.sol+ind1-1:disc.sol+ind2-1);
        outarr{(i-1)*numpart+j+1} = calc_dcs_dt(csi,phii,phimi,ci, ...
                io((i-1)*numpart+j+1), ...
                kappa((i-1)*numpart+j+1),a,b, ...
                part_steps((i-1)*numpart+1+j+1), ...
                alpha, cwet, wet_steps);
    end
%    ind1 = sum(part_steps(1:i))+1;
%    ind2 = sum(part_steps(1:i+1));
%    csi = cpcs(disc.sol+ind1-1:disc.sol+ind2-1);
%    outarr{i} = calc_dcs_dt(csi,phii,ci,io(i),kappa(i),...
%                                a,b,part_steps(i+1),alpha,cwet,wet_steps);
end

% Move to output vector
val(disc.sol:disc.sol+disc.ssteps-1) = real(cell2mat(outarr));

% Charge conservation for phim domain

% % original setup
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Construct our charge conservation the electron conducting phase
% phimtmp = zeros(2+disc.steps,1);
% phimtmp(2:end-1) = cpcs(disc.sol+disc.ssteps:end-1);
% % No flux condition
% phimtmp(1) = phimtmp(2);    
% 
% % Robin condition - set ground and flux out
% phimtmp(end-1) = currset/mcond/Nx; 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% new test
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Construct our charge conservation the electron conducting phase
phimtmp = zeros(2+disc.steps,1);

% Robin condition - set ground and flux out
% phimtmp(end-1) = currset / mcond / Nx; 

phimtmp(2:end-1) = cpcs(disc.sol+disc.ssteps:end-1);

% No flux condition
phimtmp(1) = phimtmp(2); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



val(disc.sol+disc.ssteps:end-1) = -diff(-mcond.*diff(phimtmp)*Nx).*Nx;
% Finally the current condition
val(end) = currset;

return;

function dcsdt = calc_dcs_dt(cs,phi,phim,ccath,io,kappa,a,b,ssteps,alpha,cwet,wet_steps)

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
eta = mu-phi+phim;                                      % ME: check the sign...
dcsdt = ecd.*(exp(-alpha.*eta)-exp((1-alpha).*eta));

return;

function M = genMass(disc,poros,part_steps,Nx,epsbeta,tp,pvolvec)
    
M = sparse(disc.len,disc.len);
% Electrolyte terms
M(1:disc.ss,1:disc.ss) = speye(disc.ss);
M(disc.ss+1:disc.ss+disc.steps,disc.ss+1:disc.ss+disc.steps) = poros*speye(disc.steps);

% Mass conservation (solid particles)
numpart = max(size(pvolvec))/Nx;
for i=1:Nx
%    ind1 = sum(part_steps(1:i))+1;
%    ind2 = sum(part_steps(1:i+1));
%    M(disc.ss+i, disc.sol+ind1:disc.sol+ind2) = (beta*(1-tp))/part_steps(i+1);
    % Loop through particles in each volume
    for j=0:numpart-1
        ind1 = sum(part_steps(1:(i-1)*numpart+1+j))+1;
        ind2 = sum(part_steps(1:(i-1)*numpart+1+j+1));
        M(disc.ss+i, disc.sol+ind1:disc.sol+ind2) = ...
                (epsbeta*(1-tp))/part_steps((i-1)*numpart+1+j+1) ...
                * pvolvec((i-1)*numpart+j+1,1) ...
                / sum(pvolvec((i-1)*numpart+1:i*numpart));
    end
end

% Potential terms
for i=1:Nx
%    ind1 = sum(part_steps(1:i))+1;
%    ind2 = sum(part_steps(1:i+1));
%    M(2*disc.ss+disc.steps+i, disc.sol+ind1:disc.sol+ind2) = beta/part_steps(i+1);
    % Loop through particles
    for j=0:numpart-1
        ind1 = sum(part_steps(1:(i-1)*numpart+1+j))+1;
        ind2 = sum(part_steps(1:(i-1)*numpart+1+j+1));
        M(2*disc.ss+disc.steps+i, disc.sol+ind1:disc.sol+ind2) = ...
                epsbeta/part_steps((i-1)*numpart+1+j+1) ...
                * pvolvec((i-1)*numpart+j+1,1) ...
                / sum(pvolvec((i-1)*numpart+1:i*numpart));
    end
end

% Solid particles
M(disc.sol:disc.sol+disc.ssteps-1,disc.sol:disc.sol+disc.ssteps-1) = speye(sum(part_steps));

% Potential terms
for i=1:Nx
%    ind1 = sum(part_steps(1:i))+1;
%    ind2 = sum(part_steps(1:i+1));
%    M(end, disc.sol+ind1:disc.sol+ind2) = 1/(part_steps(i+1)*Nx*Ny);
    % Loop through particles
    for j=0:numpart-1
        ind1 = sum(part_steps(1:(i-1)*numpart+1+j))+1;
        ind2 = sum(part_steps(1:(i-1)*numpart+1+j+1));
        M(end, disc.sol+ind1:disc.sol+ind2) = ...
        pvolvec((i-1)*numpart+j+1,1) ...
        / sum(pvolvec((i-1)*numpart+1:i*numpart)) ...
        / (part_steps((i-1)*numpart+1+j+1)*Nx);
    end
end

% Electric potential
for i=1:Nx
    % Loop through particles
    for j=0:numpart-1
        ind1 = sum(part_steps(1:(i-1)*numpart+1+j))+1;
        ind2 = sum(part_steps(1:(i-1)*numpart+1+j+1));
        M(disc.sol+disc.ssteps+i-1, disc.sol+ind1:disc.sol+ind2) = ...
                -epsbeta./part_steps((i-1)*numpart+1+j+1) ...
                * pvolvec((i-1)*numpart+j+1,1) ...
                / sum(pvolvec((i-1)*numpart+1:i*numpart));      
    end
end

return;

function [value, isterminal, direction] = events(t,cpcs,io,currset,kappa,a,b,alpha,cwet,wet_steps,...
                 part_steps,Nx,disc,tp,zp,zm,nDp,nDm,mcond,porosvec,pvolvec,tr,epsbeta,ffend)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)      

% Calculate the filling fraction
%totalpart = max(size(io));
%numpart = totalpart/Nx;
numpart = max(size(pvolvec))/Nx;
ffvec = 0;
for j=1:Nx
        for k=0:numpart-1
            ind1 = sum(part_steps(1:(j-1)*numpart+1+k))+1;
            ind2 = sum(part_steps(1:(j-1)*numpart+1+k+1));
%            sum(cpcs(disc.sol+ind1:disc.sol+ind2))
%            (Nx*part_steps((j-1)*numpart+1+k+1))
%            pvolvec((j-1)*numpart+k+1,1)
%            sum(pvolvec((j-1)*numpart+1:j*numpart))
            ffvec = ffvec + ...
                    sum(cpcs(disc.sol+ind1:disc.sol+ind2)) ...
                    / (Nx*part_steps((j-1)*numpart+1+k+1)) ...
                    * pvolvec((j-1)*numpart+k+1,1) ...
                    / sum(pvolvec((j-1)*numpart+1:j*numpart));
%                    sum(cpcs(disc.sol+ind1-1:disc.sol+ind2-1))/(Nx*Ny*part_steps(j+1));
        end
%        ind1 = sum(part_steps(1:j))+1;
%        ind2 = sum(part_steps(1:j+1));
%        ffvec = ffvec + ...
%                sum(cpcs(disc.sol+ind1-1:disc.sol+ind2-1))/(Nx*Ny*part_steps(j+1));
end
value = ffvec - ffend;
isterminal = 1;
direction = 0;

return;
