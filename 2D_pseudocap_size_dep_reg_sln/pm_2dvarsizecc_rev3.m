function [t,cpcs,cmat,csmat,ffvec,vvec,disc,param,szs] = pm_2dvarsizecc_rev3(dim_crate,param,useparam,cpcs0,FFend)

% This script can simulate any direction (charge/discharge) and accepts a
% given starting point and ending filling fraction

% This function simulates a 2D pseudocapacitor cell with the ability to
% vary local particle size in each volume.  The current is modeled via the
% Modified Butler-Volmer Equation with a transfer coefficient of 0.5. Solid
% diffusion is not modeled, but will be in future models.

% This model assumes constant diffusivities in the electrolyte, which are
% given by the dilute limit diffusivities of the cation and anion.  The
% next revision will address CST

% DISCRETIZATION SETUP
%   The matrix will be represented as a vector, which is reshaped to a
%   matrix before gradient/divergence operations, then converted back to a
%   vector before passing to the ODE solver

% CONSTANTS
k = 1.381e-23;
e = 1.602e-19;
T = 298;
Na = 6.02e23;
F = Na*e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DIMENSIONAL VALUES HERE

% Discharge settings
% dim_crate = 0.1;                  % C-rate (electrode capacity per hour)
dim_io = 0.01;                      % Exchange current density, A/m^2 (0.1 for H2/Pt)
Vinit = 3.5;                        % Initial voltage, (V)

% Electrode properties
Lx = 50e-6;                         % electrode thickness, m
Lp = 0.69;                          % Volume loading percent active material
poros = 0.3;                        % Porosity
c0init = 1200;                          % Initial electrolyte conc., mol/m^3
zp = 1;                             % Cation charge number
zm = 1;                             % Anion charge number
Dp = 2.2e-10;                       % Cation diff, m^2/s, LiPF6 in EC/DMC
Dm = 2.94e-10;                      % Anion diff, m^2/s, LiPF6 in EC/DMC
Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm);   % Ambipolar diffusivity
tp = zp*Dp / (zp*Dp + zm*Dm);       % Cation transference number

% Particle size distribution
avg = 160e-9;                      % Average particle size, m
stddev = 30e-9;                    % Standard deviation, m

% Discretization settings
sf = 0.5;      % Length of separator as fraction of cathode
xs = 20;
ys = 10;
ss = sf*xs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Material properties
rhos = 1.3793e28;                   % site density, 1/m^3
csmax = rhos/Na;                    % maximum concentration, mol/m^3

% Generate our particle sizes from a normal distribution with a mean and
% std deviation of size
szs = abs(avg + stddev*randn(xs*ys,1));
ap = 3.475./(reshape(szs,xs*ys,1));

% Non-dimensional discharge settings
td = Lx^2 / Damb;       % Diffusive time
nDp = Dp / Damb;
nDm = Dm / Damb;
currset = dim_crate * (td/3600);
io = (ap .* dim_io .* td) ./ (F .* csmax);

% MATERIAL PROPERTIES 
Vo = 3.422;
ndVinit = (Vinit - Vo)*(e/(k*T));

% Convert sizes to nm units for rest of simulation
szs = szs/1e-9;

% MESH DIMENSIONS

tinit = 0;
if currset ~= 0
    tfinal = 1/abs(currset);
else
    tfinal = 1000;
end
tsteps = 200;
tr = linspace(tinit,tfinal,tsteps);
sl = ss*ys;
el = xs*ys;

% NON-DIMENSIONAL QUANTITIES
beta = ((1-poros)*Lp*csmax) / (poros*c0init);
betamat = beta*ones(ys,xs);

% GET OUR MASS MATRIX
M = genMassMatrix(ss,xs,ys,tp,betamat,poros);

% CONSTRUCT POROSITY MATRICES FOR ELECTRODE
%   Dx, Dy - diffusivity corrections, D_eff = D*poros^.5
%   Kx, Ky - conductivity corrections, K_eff = K*poros^1.5
xfmat = ones(ys,xs+ss+1);
yfmat = ones(ys+1,xs+ss);
% First Dx
Dxmat = xfmat;
Dxmat(:,ss+1:end) = poros^1.5;
% Now Dy
Dymat = yfmat;
Dymat(:,ss+1:end) = poros^1.5;
% Then insert into our struct
porosmats = struct('Dx',Dxmat,'Dy',Dymat);

% SETUP SOLVER AND PASS FUNCTION
options = odeset('MASS',M,'MassSingular','yes','MStateDependence','none','Events',@events);
%options = odeset('MASS',M,'Jacobian',@calcJ,'MassSingular','yes','MStateDependence','none');
%options = odeset('MASS',M,'Jacobian',@calcJ,'MassSingular','yes','MStateDependence','none','RelTol',1e-7);

% INITIALIZE VECTORS
c0 = c0init / 1000;
% These are guesses for the solver
if (currset>0)
    cs0 = .01;      % Start empty for discharge
else
    cs0 = .99;      % Start full for charging
end

% Now we generate a vector of regular solution parameters
if ~useparam
    param = size2regsoln(szs);
end
parammat = reshape(param,ys,xs);
% Now we generate initial concentrations based on these regular solution
% parameters and our initial concentration (assumed)

cs0vec = calcinitcs(param,cs0,ndVinit);

% If we have one constant supplied for our initial condition, let the
% script generate an initial condition
if max(size(cpcs0)) == 1
    cpcs0 = zeros(2*sl+3*el+1,1);
    cpcs0(1:sl+el) = c0;
    cpcs0(sl+el+1:2*(sl+el)) = -ndVinit;
    cpcs0(2*(sl+el)+1:2*sl+3*el) = cs0vec;
    cpcs0(end) = -ndVinit; 
end

tic     % Start timer
disp('Solving ODE...')
[t,cpcs] = ode15s(@petsolver,tr,cpcs0,options,ss,xs,ys,zp,zm,nDp,nDm, ...
                        parammat,io,currset,betamat,tp,sl,el,porosmats,tr,FFend);
toc     % End timer

tlen = max(size(t));
% PREPARE DATA FOR OUTPUT
% disp('Fomatting data for output...')
cmat = zeros(tlen,ys,xs+ss);
csmat = zeros(tlen,ys,xs);
ffvec = zeros(tlen,1);
for i=1:tlen
    cmat(i,:,:) = real(reshape(cpcs(i,1:sl+el),ys,xs+ss));
    csmat(i,:,:) = real(reshape(cpcs(i,2*(sl+el)+1:2*sl+3*el),ys,xs));
    ffvec(i) = real(sum(cpcs(i,2*(sl+el)+1:2*sl+3*el))/(xs*ys));
end
vvec = -(real(cpcs(:,end)))*(k*T/e);
vvec = vvec + Vo;

disc = struct('xs',xs,'ys',ys,'ss',ss,'sf',sf);

return;

function val = petsolver(t,cpcs,ss,xs,ys,zp,zm,nDp,nDm,parammat,io,currset, ...
                                betamat,tp,sl,el,porosmats,tr,FFend)

% This function is passed to ode15s to solve the PET equations for the
% cell. The cpcs vector is defined as the following:
%   
%   cpcs = [c_sep; c_elec; phi_sep; phi_elec; cs_elec; c0; phi0]
%   Disc.   ss*ys   xs*ys    ss*ys    xs*ys     xs*ys    2 
%
% From there, the PET equations for mass and charge conservation inside the
% cell are solved for the RHS.  All sink terms/relations are handled by the
% mass matrix.  The only explicitly needed term is the flux condition which
% yields the lithium concentration BC at the separator inlet.

%Initialize output
val = zeros(2*sl+3*el+1,1);

% Get dx and dy
dx = 1/xs;
dy = 1/ys;

% First let's get the values we need the matrices to simplify our x/y
% derivatives and separate c, phi, and cs
cvec =  cpcs(1:sl+el);
phivec = cpcs(sl+el+1:2*(sl+el));
csvec = cpcs(2*(sl+el)+1:2*sl+3*el);
phi0 = cpcs(end);

% SOLID REACTION RATE
% Extract our potential and conc vectors for the electrode
celecvec = cpcs(sl+1:sl+el);
celecmat = reshape(celecvec,ys,xs);
phielecvec = cpcs(sl+el+sl+1:2*(sl+el));
phielecmat = reshape(phielecvec,ys,xs);
iomat = reshape(io,ys,xs);
% Get our solid concentration matrix
csmat = reshape(csvec,ys,xs);
deltaphi = 0-phielecmat;
deltaphieq = calc_phieq(csmat,parammat);
eta = deltaphi-deltaphieq;
dcsdtmat = calc_R(iomat,parammat,celecmat,csmat,eta);
% Reshape and insert into output
dcsdtvec = reshape(dcsdtmat,el,1);
val(2*(sl+el)+1:2*sl+3*el) = dcsdtvec;

% MASS CONSERVATION
% First we need our flux condition into the separator
fluxperarea = sum(sum(betamat.*dcsdtmat)).*(1-tp).*dx.*dy;
% d^2c / dx^2 Calculation
cmatx = zeros(ys,ss+xs+2);
cmatx(:,2:end-1) = reshape(cvec,ys,ss+xs);
cmatx(:,1) = cmatx(:,2)+(fluxperarea*dx);       % Assumes flux equally distributed
cmatx(:,end) = cmatx(:,end-1);
% Correction for porosity
cmatxflux = (-diff(cmatx,1,2)./(dx)).*porosmats.Dx;
% Now mass conservation
d2cdx2 = -diff(cmatxflux,1,2)./(dx);
% d^2c / dy^2 Calculation
cmaty = zeros(ys+2,ss+xs);
cmaty(2:end-1,:) = reshape(cvec,ys,ss+xs);
cmaty(1,:) = cmaty(2,:);
cmaty(end,:) = cmaty(end-1,:);
cmatyflux = (-diff(cmaty,1,1)./(dy)).*porosmats.Dy;
% Need porosity correction for diffusivity in y flux
d2cdy2 = -diff(cmatyflux,1,1)./(dy);
% Sum the two terms, shift to vector, insert into output
d2csum = d2cdx2+d2cdy2;
val(1:sl+el) = reshape(d2csum,sl+el,1);

% CHARGE CONSERVATION
% Average concentration matrix in x direction
cmatavgx = (cmatx(:,1:end-1)+cmatx(:,2:end))/2;
% Average concentration matrix in y direction
cmatavgy = (cmaty(1:end-1,:)+cmaty(2:end,:))/2;
% Potential field for x derivatives
phimatx = zeros(ys,ss+xs+2);
phimatx(:,2:end-1) = reshape(phivec,ys,ss+xs);
phimatx(:,1) = phi0;
phimatx(:,end) = phimatx(:,end-1);
% Potential field for y derivatives
phimaty = zeros(ys+2,ss+xs);
phimaty(2:end-1,:) = reshape(phivec,ys,ss+xs);
phimaty(1,:) = phimaty(2,:);
phimaty(end,:) = phimaty(end-1,:);
% Calculate current density for x and y
currdensx = porosmats.Dx.*((-(nDp-nDm).*(diff(cmatx,1,2))- ...
                (zp*nDp+zm*nDm).*cmatavgx.*diff(phimatx,1,2))./dx);
currdensy = porosmats.Dy.*((-(nDp-nDm).*(diff(cmaty,1,1))- ...
                (zp*nDp+zm*nDm).*cmatavgy.*diff(phimaty,1,1))./dy);
% Finally, we take the divergence of our current densities
divcdx = -diff(currdensx,1,2)./dx;
divcdy = -diff(currdensy,1,1)./dy;
% Now we have our two contributions, add them and insert into output
divcdsum = divcdx + divcdy;
val(sl+el+1:2*(sl+el)) = reshape(divcdsum,sl+el,1);

% Finally we need our boundary conditions
% CHARGE CONSERVATION BC
val(end) = currset;

return;

function R = calc_R(io,parammat,cmat,csmat,eta)

% This function calculates the reaction rate of the particles based on the
% overpotential.  A transfer coefficient of 0.5 is assumed
ecd = -2.*io.*sqrt(cmat.*csmat.*(1./csmat)).*exp((parammat/2).*(1-2.*csmat));
R = ecd.*sinh(eta);

return;

function phieq = calc_phieq(cs,avec)

% Calculates equilibrium potential
phieq = -log(cs./(1-cs))-avec.*(1-2.*cs);

return

function M = genMassMatrix(ss,xs,ys,tp,betamat,poros)

% Need to formulate a sparse mass matrix
% Get some numbers to make things easier
sl = ss*ys;
el = xs*ys;

% Initialize the matrix
M = sparse(2*sl+3*el+1,2*sl+3*el+1);

% Insert our values into the mass matrix
M(1:sl,1:sl) = speye(sl,sl);
M(sl+1:sl+el,sl+1:sl+el) = speye(el,el)*poros;
M(sl+1:sl+el,2*(sl+el)+1:2*sl+3*el) = ...
    sparse(1:el,1:el,reshape(betamat,el,1),el,el)*(1-tp);
M(2*sl+el+1:2*(sl+el),2*(sl+el)+1:2*sl+3*el) = ...
    sparse(1:el,1:el,reshape(betamat,el,1),el,el);
M(2*(sl+el)+1:2*sl+3*el,2*(sl+el)+1:2*sl+3*el) = speye(el,el);
M(end,2*(sl+el)+1:2*sl+3*el) = 1/(xs*ys);

return;

function param = size2regsoln(sz)

% This function returns the regular solution parameter for a given particle
% size.  The particles are assumed to be spherical, from which an
% Area:Volume ratio is calculated and then used to produce a regular
% solution model parameter.

% Parameters from polynomial curve fit
p1 = -54.17;
p2 = -1.532;
p3 = -5.68;
p4 = 3.459;

% f(x) = p1*x^3 + p2*x^2 + p3*x + p4
 
AV = 3.475./sz;     % A:V = (4 pi r^2 / (4/3) pi r ^3)
param = p1.*AV.^3 + p2.*AV.^2 + p3.*AV + p4;

len = max(size(param));
for i=1:len
    if param(i) < 2
        param(i) = 2;
    end
end

return;

function val = calcinitcs(paramvec,cs0,ndVinit)

% This function returns a vector of particle concentrations for various
% particles with different regular solution parameters based on an assumed
% mean particle radius with initial concentration cs0

% This script facilitates the initial condition of the DAE solver
% Assume 4.5kT for large particles

len = max(size(paramvec));
csinit = cs0*ones(len,1);

opt = optimset('Display','off');
val = fsolve(@calcerr,csinit,opt,paramvec,ndVinit);


return;

function mu = calcmu(cs,params)

mu = -log(cs./(1-cs))-params.*(1-2.*cs);

return;

function err = calcerr(csvec,paramvec,ndVinit)

err =   ndVinit - calcmu(csvec,paramvec);

return;

function [value, isterminal, direction] = events(t,cpcs,ss,xs,ys,zp,zm,nDp,nDm,parammat,io,currset, ...
                                betamat,tp,sl,el,porosmats,tr,FFend)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
perc = (t/tfinal) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)                       

ff = real(sum(cpcs(2*(sl+el)+1:2*sl+3*el))/(xs*ys));
value = ff - FFend;
isterminal = 1;
direction = 0;
                        
return;