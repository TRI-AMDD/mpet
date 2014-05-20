function [t,cpcs,cmat,csmat,ffvec,vvec,disc] = ...
                pm_2dvarsizecc_rev1(io,currset,szs,cpcs0,FFend,Rcont,cp,xs,ys,tsteps)

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

% TO DO:
%   Add porosity effects in conductivity and diffusivity
%   Variable diffusivities (CST)  (rev2?)

% CONSTANTS
k = 1.381e-23;
e = 1.602e-19;
T = 298;

% CELL DIMENSIONS/LOADING
poros = cp.poros;       % Porosity
c0 = cp.c0;             % mol/L - initial concentration of electrolyte
Lp = cp.Lp;             % Loading percent (vol. fraction active material)

% ELECTROLYTE PROPERTIES
zp = cp.zp;
zm = cp.zm;
Dp = cp.Dp;
Dm = cp.Dm;
Damb = cp.Damb;
nDp = Dp/Damb;
nDm = Dm/Damb;
tp = cp.tp;

% CALCULATE NON-DIMENSIONAL CURRENT FROM C-RATE
% td = L^2/Damb;
% onec = td/3600;
% currset = onec*crate;

% MATERIAL PROPERTIES 
a0 = 4.51;              % reference regular soln parameter - ONLY FOR INITIAL VALUE CALCULATION
csmax = cp.csmax;           % mol/L - maximum solid concentration
Vo = cp.Vo;

% MESH DIMENSIONS
sf = 0.5;      % Length of separator as fraction of cathode
% xs = 20;
% ys = 10;
ss = sf*xs;
tinit = 0;
if currset ~= 0
    tfinal = 1/abs(currset);
else
    tfinal = 1000;
end
tr = linspace(tinit,tfinal,tsteps);
sl = ss*ys;
el = xs*ys;

% NON-DIMENSIONAL QUANTITIES
beta = ((Lp*(1-poros)))*(csmax/c0);
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
Dxmat(:,ss+1) = (1+poros^1.5)/2;
Dxmat(:,ss+2:end) = poros^1.5;
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
c0 = 1;
if (currset>0)
    cs0 = .01;      % Start empty for discharge
else
    cs0 = .99;      % Start full for charging
end


% Use our passed size vector
param = size2regsoln(szs);
parammat = reshape(param,ys,xs);
iomat = reshape(io,ys,xs);
% Now we generate initial concentrations based on these regular solution
% parameters and our initial concentration (assumed)
cs0vec = calcinitcs(param,cs0,a0);

% If we have one constant supplied for our initial condition, let the
% script generate an initial condition
if max(size(cpcs0)) == 1
    cpcs0 = zeros(2*sl+3*el+1,1);
    cpcs0(1:sl+el) = c0;
    cpcs0(sl+el+1:2*(sl+el)) = -calc_phieq(cs0,a0);
    cpcs0(2*(sl+el)+1:2*sl+3*el) = cs0vec;
    cpcs0(end) = -calc_phieq(cs0,a0);
end

% tic     % Start timer
% disp('Solving ODE...')
[t,cpcs] = ode15s(@petsolver,tr,cpcs0,options,ss,xs,ys,zp,zm,nDp,nDm, ...
                        parammat,iomat,currset,betamat,tp,sl,el,porosmats,tr,FFend);
% toc     % End timer

tlen = max(size(t));
% PREPARE DATA FOR OUTPUT
% disp('Fomatting data for output...')
cmat = zeros(tlen,ys,xs+ss);
csmat = zeros(tlen,ys,xs);
ffvec = zeros(tlen,1);
vvec = zeros(tlen,1);
for i=1:tlen
    cmat(i,:,:) = real(reshape(cpcs(i,1:sl+el),ys,xs+ss));
    csmat(i,:,:) = real(reshape(cpcs(i,2*(sl+el)+1:2*sl+3*el),ys,xs));
    ffvec(i) = real(sum(cpcs(i,2*(sl+el)+1:2*sl+3*el))/(xs*ys));
end
vvec = -(real(cpcs(:,end)))*(k*T/e);
vvec = vvec + Vo - currset * Rcont;
% vvec = vvec + Vo;

disc = struct('xs',xs,'ys',ys,'ss',ss,'sf',sf);
bm = betamat;

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
c0 = cpcs(end-1);
phi0 = cpcs(end);

% SOLID REACTION RATE
% Extract our potential and conc vectors for the electrode
celecvec = cpcs(sl+1:sl+el);
celecmat = reshape(celecvec,ys,xs);
phielecvec = cpcs(sl+el+sl+1:2*(sl+el));
phielecmat = reshape(phielecvec,ys,xs);
% Get our solid concentration matrix
csmat = reshape(csvec,ys,xs);
deltaphi = 0-phielecmat;
deltaphieq = calc_phieq(csmat,parammat);
eta = deltaphi-deltaphieq;
dcsdtmat = calc_R(io,parammat,celecmat,csmat,eta);
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

% ERROR FROM PREVIOUS MODEL, CORRECTED 20130524
% ecd = -2.*io.*sqrt(cmat.*csmat.*(1./csmat)).*exp((parammat/2).*(1-2.*csmat));

ecd = -2.*io.*sqrt(cmat.*csmat.*(1-csmat)).*exp((parammat/2).*(1-2.*csmat));
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
% size (C3 particle, measured in nm in the 100 direction).
% The particles are assumed to be spherical, from which an
% Area:Volume ratio is calculated and then used to produce a regular
% solution model parameter.

% Parameters from polynomial curve fit
p1 = -1.168e4;
p2 = 2985;
p3 = -208.3;
p4 = -8.491;
p5 = -10.25;
p6 = 4.516;

% f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5^x + p6
 
AV = 3.6338./sz;     % From Dan's paper with particle sizes
param = p1.*AV.^5 + p2.*AV.^4 + p3.*AV.^3 + p4.*AV.^2 + p5.*AV + p6;

len = max(size(param));
for i=1:len
    if param(i) < 2
        param(i) = 2;
    end
end

return;

function val = calcinitcs(paramvec,cs0,a0)

% This function returns a vector of particle concentrations for various
% particles with different regular solution parameters based on an assumed
% mean particle radius with initial concentration cs0

% This script facilitates the initial condition of the DAE solver
% Assume 4.5kT for large particles
phi0 = calc_phieq(cs0,a0);
len = max(size(paramvec));
csinit = cs0*ones(len,1);

opt = optimset('Display','off');
val = fsolve(@calcerr,csinit,opt,paramvec,phi0);


return;

function mu = calcmu(cs,params)

mu = -log(cs./(1-cs))-params.*(1-2.*cs);

return;

function err = calcerr(csvec,paramvec,phi0)

err = calcmu(csvec,paramvec) - phi0;

return;

function [value, isterminal, direction] = events(t,cpcs,ss,xs,ys,zp,zm,nDp,nDm,parammat,io,currset, ...
                                betamat,tp,sl,el,porosmats,tr,FFend)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)                       

ff = real(sum(cpcs(2*(sl+el)+1:2*sl+3*el))/(xs*ys));
value = ff - FFend;
isterminal = 1;
direction = 0;
                        
return;
