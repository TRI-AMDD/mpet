function [t,cpcs,clv,plv,ffvec,vvec,disc] = pm_1Dcc_electcond_rev1(io,currset,mcond,Rcont)

% This script simulates the pseudocapacitor discharge for phase separating
% particles.  This script includes proper handling of the current density
% as well as proper handling of the porosity.

% CONSTANT CURRENT DISCHARGE

% Electrolyte concentration  		1 : ss + steps
% Electrolyte potential         	ss + steps + 1 : 2*(ss + steps)
% Solid conc.                       2*(ss + steps) + 1 : 2*ss + 3*steps
% Metal potential                   2*ss+3*steps+1 : 2*ss+4*steps
% Potential boundary condition      end

tic

% Discharge settings
% io = .11;
% currset = .001;
L = 1;

% Discretization setup
sf = 0.5;
steps = 20;
ss = sf * steps;
tfinal = L/abs(currset);
tsteps = 400;
tr = linspace(0,tfinal,tsteps);
dt = tr(2)-tr(1);

% Constants
k = 1.381e-23;
T = 298;
e = 1.602e-19;

% Material properties
a = 4.5;
Lp = 0.5;       % Loading percent by volume
csmax = 22.8;   % mol/L
poros = .3;
beta = (1-poros)*Lp*csmax;
% mcond = .03;
V0 = 3.422;

% Electrolyte properties
zp = 1;
zm = 1;
Dp = 1.25e-10;      % m^2/s
Dm = 4e-10;         % m^2/s
Damb = (Dp*Dm*(zp+zm))/(zm*Dm+zp*Dp);
nDp = Dp/Damb;
nDm = Dm/Damb;
td = L^2/Damb;
tp = (zp*Dp)/(zp*Dp+zm*Dm);

% Initialize vectors
if currset > 0
    cs0 = 0.01;
else
    cs0 = 0.90;
end
cpcs0 = zeros(2*ss+4*steps+1,1);            % Empty vector
cpcs0(1:ss+steps) = ones(ss+steps,1);       % Elect. conc.
cpcs0(2*(ss+steps)+1:2*ss+3*steps) = cs0*(ones(steps,1)+.1*rand(steps,1));   % Solid conc.
phip0 = calcmu(cs0,a);            % Get equilibrium potential
cpcs0(ss+steps+1:2*(ss+steps)) = phip0*ones(ss+steps,1);    % Elect. pot.
cpcs0(end) = -phip0(1);     
M = genMass(ss,steps,beta,tp,L);

options = odeset('MASS',M,'MassSingular','yes','MStateDependence','none','RelTol',1e-7,'AbsTol',1e-8);
%options = odeset('MASS',M,'MassSingular','yes','MStateDependence','none');

[t,cpcs] = ode15s(@calcrhs,tr,cpcs0,options,ss,steps,nDp,nDm, ...
                                poros,zp,zm,currset,a,io,tp,beta,L,mcond);
                 
cpcs = real(cpcs);              % Only real values
                            
% generate our filling fraction and voltage outputs
tlen = max(size(t));                
ffvec = zeros(tlen,1);
vvec = V0-(k*T/e).*cpcs(:,end)' - currset*Rcont;
for i=1:tlen
    ffvec(i) = sum(cpcs(i,2*(ss+steps)+1:2*ss+3*steps))/(steps);
end

currvec = zeros(tlen,1);

clv = linspace(-sf,1,ss+steps);
plv = linspace(0,1,steps);

disc = struct('ss',ss,'steps',steps,'L',L);

toc

return;

function val = calcrhs(t,cpcs,ss,steps,nDp,nDm,poros,zp,zm,currset,a,io,tp,beta,L,mcond)

val = zeros(2*ss+4*steps+1,1);
dx = L/steps;

% Need a porosity vetor
porosvec = ones(ss+steps+1,1);
porosvec(ss+1) = (poros+1)/2;
porosvec(ss+2:end) = poros;

% Apply our Bruggeman relation
porosvec = porosvec.^(3/2);

% Calculate our concentration outside the cell.  The last value in the cpcs 
% vector is concentration
tmpc = zeros(ss+steps+2,1);
tmpc(2:end-1) = cpcs(1:ss+steps);
tmpc(1) = dx*(1-tp)*beta*currset + tmpc(2);
tmpc(end) = tmpc(end-1);
cflux = porosvec.*(-diff(tmpc)./dx);
val(1:ss+steps) = -diff(cflux)./dx;

% Then we take the divergence of our current density
tmpp = zeros(ss+steps+2,1);
tmpp(2:end-1) = cpcs(ss+steps+1:2*(ss+steps));
tmpp(1) = cpcs(end); tmpp(end) = tmpp(end-1);
cavg = (tmpc(1:end-1)+tmpc(2:end))/2;
currdens = porosvec .* (-(nDp-nDm).*(diff(tmpc)./dx)- ...
                (zp*nDp+zm*nDm).*cavg.*(diff(tmpp)/dx));
val(ss+steps+1:2*(ss+steps)) = -diff(currdens)./dx;

% Extract our phi_metal potentials
phim = cpcs(2*ss+3*steps+1:end-1);
% Construct our charge conservation the electron conducting phase
phimtmp = zeros(2+steps,1);
phimtmp(2:end-1) = phim;
phimtmp(1) = phimtmp(2);        % No flux condition
% Robin condition - set ground and flux out
phim(end-1) = beta*currset*mcond*dx;
phimflux = -mcond.*diff(phimtmp)/dx;
diffphim = -diff(phimflux)./dx;
val(2*ss+3*steps+1:end-1) = diffphim;

% Extract our solid concentrations
cstmp = cpcs(2*(ss+steps)+1:2*ss+3*steps);

% Next we calculate our exchange current and overpotential values
io1 = calcexcurr(cstmp,cpcs(ss+1:ss+steps),a,io);
eta1 = (phim-cpcs(ss+steps+ss+1:2*(ss+steps)))-(-calcmu(cstmp,a));

% Insert into our vector
val(2*(ss+steps)+1:2*ss+3*steps) = io1.*sinh(eta1/2);

% Finally our total charge conservation
val(end) = currset;

return;

function val = calcexcurr(cs,celec,a,io)

val = -2*io*sqrt(celec.*cs.*(1-cs)).*exp((a/2).*(1-2.*cs));

return;

function val = calcmu(cs,a)

% Calculates the chemical potential for one of the layers

val = log(cs./(1-cs)) + a*(1-2*cs);

return;

function M = genMass(ss,steps,beta,tp,L)

% Function generates the Mass matrix for our DAE
M = sparse(2*ss+4*steps+1,2*ss+4*steps+1);

% Separator electrolyte
M(1:ss,1:ss) = speye(ss,ss);

% Electrode electrolyte conc.
M(ss+1:ss+steps,ss+1:ss+steps) = speye(steps,steps);
M(ss+1:ss+steps,2*(ss+steps)+1:2*ss+3*steps) = beta*(1-tp)*speye(steps,steps);

% Potential
M(2*ss+steps+1:2*(ss+steps),2*(ss+steps)+1:2*ss+3*steps) = beta*speye(steps,steps);

% Solid concentration
M(2*ss+2*steps+1:2*ss+3*steps,2*ss+2*steps+1:2*ss+3*steps) = speye(steps,steps);

% Metal potential
M(2*ss+3*steps+1:2*ss+4*steps,2*ss+2*steps+1:2*ss+3*steps) = -beta*speye(steps,steps);

% Current conservation
M(end,2*ss+2*steps+1:2*ss+3*steps) = L/(steps);

return;
