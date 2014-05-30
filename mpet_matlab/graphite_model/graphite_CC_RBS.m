function [t,cpcs,clv,plva,plvc,ffvec,vvec,disc] = graphite_CC_RBS(cporos)

% This script solves a full cell with an anode/separator/cathode using the 
% pseudocapacitor approximation.  The anode and cathode particles each 
% have their own OCP functions and separate exchange current density values
% 
% Materials assumed here are:
%     Anode: Li (m)
%     Cathode: LiC6
%     Electrolyte: LiPF6 1M in EC/DMC

% DISCRETIZATION
% Electrolyte concentration  		1 : len
% Electrolyte potential         	len + 1 : 2 * len 
% Anode solid conc.                 2*len+1 : 2*len + asteps
% Cathode layer 1 solid conc.       2*len+asteps+1 : 2*len+asteps+csteps
% Cathode layer 2 solid conc.       2*len+asteps+csteps+1 : 2*len+asteps+2*csteps
% Phi_(m,a)                     	end

tic

% Discharge settings
ioa = 1;
ioc = 1;
currset = .01;
phima = 4.63;

% Discretization setup
catlen = 10;
seplen = 1;
anlen = 10;
csteps = 80;
asteps = (anlen/catlen).*csteps;
ss = ceil((seplen/catlen).*csteps);
% asteps = 50;
% csteps = 50;
% ss = 20;
len = asteps + ss + csteps;
% Lengths
dxa = catlen / csteps;
dxs = seplen / ss;
dxc = anlen / asteps;
% asteps = 50;
% csteps = 50;
% ss = 20;
% len = asteps + ss + csteps;
% % Lengths
% catlen = 20;
% seplen = 1;
% anlen = 12;
% dxa = catlen / csteps;
% dxs = seplen / ss;
% dxc = anlen / asteps;
% Discretization vector to handle different mesh/lengths
disc = zeros(len+2,1);
disctmp = zeros(len,1);
disctmp(1:asteps) = dxa;
disctmp(asteps+1:asteps+ss) = dxs;
disctmp(asteps+ss+1:len) = dxc;
disc(2:end-1) = disctmp;
disc(1) = disc(2); disc(end) = disc(end-1);
% Need these vectors for our flux and 
fdisc = (disc(1:end-1)+disc(2:end))/2;
divdisc = (fdisc(1:end-1) + fdisc(2:end))/2;

% Time discretization
tfinal = 10000;
tsteps = 1000;
tr = linspace(0,tfinal,tsteps);
dt = tr(2)-tr(1);

% Constants
k = 1.381e-23;
T = 298;
e = 1.602e-19;

% Material properties
cc = 30;
ac = 3.4;
bc = 1.4;
V0 = 5.3174;

% [ac, bc, V0] = calcabc(cc);
% ac = 3.3955;             % Graphite regular solution parameter
% bc = 0.7065;             % Graphite interaction parameter
% cc = 20;

aa = 1; % NOT USED
aporos = 1;                 % Assume volume averaged anode porosity - high to limit anode transport effects
% cporos = 0.3;
Lpa = 1;
Lpc = .8;
abeta0 = 76.9;              % cs_max (anode)
cbeta0 = 28.35/2;           % cs_max (cathode)
abeta = Lpa*abeta0; % Inconsistent, anode porosity = 0 and 1
cbeta = Lpc*(1-cporos)*cbeta0;
% V0 = 4.7071;     % Standard potential (offset, reference to Li metal)

% Electrolyte properties
cp.Damb = 1.48e-10;                 % Ambipolar diffusivity, m^2/s
cp.tp = 0.35;                       % Transference number
cp.zp = 1;                          % Cation charge number
cp.zm = 1;                          % Anion charge number
cp.Dm = cp.Damb/(2*cp.tp);          % Anion diffusivity
cp.Dp = (cp.tp*cp.Dm)/(1-cp.tp);    % Cation diffusivity

zp = cp.zp;
zm = cp.zm;
Dp = cp.Dp;      
Dm = cp.Dm;         

% nDp = 0.5*(Dp+Dm)/Dm;
% nDm = 0.5*(Dp+Dm)/Dp;
nDp = 1/(2*(1-cp.tp));
nDm = 1/(2*cp.tp);
tp = cp.tp;

% Initialize vectors
cs0 = 0.01;
phieqc0 = calc_phieqc(cs0,cs0,ac,bc,cc,1);
cpcs0 = zeros(2*len + asteps + 2*csteps + 1,1);       % Empty vector
cpcs0(1:len) = ones(len,1);
cpcs0(len+1:2*len) = -phieqc0.*ones(len,1);
cpcs0(2*len+1:2*len+asteps) = ones(asteps,1);
cpcs0(2*len+asteps+1:2*len+asteps+2*csteps) = cs0.*ones(2*csteps,1);

% Initialize a porosity vector
% porvec = ones(len+2,1);
% porvec(1:asteps) = aporos;
% porvec(asteps+ss+1:len) = cporos;
% porvec(1) = porvec(2); porvec(end) = porvec(end-1);
% porvec = (porvec(1:end-1) + porvec(2:end))/2;
% porvec = porvec.^(3/2);

% POROSITY VECTOR WITH INDICES MOVED TO 1 INSIDE SEPARATOR (RIGHT OR WRONG?)
% porvec = ones(len+1,1);
% porvec(1:asteps) = aporos;
% porvec(end-csteps:end) = cporos;
% porvec = porvec.^(3/2);     % Bruggeman

% ORIGINAL POROSITY VECTOR (SEEMS OK IN VIDEO)
porvec = ones(len+1,1);
porvec(1:asteps-1) = aporos;
porvec(end-csteps+1:end) = cporos;
porvec = porvec.^(3/2);     % Bruggeman

M = genMass(len,asteps,ss,csteps,abeta,cbeta,tp,anlen,catlen,aporos,cporos,dxc);
%spy(M)

%options = odeset('MASS',M,'MassSingular','yes','MStateDependence','none','RelTol',1e-7);
% options = odeset('MASS',M,'MassSingular','yes','MStateDependence','none');
options = odeset('MASS',M,'MassSingular','yes','MStateDependence','none','Events',@events);

[t,cpcs] = ode15s(@calcrhs,tr,cpcs0,options,len,ss,asteps,csteps,porvec,fdisc, ...
                        divdisc,nDp,nDm,zp,zm,ioa,ioc,ac,bc,cc,aa,phima,tr,currset);

% generate our filling fraction and voltage outputs
tlen = max(size(t));                
ffvec = zeros(tlen,1);
vvec = (k*T/e).*(V0-cpcs(:,end)');
%vvec = zeros(tlen,1);
for i=1:tlen
    ffvec(i) = sum(cpcs(i,2*len+asteps+1:2*len+asteps+2*csteps))/(2*csteps);
end

%t = 0;
%cpcs = 0;
%ffvec = 0;
%vvec = 0;
clv = [linspace(-anlen,0,asteps),linspace(0,seplen,ss),linspace(seplen,seplen+catlen,csteps)];
plva = linspace(-anlen,0,asteps);
plvc = linspace(seplen,seplen+catlen,csteps);
disc = struct('ss',ss,'asteps',asteps,'csteps',csteps, ...
                'catlen',catlen,'seplen',seplen,'anlen',anlen,...
                'cbeta',cbeta0,'abeta',abeta0,'aporos',aporos,'cporos',cporos);
toc

return;

function val = calcrhs(t,cpcs,len,ss,asteps,csteps,porvec,fdisc,divdisc,nDp,nDm,zp, ...
                            zm,ioa,ioc,ac,bc,cc,aa,phima,tr,currset)
                  
val = zeros(2*len+asteps+2*csteps+1,1);

% First we handle the concentration                  
ctmp = zeros(len+2,1);
ctmp(2:end-1) = cpcs(1:len);
ctmp(1) = ctmp(2); ctmp(end) = ctmp(end-1);         % No flux
cflux = -porvec.*(diff(ctmp)./fdisc);
divc = -diff(cflux)./divdisc;
val(1:len) = divc;

% Next the current density
tmpp = zeros(len+2,1);
tmpp(2:end-1) = cpcs(len+1:2*len);
tmpp(1) = tmpp(2); tmpp(end) = tmpp(end-1);         % Insulation
cavg = (ctmp(1:end-1)+ctmp(2:end))./2;
currdens = porvec.*( -(nDp-nDm).*diff(ctmp)./fdisc - ...
                (zp*nDp+zm*nDm).*cavg.*diff(tmpp)./fdisc);
val(len+1:2*len) = -diff(currdens)./divdisc;

% Now the reaction rates
acstmp = cpcs(2*len+1:2*len+asteps);
ccs1tmp = cpcs(2*len+asteps+1:2*len+asteps+csteps);
ccs2tmp = cpcs(2*len+asteps+csteps+1:2*len+asteps+2*csteps);
% Then our exchange current densities
ioan = calcexcurra(cpcs(1:asteps),acstmp,ioa);
iocat1 = calcexcurrc(ccs1tmp,ccs2tmp,cpcs(asteps+ss+1:len),ioc,ac,bc,cc);
iocat2 = calcexcurrc(ccs2tmp,ccs1tmp,cpcs(asteps+ss+1:len),ioc,ac,bc,cc);

% Get our overpotentials
phima = cpcs(end);
etaa = (phima - cpcs(len+1:len+asteps)) - calc_phieqa(acstmp,aa);
etac1 = (0 - cpcs(len+asteps+ss+1:2*len)) - calc_phieqc(ccs1tmp,ccs2tmp,ac,bc,cc,cpcs(asteps+ss+1:len));
etac2 = (0 - cpcs(len+asteps+ss+1:2*len)) - calc_phieqc(ccs2tmp,ccs1tmp,ac,bc,cc,cpcs(asteps+ss+1:len));

% Calculate the reaction rates
val(2*len+1:2*len+asteps) = ioan.*sinh(etaa/2);
val(2*len+asteps+1:2*len+asteps+csteps) = iocat1.*sinh(etac1/2);
val(2*len+asteps+csteps+1:2*len+asteps+2*csteps) = iocat2.*sinh(etac2/2);

val(end) = currset;
                
return;

function excurra = calcexcurra(celec,csa,ioa)

% Calculates the exchange current density of the anode material
excurra = -2.*ioa.*sqrt(celec);

return;

function excurrc = calcexcurrc(csc1,csc2,celec,ioc,ac,bc,cc)

% Calculates the exchange current density of the cathode material
excurrc = -2.*ioc.*sqrt(celec.*csc1.*(1-csc1)).*exp((ac/2).* ...
            (1-2.*csc1)+(bc/2).*csc2+(cc/2).*(csc2.*(1-csc2).*(1-2.*csc1)));

return;

function phieqa = calc_phieqa(cs,aa)

% Returns equilibrium potential for anode material
% For lithium metal use constant value
sz = max(size(cs));
phieqa = zeros(sz,1);

return;

function phieqc = calc_phieqc(cs1,cs2,ac,bc,cc,celec)

% Returns equilibrium potential for cathode material
% Use regular solution model for potential
phieqc = -log(cs1./(1-cs1)) - ac.*(1-2.*cs1) - bc.*cs2 - ...
            cc.*cs2.*(1-cs2).*(1-2.*cs1) + log(celec);

return;

function M = genMass(len,asteps,ss,csteps,abeta,cbeta,tp,anlen,catlen,aporos,cporos,dxc)

% Returns mass matrix 
% Initialize mass matrix
M = sparse(2*len+asteps+2*csteps+1,2*len+asteps+2*csteps+1);

% Electrolyte concentration
%M(1:len,1:len) = speye(len,len);
M(1:asteps,1:asteps) = aporos.*speye(asteps);
M(asteps+1:asteps+ss,asteps+1:asteps+ss) = speye(ss);
M(asteps+ss+1:len,asteps+ss+1:len) = cporos.*speye(csteps);

M(1:asteps,2*len+1:2*len+asteps) = abeta.*(1-tp).*speye(asteps);
M(asteps+ss+1:len,2*len+asteps+1:2*len+asteps+csteps) = cbeta.*(1-tp).*speye(csteps);
M(asteps+ss+1:len,2*len+asteps+csteps+1:2*len+asteps+2*csteps) = cbeta.*(1-tp).*speye(csteps);

% Electrolyte potential
M(len+1:len+asteps,2*len+1:2*len+asteps) = abeta.*speye(asteps);
M(len+asteps+ss+1:2*len,2*len+asteps+1:2*len+asteps+csteps) = cbeta.*speye(csteps);
M(len+asteps+ss+1:2*len,2*len+asteps+csteps+1:2*len+asteps+2*csteps) = cbeta.*speye(csteps);

% Solids
M(2*len+1:2*len+asteps+csteps,2*len+1:2*len+asteps+csteps) = speye(asteps+csteps);
M(2*len+asteps+csteps+1:2*len+asteps+2*csteps,2*len+asteps+csteps+1:2*len+asteps+2*csteps) = speye(csteps,csteps);

% Current definition
M(end, 2*len+asteps+1:2*len+asteps+csteps) = dxc;
M(end, 2*len+asteps+csteps+1:2*len+asteps+2*csteps) = dxc;

% Charge balance
% M(end,2*len+1:2*len+asteps) = anlen/asteps;
% M(end,2*len+asteps+1:2*len+asteps+2*csteps) = catlen/csteps;

return;

function [value, isterminal, direction] = events(t,cpcs,len,ss,asteps,csteps,porvec,fdisc,divdisc,nDp,nDm,zp, ...
                            zm,ioa,ioc,ac,bc,cc,aa,phima,tr,currset)
                        
value = 0;
isterminal = 0;
direction = 0;
tfinal = tr(end);
tsteps = max(size(tr));
perc = ((t/tsteps) / (tfinal/tsteps)) * 100;
dvec = [num2str(perc),' percent completed'];
disp(dvec)                       
                        
return;
