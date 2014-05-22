function [t,cpcs1e,cmat,csmat,ffvec1e,vvec1e,disc] = runsim_RBS()

% This script simulates the Gaberscek experiment as closely as possible.
% The order of steps is as follows:

%  1. Start with electrode empty (SOC = .01), charge at C/10 to SOC = .6
%  2. Equilibrate cell
%  3. Discharge at C-rate to 40% SOC
%  4. Equilibrate cell
%  5. Charge at C-rate to 60% SOC

% We will perform this with multiple discharge rates then save the figures
% for charge and discharge (on top of each other) - from these figures the
% voltage gap data will be read (manually) and plotted. 

% Ideally, the contact resistance needs to be fit to the data, but it
% should be acceptable to use a similar value to was initially fit (scaled
% by the new current density for C/10).

% Set our discretization
xs = 80;
ys = 1;
tsteps = 200;

crate = 5.0;

% Construct our particle size distribution
% szs = abs(normrnd(mean,stddev,xs*ys,1));  

% LOG NORMAL
mean = 160; % these need to be in nm for the size2regsln function
stddev = 20;
var = stddev.^2;
%var = 10;
m = mean;
v = var;	
mu = log((m^2)/sqrt(v+m^2));
sig = sqrt(log(v/(m^2)+1));
szs = lognrnd(mu,sig,xs*ys,1);

% io = .364;
% onec = .00364;
Rcont = 0;
FFend = 0.95;

% Calculate our exchange current density and 1C dimensionless currents
cp.L = 50e-6;                       % electrode length, m
cp.poros = 0.4;                     % porosity, dimensionless
cp.Ds = 50e-9;                      % particle length, m
cp.Lp = 0.69;                        % Loading percent
cp.c0 = 1;                          % Initial electrolyte concentration, mol/L

% UPDATED 9/18/13 - CORRECTED THICKNESS, cp.ap ACTIVE AREA CALCULATION, io
cp.thickness = 20e-9*ones(xs*ys,1);     % Particle thickness
cp.ap = 2./(cp.thickness);
cp.csmax = 22.6;                    % max LiFePO4 lithium, mol/L
cp.Vo = 3.422;                      % std. potential, V
cp.F = 96485;                       % Faraday's constant
% cp.io = 0.741;                      % exchange current density, A/m^2  (guess)
%cp.io = 0.00741;                   % Exchange current density, A/m^2
cp.io = 0.16;                   % Exchange current density, A/m^2

% Order of these values inserted depends on what's provided
% This assumes tp and Damb provided
%cp.Damb = 1.48e-10;                 % Ambipolar diffusivity, m^2/s
%cp.tp = 0.35;                       % Transference number
cp.zp = 1;                          % Cation charge number
cp.zm = 1;                          % Anion charge number
%cp.Dm = cp.Damb/(2*cp.tp);          % Anion diffusivity
%cp.Dp = (cp.tp*cp.Dm)/(1-cp.tp);    % Cation diffusivity
cp.Dp = 2.2e-10;
cp.Dm = 2.94e-10;
cp.Damb = ((cp.zp+cp.zm)*cp.Dp*cp.Dm)/(cp.zp*cp.Dp + cp.zm*cp.Dm);
cp.tp = cp.zp*cp.Dp / (cp.zp*cp.Dp + cp.zm*cp.Dm);

% Dimensionless exchange current density
io = (cp.L^2/cp.Damb) .* cp.ap .* (cp.io./(cp.csmax*1000*cp.F));
onec = (cp.L^2/cp.Damb) / 3600;

%crates = [1/50, 1/100, 1/131, 1/200];
% Need other data sets (2013-01-03)
%crates = [1/1000, 1/200, 1/131, 1/100, 1/50, 1/30];
% crates = [1/1000];
%
%numc = max(size(crates));
%
%% Initialize our electrode
%[t,cpcs1,cmat,csmat,ffvec1,vvec1,disc] = ...
%                pm_2dvarsizecc_rev1(io,onec/10,szs,0,.2,Rcont,cp,xs,ys,tsteps);
%% Equilibrate the particles
%[t,cpcs1e,cmat,csmat,ffvec1e,vvec1e,disc] = ...
%                pm_2dvarsizecc_rev1(io,0,szs,cpcs1(end,:)',.5,Rcont,cp,xs,ys,tsteps);
%
%for i=1:numc
%    
%% io,currset,param,cpcs0,FFend,Rcont           
%[t,cpcs2,cmat,csmat,ffvec2,vvec2,disc] = ...
%                pm_2dvarsizecc_rev1(io,onec*crates(i),szs,cpcs1e(end,:)',.7,Rcont,cp,xs,ys,tsteps);            
%[t,cpcs2e,cmat,csmat,ffvec2e,vvec2e,disc] = ...
%                pm_2dvarsizecc_rev1(io,0,szs,cpcs2(end,:)',.5,Rcont,cp,xs,ys,tsteps);                        
%[t,cpcs3,cmat,csmat,ffvec3,vvec3,disc] = ...
%                pm_2dvarsizecc_rev1(io,-onec*crates(i),szs,cpcs2e(end,:)',.2,Rcont,cp,xs,ys,tsteps);
%[t,cpcs1e,cmat,csmat,ffvec1e,vvec1e,disc] = ...
%                pm_2dvarsizecc_rev1(io,0,szs,cpcs1(end,:)',.5,Rcont,cp,xs,ys,tsteps);
%
%% Construct the output file name
%str = ['varsz_' num2str(mean) '_' num2str(var) '_' num2str(crates(i)) '.mat'];
%save(str);
%
%end 

[t,cpcs1e,cmat,csmat,ffvec1e,vvec1e,disc] = ...
                pm_2dvarsizecc(io,onec*crate,szs,0,FFend,Rcont,cp,xs,ys,tsteps);
          
            
return;
