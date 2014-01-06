function trackint4(filename)

% This function tracks the red-gold interface in the graphite electrode as
% it is discharged at constant potential.  The function returns the
% interface position versus position.

% INPUTS
% t = time vector
% cpcs = output vector from simulation
% disc = discretization settings
% csint = concentration at interface
load(filename)

ss = disc.ss;
asteps = disc.asteps;
csteps = disc.csteps;
len = ss+asteps+csteps;
L = disc.catlen;
dx = L/csteps;

csint = 0.5;

dims = size(cpcs);
tsteps = dims(1);
rgpos = zeros(tsteps,1);

cs1vt = cpcs(:,2*len+asteps+1:2*len+asteps+csteps);
cs2vt = cpcs(:,2*len+asteps+csteps+1:2*len+asteps+2*csteps);
cs = ((cs1vt+cs2vt)/2)';

% Cycle through time steps
for i=tsteps:-1:1
   % Cycle through  
   for j=1:csteps
       if cs(j,i) >= csint
           rgpos(i) = j*dx;
       end       
   end  
end

% Get data from excel spreadsheet
data = xlsread('gold-red vs time.xls','A2:B54');
xdata = data(:,1);
ydata = data(:,2);

xsc = .077;
xsh = .82;
ysc = 1.23;

% % % % % % % % % poros = 0.3
% xsc = .077;
% xsh = .82;
% ysc = 1.23;


% Plot vs the data

close all
figure
plot(xdata,ydata,'ro', ...
        t*xsc+xsh, rgpos*ysc,'b-','LineWidth',2)

xlabel('t/t_0','FontSize',14)
ylabel('Red/Gold interface position (mm)','FontSize',14)
legend('Data','Model','Location','Northwest')
axis([.8 3 0 1.2])
set(gca,'FontSize',14)

L1 = ysc*.001;
T1 = xsc*28694;
D1 = L1^2/T1;
disp(['D1 = ',num2str(D1)])







return;


function val = calcerr(var,xdata,ydata,t,rgpos,yscale)

val = 0;
xsc = var(1);
xsh = var(2);
xlen = max(size(xdata));

rgpossc = yscale.*rgpos;
tshsc = t.*xsc+xsh;

for i=1:xlen
    rgpostmp = interp1q(tshsc,rgpossc,xdata(i));
    val = val + (rgpostmp-ydata(i))^2;
end

return;