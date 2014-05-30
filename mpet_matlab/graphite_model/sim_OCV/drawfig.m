function drawfig()

clear

c = 30;
a = 3.4;
b = 1.4;
V0 = 5.3174;
params = struct('cc',c,'ac',a,'bc',b,'V0',V0);

% [t,cpcs,clv,plva,plvc,ffvec,vvec,disc] = pm_gr1dasc_cc_rev2(1,.1,.0001,params);
% save('graphite_ocp.mat')
load('graphite_ocp.mat')

plotvoltageswithfreeenergy(a,b,c,V0,ffvec,vvec);

return;

function plotvoltageswithfreeenergy(a,b,c,V0,ffvec,vvec)

e = 1.602e-19;
k = 1.381e-23;
T = 298;

% Make some plots
x = [0:.001:1];
V = V0 -(2*log(x./(1-x))+(2*b-4*a).*x+2*a+c.*(2.*x.*(1-x).^2-2.*x.^2.*(1-x)));

num = size(x);
Vsep = zeros(1,num(2));
xeq = exp(-a)*(1+(2*a-b-c)*exp(-a));

for i=1:num(2)
    if (x(i)<xeq)
        Vsep(i) = V0 - (2*log(x(i)/(1-x(i)))+(2*b-4*a)*x(i)+2*a ...
                         +c*(2*x(i)*(1-x(i))^2-2*x(i)^2*(1-x(i))));
    elseif ((x(i)>=xeq) && (x(i)<=.5))
        Vsep(i)= V0 - (2*(1+2*exp(-a)*a-exp(-a)*b-exp(-a)*c)*exp(-a)*b);
    elseif ((x(i)>.5) && (x(i) <= (1-xeq)))
        Vsep(i) = V0 - (-2*(-1+exp(-a)+2*exp(-2*a)*a-exp(-2*a)*b-exp(-2*a)*c)*b);
    else
        Vsep(i) = V0 - (2*log(x(i)/(1-x(i)))+(2*b-4*a)*x(i)+2*a ...
                         +c*(2*x(i)*(1-x(i))^2-2*x(i)^2*(1-x(i))));
    end
end

% Contour plot
x1 = [0:.01:1];
num2 = max(size(x1));
x2 = x1;
ggrid = zeros(num2,num2);
for i=1:num2
    for j=1:num2
        ggrid(i,j) = x1(i)*log(x1(i)) + (1-x1(i))*log(1-x1(i)) + ...
                        x2(j)*log(x2(j)) + (1-x2(j))*log(1-x2(j)) + ...
                        a*x1(i)*(1-x1(i)) + a*x2(j)*(1-x2(j)) + ...
                        b*x1(i)*x2(j) + c*x1(i)*(1-x1(i))*x2(j)*(1-x2(j));
    end
end

V2 = V*((k*T)/e);
Vsep2 = Vsep*((k*T)/e);

% Real data
rd = xlsread('ocp.xls','A2:B113');

% figure
% plot(x,V,'--b',x,Vsep,'-b')
% axis([0 1 2 10])
% legend('Homogeneous','Phase Separating','Location','North')
% title('Homogeneous and Phase Separating Voltage Profiles vs. Filling Fraction')
% ylabel('Dimensionless Voltage, Ve/kT')
% xlabel('Filling Fraction, x')

scrsz = get(0,'ScreenSize');
% Position -> Left Bottom Width Height
figure('Position',[1 1 3*scrsz(3)/8 scrsz(4)])
subplot(2,1,1)
hold on
plot(rd(:,1),rd(:,2),'or','LineWidth',2)
plot(x,V2,'--b','LineWidth',2)
plot(ffvec,vvec,'-b','LineWidth',2)
legend('Homogeneous','Phase Separating','Data','Location','North')
ylabel('Voltage, V','FontSize',12)
xlabel('Filling Fraction, x','FontSize',12)
set(gca,'FontSize',12)
hold off
subplot(2,1,2)
% set(gcf, 'renderer', 'zbuffer');
% contourf(x1,x2,ggrid,40)
surf(x1,x2,ggrid)
view([58 7])
% caxis([-.01 .75])
xlabel('Layer 1 Filling Fraction','FontSize',12)
ylabel('Layer 2 Filling Fraction','FontSize',12)
colorbar('Location','EastOutside')
set(gca,'FontSize',12)

return;
