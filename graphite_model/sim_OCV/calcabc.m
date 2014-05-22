function [a, b, V0] = calcabc(c)

close all

e = 1.602e-19;
k = 1.381e-23;
T = 298;

% Constants to fit 
xm = .04;
dV = 1.36245013;
% c = 40;

% Initial guess
x0 = [3.4; 0.5];
% VOltage plateaus
V12 = 4.67125758;
V24 = 3.30880745;

ab = fsolve(@solver,x0,[],dV,xm,c);

% Now get V0
a = ab(1);
b = ab(2);
mu1 = 2*(1+2*exp(-a)*a-exp(-a)*b-exp(-a)*c)*exp(-a)*b;
mu2 = -2*(-1+exp(-a)+2*exp(-2*a)*a-exp(-2*a)*b-exp(-2*a)*c)*b;
V0 = V12 + mu1;
%V02 = V24 + mu2;

% Plot the voltage curves
% plotvoltages(a,b,c,V0);
plotvoltageswithfreeenergy(a,b,c,V0)


return;

function val = solver(x,dV,xm,c)

a = x(1);
b = x(2);
mu1 = 2*(1+2*exp(-a)*a-exp(-a)*b-exp(-a)*c)*exp(-a)*b;
mu2 = -2*(-1+exp(-a)+2*exp(-2*a)*a-exp(-2*a)*b-exp(-2*a)*c)*b;

val = zeros(2,1);
val(1) = exp(-a)*(1+(2*a-b-c)*exp(-a)) - xm;
val(2) = mu2-mu1-dV;

return;

function plotvoltages(a,b,c,V0)

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
%figure('Position',[1 1 3*scrsz(3)/8 scrsz(4)])
% subplot(2,1,1)
hold on
plot(x,V2,'--b','LineWidth',2)
plot(x,Vsep2,'-b','LineWidth',2)
plot(rd(:,1),rd(:,2),'or')
legend('Homogeneous','Phase Separating','Data','Location','North')
ylabel('Voltage, V','FontSize',12)
xlabel('Filling Fraction, x','FontSize',12)
set(gca,'FontSize',12)
hold off
% subplot(2,1,2)
% set(gcf, 'renderer', 'zbuffer');
%contourf(x1,x2,ggrid,40)
% surf(x1,x2,ggrid)
%caxis([-.01 .75])
% xlabel('Layer 1 Filling Fraction','FontSize',12)
% ylabel('Layer 2 Filling Fraction','FontSize',12)
% colorbar('Location','EastOutside')
% set(gca,'FontSize',12)

return;

function plotvoltageswithfreeenergy(a,b,c,V0)

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
plot(x,V2,'--b','LineWidth',2)
plot(x,Vsep2,'-b','LineWidth',2)
plot(rd(:,1),rd(:,2),'or')
legend('Homogeneous','Phase Separating','Data','Location','North')
ylabel('Voltage, V','FontSize',12)
xlabel('Filling Fraction, x','FontSize',12)
set(gca,'FontSize',12)
hold off
subplot(2,1,2)
set(gcf, 'renderer', 'zbuffer');
% contourf(x1,x2,ggrid,40)
surf(x1,x2,ggrid)
% caxis([-.01 .75])
xlabel('Layer 1 Filling Fraction','FontSize',12)
ylabel('Layer 2 Filling Fraction','FontSize',12)
colorbar('Location','EastOutside')
set(gca,'FontSize',12)

return;