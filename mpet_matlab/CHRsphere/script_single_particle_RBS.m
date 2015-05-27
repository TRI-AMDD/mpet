tic;
clear all;
% close all;

sdisc = 1e-9; % nm
%Grid_Num = 16;

Rp = 100e-9; % m
Grid_Num = 1 + ceil(Rp/sdisc);
B = 0;
Omega = 115.855*10^-3; % eV
kappa = 3.13034*10^9; % eV/m
Va = 3.422; % V
D = 5.3e-19; % m^2/s
%kBT = 8.617*10^-5 * 298.15; % eV
kBT = 8.62047e-5 * 298; % eV
e = 1.602e-19; % C
NA = 6.022 * 10^23; 
Rho_s = 1.3793 * 10^28; % m^-3
Cf = Rho_s / NA; % mol/m^3
C0 = 0.01*Cf; % mol/m^3
capfrac = 0.985;

c_rate = 1E-3;
I = Cf * 4/3* pi *Rp^3 * 96485 / 3600; % A
I = I * c_rate;
k = 1.6E-1; % A/m^2

timestep = 200 + 0;
%dtau = 5 * 10^2; 
%dtau = 5 * 10^-3; 


if I < 0
    C0 = Cf - C0;
end

kappa = kappa / (Rp^2 * Rho_s * kBT);
Omega = Omega / kBT;



tau = Rp^2 / D;
I = I /((96485) / tau);
ie = k;
k = 4 * pi * Rp^4 /96485/D * k;

SF = I / (4 * pi * Rp^3 * Cf);
x = linspace(0, 1, Grid_Num)' ;
%c = C0 * ones(size(x)) / Cf; 
c = 0.01 + 0.005*sin(x*2*pi/0.5)
zzz
D = D/Rp^2*tau;

%if I > 0
%timestep = floor((1-C0/Cf)/(SF * 3)/dtau)-1
%else if I < 0
%        timestep = floor((-C0/Cf)/(SF * 3)/dtau)-1
%    end
%end
%if I == 0
%    timestep = 1000;
%end
dtau = (capfrac)/(3*SF*(timestep));
%dtau*tau*timestep
%zzz


FF = zeros(timestep + 1, 1);
Phi = zeros(timestep + 1, 1);

[~, mu] = ChemicalPotential( c, x, B, Omega, kappa, x(2)-x(1) );
del_Phi = -2 * asinh(I / (k * (1 - c(end)) * exp(mu(end) / 2)) / 2) - mu(end);

% display(['The nondimensional current is ' num2str(I /(k * 0.5))])

concentration = zeros(length(c),timestep + 1);
concentration( : , 1) = c;
FF(1) = C0/Cf;
Phi(1) = del_Phi;

%%

dr = [];
tvec = zeros(timestep + 1, 1);
tvec(1) = 0;
for i = 1 : timestep 

    
    if (del_Phi < -62 || del_Phi > 62) && (i > 2)
        i = i - 1;
        break;
    end
    
    [~, mu] = ChemicalPotential( c, x, B, Omega, kappa, x(2)-x(1) );
    del_Phi = -2 * asinh((I) / (k * (1 - c(end)) * exp(mu(end) / 2)) / 2) - mu(end);
   
    FF(i + 1) = C0/Cf + i * 3 * dtau * SF;
    Phi(i + 1) = del_Phi;
    c = SingleParticle( dtau, c, x, SF, D, B, Omega, kappa );  
    t = i * dtau * tau;
    concentration( : , i + 1) = c;
    tvec(i + 1) = t;
    
    
%     figure(3)
%     if c(end) < 0.05
%         dr = [dr; 0];
%     else
%         iii = 0;
%         while c(Grid_Num - iii) > 0.05 && iii < Grid_Num - 1
%                 iii = iii + 1;
%         end
%         iii
%         dr = [dr ; iii / Grid_Num];
%     end
%     tvec = [tvec;t];
%     plot(sqrt(tvec) , dr , '-b','linewidth', 2);
%     hold off;
%     refresh;
%     pause(.01)   
    

%     figure(1)
%     plot(x , c , '-','linewidth', 2);
%     ylim([0, 1])
%     title(['t = ' num2str(dtau * i * tau) 's, C = ' num2str(c(end) * Cf)],'fontsize',14)
%     xlabel('Radius','fontsize',14);
%     ylabel('Filling Fraction','fontsize',14);
%     refresh;
%     pause(.01)
    disp([num2str(i/timestep*100) ' percent done']);

 
end
%%
%save('out0.mat', 'Phi', 'concentration', 'FF', 'c_rate')
%save('out0.mat')
voltage = Va + kBT * Phi(1:i+1);
concentration = concentration(:, 1 : i + 1);
FF = FF(1 : i + 1);
save('out1.mat')

%f1 = figure (2);
%plot(FF, voltage, 'linewidth', 2);
%xlim([0, 1])
%xlabel('Filling Fraction X','fontsize',14);
%ylabel('Voltage (V)','fontsize',14);
%% title('Voltage-Filling Fraction Plot', 'fontsize',14)
%% title('Voltage-Filling Fraction Plot with Different Surface Wetting', 'fontsize',14)
%hold all;

%%%
%% close all;
%% 
%f2 = figure(3);
%subplot(2,2,4)
%colormap('default')
%% xset('colormap',jetcolormap(128));
%contourf(x, FF, concentration',200);
%shading flat
%grid on;
%% colorbar;
%xlabel('Position r','fontsize',14);
%ylabel('Filling Fraction X','fontsize',14);
%% title(['Concentration Profile during Ion Filling ( I/I_0 = ',...
%%     num2str(I /(ie * 0.5)), ' )']...
%%     ,'fontsize',14)
%
%title(['C rate = 1E3'], 'fontsize',14)
%% title(['Concentration Profile during Ion Filling ( I/I_0 = ',...
%%     num2str(I /(k * 0.5)), ', \nablac|_{surface} = ' num2str(surface_strain(c))...
%%     ', \Omega = ' num2str(Omega) ' )']...
%%     ,'fontsize',24)

% hold on
% hBar = colorbar;
% get(hBar, 'Position');
% set(hBar, 'Position',[0.935    0.08    0.015    0.88])

% % 
% % 
% % %%
% f4 = figure(4);
% 
% angle = linspace(0, 2 * pi , 200);
% 
% X = x * sin(angle);
% Y = x * cos(angle);
% 
% subplot(2,2,1)
% index = 1;
% Z = concentration(:, index) * ...
%     ones(1,length(angle));
% surf(X,Y,Z);
% zlabel('Concentration','fontsize',14)
% shading flat
% zlim([0 1]);
% caxis([0 1])
% title(['X = ', num2str(FF(index))], 'fontsize',14);
% 
% 
% subplot(2,2,2)
% index = floor(size(concentration,2) * 1/4);
% Z = concentration(:, index) * ...
%     ones(1,length(angle));
% surf(X,Y,Z);
% zlabel('Concentration','fontsize',14)
% shading flat
% zlim([0 1]);
% caxis([0 1])
% colormap('default')
% title(['X = ', num2str(FF(index))], 'fontsize',14);
% 
% subplot(2,2,3)
% index = floor(size(concentration,2) * 3/4);
% Z = concentration(:, index) * ...
%     ones(1,length(angle));
% surf(X,Y,Z);
% zlabel('Concentration','fontsize',14)
% shading flat
% zlim([0 1]);
% caxis([0 1])
% title(['X = ', num2str(FF(index))], 'fontsize',14);
% 
% subplot(2,2,4)
% index = size(concentration,2);
% Z = concentration(:, index) * ...
%     ones(1,length(angle));
% surf(X,Y,Z);
% zlabel('Concentration','fontsize',14)
% shading flat
% zlim([0 1]);
% caxis([0 1])
% title(['X = ', num2str(FF(index))], 'fontsize',14);
% 
% h=axes('position',[0.85,0.12,0.05,0.8]);
% caxis([0,1]);
% colorbar
% set(h,'visible','off')
toc
