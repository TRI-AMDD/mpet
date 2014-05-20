function animate(t,cpcs,csmat,disc,psd,ffvec,vvec,fig,output)

close all
tlen = max(size(t));

% total length vector
Lcath = 1;
Ltot = 1+disc.ss/disc.steps;
tlv = linspace(0, Ltot, disc.ss+disc.steps);
% cathode length vector
clv = linspace(0, Lcath, disc.steps);

figure
dims = size(csmat);

if strcmp(fig,'s')
    % Plot the solid concentrations
    for i=1:tlen
        csmattmp = squeeze(csmat(i,:,:));
        surf(csmattmp);    
        axis([0 dims(3) 0 dims(2) 0 1])
        M(i) = getframe(gcf);
    end
elseif strcmp(fig,'c')
    % Plot the electrolyte concentration
    for i=1:tlen
        plot(tlv, cpcs(i,1:disc.ss+disc.steps));
%        axis([0 disc.ss+disc.steps 0 2])
        axis([0 Ltot 0 2])
        M(i) = getframe(gcf);
    end
elseif strcmp(fig, 'cathp')
    % Plot the cathode potential
    for i=1:tlen
        ind1 = disc.sol+disc.totalpart;
        ind2 = disc.sol+disc.totalpart + disc.steps - 1;
        plot(clv,cpcs(i,ind1:ind2),'LineWidth',1)
        axis([clv(1) clv(end) -1 80])
        ylabel('Dimensionless Cathode Potential')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
    end
end

if(output)
    movie2avi(M,'/home/raymond/movie_homog_5C.avi')    
end


return;
