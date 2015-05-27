function res = drawplates(soc, t,cpcs,ffvec,vvec,disc,part, fignum)
% load rotationmatrix
partsize=sqrt(length(part.sizes));
rotationangle=zeros(partsize,partsize)
for index = 1:length(soc)
    clf(figure(fignum))
        hold on

    socindex=soc(index)*2;
    tlen = max(size(t));
    
    Nx = disc.Nx;
    numpart = disc.numpart;
    Ny = numpart;
    
    ss = disc.ss;
    ssx = disc.ss / 1;
    
    e = 1.6e-19;
    k = 1.381e-23;
    T = 298;
    
     x=[0, .75, 1.18, 1.92, 1.18, .75]-1.92/2;
    y=[.5, 1, 1, .5, 0, 0]-.5;
    
    phi = zeros(tlen,Ny,ssx+Nx);
    cs = cell(numpart,1);
    % Put each particle in its own entry in a cell array
    for i=1:numpart*Nx
        % Get the indices
        ind1 = sum(part.steps(1:i))+1;
        ind2 = sum(part.steps(1:i+1));
        
        csmean(i) = mean(cpcs(socindex,disc.sol+ind1-1:disc.sol+ind2-1));
    end
    
    sizenorm = .7+0*part.sizes/max(part.sizes);
    for j = 1:Nx
        for k = 1:numpart
            color = [csmean((j-1)*numpart+k), 1-csmean((j-1)*numpart+k), 0]/(max(csmean((j-1)*numpart+k), 1-csmean((j-1)*numpart+k)));
            
            %rotation matrix
            angle = rotationangle(j,k)*0.5;
% angle = randn(1)
            R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
            matrix = R*[x; y];
            xrot = matrix(1,:);
            yrot=matrix(2,:);
            
            if (csmean((j-1)*numpart+k) < 0.8) && (csmean((j-1)*numpart+k)>0.2)
                fill(xrot+2*j, yrot+k*2, color, 'edgecolor', 'y', 'linewidth', 2)
%                 rectangle('Position', [1+j 1+k 1.2*sizenorm((j-1)*numpart+k) 1*sizenorm((j-1)*numpart+k)], 'curvature', [.75 .75]', 'Facecolor', color, 'edgecolor', 'y', 'linewidth', 2);
            else
                                fill(xrot+2*j, yrot+k*2, color, 'linewidth', 2)

%                 rectangle('Position', [1+j 1+k 1.2*sizenorm((j-1)*numpart+k) 1*sizenorm((j-1)*numpart+k)], 'curvature', [.75 .75]', 'Facecolor', color, 'linewidth', 2);
            end
            
        end
    end
        xlim([0, 22])
    ylim([0, 22])
    pause(0.1)

    %             thresholded =csmattmp;
end