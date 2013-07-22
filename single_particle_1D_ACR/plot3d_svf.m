function plot3d_svf()

% This function generates a 3D plot of voltage, size and filling fraction
% for the 1D surface wetted particles with stress effects

szs = 20:20:300;
numszs = max(size(szs));

data = cell(numszs,2);
% for i=1:numszs
%     [t,cs,ffvec,vvec] = acr_sp_swse_rev1(.001,szs(i));
%     data{i,1} = real(ffvec);
%     data{i,2} = real(vvec);
% end
% save('data_struct.mat','data');

load('data_struct.mat')

% Now we create the 3D plot
ffvec = .01:.01:.98;
vgrid = zeros(max(size(ffvec)),numszs);

for i=1:max(size(ffvec))
    for j=1:numszs
        vgrid(i,j) = interp1q(data{j,1},data{j,2},ffvec(i));
    end
end

surf(szs,ffvec,vgrid)

return;