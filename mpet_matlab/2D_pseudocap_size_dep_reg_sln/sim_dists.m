function sim_dists()

% We generate data for multiple size distributions

mu = [28];
sigma = [50,80,100];
%sigma = [0];

mulen = max(size(mu));
siglen = max(size(sigma));

for i=1:mulen
    for j=1:siglen
        gaberscek_sim2(mu(i),sigma(j));
    end
end

return;
