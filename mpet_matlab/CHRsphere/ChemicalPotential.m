function [mu1, mu2] = ChemicalPotential( c, x, B, Omega, kappa, dx )

mu1 = zeros(size(c));
mu2 = zeros(size(c));

boundary = [x(1); (x(1 : end - 1) + x(2 : end)) / 2; x(end)];
c_mean = c' * (boundary(2 : end).^3 - boundary(1 : end - 1).^3) / x(end)^3;

mu1(1) = B * (c(1) - c_mean) + Omega * (1 - 2 * c(1)) - ...
     3 * kappa * (2 * c(2) - 2 * c(1) ) / (dx^2);
mu1(2 : end - 1) = B * (c(2 : end - 1) - c_mean) + Omega * (1 - 2 * c(2 : end - 1)) - ...
    kappa * ((c(1 : end - 2) - 2 * c(2 : end - 1) + c(3 : end)) / dx^2 + ...
    (c(3 : end) - c(1 : end - 2)) / dx ./ x(2 : end - 1));
mu1(end) = B * (c(end) - c_mean) + Omega * (1 - 2 * c(end)) - ...
    kappa * (2/x(end)* surface_strain(c(end)) +...
    (2*c(end-1) -2*c(end) + 2*dx*surface_strain(c(end)))/dx^2);


mu2(1) = log(c(1) / (1 - c(1))) + mu1(1);
mu2(2 : end - 1) = log(c(2 : end - 1) ./ (1 - c(2 : end - 1))) + mu1(2 : end - 1);
mu2(end) = log(c(end) / (1 - c(end))) + mu1(end);




end

