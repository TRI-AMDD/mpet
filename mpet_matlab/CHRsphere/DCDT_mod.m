function dc = DCDT_mod( t, c, x, SF, D, B, Omega, kappa, tau )


dx = x(2) - x(1);
boundary = [x(1); (x(2 : end) + x(1 : end - 1)) / 2; x(end)];
F = zeros(size(boundary));

mode = 2;
if mode == 1
[mu, ~] = ChemicalPotential( c, x, B, Omega, kappa, dx );

F(1) = 0;
F(2 : end - 1) = D * (c(2 : end) - c(1 : end - 1)) / dx + ...
    D * (1- ((c(2 : end) + c(1 : end - 1)) / 2)) .* ...
    ((c(2 : end) + c(1 : end - 1)) / 2) .* (mu(2 : end) - mu(1 : end - 1)) / dx;
F(end) = SF;

end

if mode == 2
[~, mu] = ChemicalPotential_mod( c, x, B, Omega, kappa, dx );

F(1) = 0;
F(2 : end - 1) = D * (1- ((c(2 : end) + c(1 : end - 1)) / 2)) .* ...
    ((c(2 : end) + c(1 : end - 1)) / 2) .* (mu(2 : end) - mu(1 : end - 1)) / dx;
F(end) = SF;
end

Fval = F(end)/tau
F = F .* boundary.^2 * 4 * pi;
muval = mu(end)
%F(end-1)/tau
%boundary * 1e3
%boundary.^2 * 1e3
%zzz

dc = (F(2 : end) - F(1 : end - 1));
%dc(end) = log(-1);


end

