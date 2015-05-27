function cnew = SingleParticle( dt, c, x, SF, D, B, Omega, kappa )

boundary = [x(1); (x(1 : end - 1) + x(2 : end)) / 2; x(end)];
V = 4 * pi * (boundary(2 : end).^3 - boundary(1 : end - 1).^3) / 3;
diagonal = sparse(3/4 * V);
subdiag = sparse([1/4 * V(1); 1/8 * V(2 :end - 1)]);
supdiag = sparse([1/8 * V(2 :end - 1); 1/4 * V(end)]);

MassMatrix = diag(diagonal) + diag(subdiag, -1) + diag(supdiag, 1);
%tmp = full(MassMatrix/(4*pi)*1e3)
%tmp = V/(4*pi)
%zzz



%options = odeset('Mass', MassMatrix, 'MaxStep', dt/50, 'MaxOrder', 1, ...
%options = odeset('Mass', MassMatrix, 'MaxOrder', 1, ...
options = odeset('Mass', MassMatrix, ...
    'RelTol', 1e-6,'MassSingular','no');

[~,cnew] = ode15s(@(t,c)DCDT(t, c, x, SF, D, B, Omega, kappa), [0 dt], c, options);


cnew = cnew(end, :)';

end
