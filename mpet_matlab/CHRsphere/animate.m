function animate(type)
load('out1.mat')
numt = max(size(FF));

if strcmp(type, 'csld')
    for i=1:numt
        cvec = concentration(:, i);
        plot(cvec)
        axis([0 max(size(cvec)) 0 1])
        pause(0.02)
    end

elseif strcmp(type, 'v')
    plot(FF, voltage)
%    plot(FF, Phi)

elseif strcmp(type, 'csurf')
    csurf = concentration(end, :);
    format long e
    tvec(end-6:end)
    csurf(end-6:end)
    tvec
    csurf'
    plot(tvec, csurf)

elseif strcmp(type, 'musurf')
    musurf = zeros(numt, 1);
    for i=1:numt
        cvec = concentration(:, i);
        [~, muvec] = ChemicalPotential(cvec, x, B, Omega, kappa, x(2) - x(1));
        musurf(i) = muvec(end);
    end
    plot(tvec, musurf)

elseif strcmp(type, 'soc')
    FF;
    cmeanvec = zeros(numt, 1);
    for i=1:numt
        cvec = concentration(:, i);
        boundary = [x(1); (x(1 : end - 1) + x(2 : end)) / 2; x(end)];
        c_mean = cvec' * (boundary(2 : end).^3 - boundary(1 : end - 1).^3) / x(end)^3;
        cmeanvec(i) = c_mean;
    end

elseif strcmp(type, 'rhs')
    i = 1;
%    D/tau
%    Omega
    kappaval = kappa
    cvec = concentration(:, i);
    cend = cvec(end-5:end)*1e3
    t = tvec(i);
%    tau
    tnow = t*tau
    dcvec = DCDT_mod(t, cvec, x, SF, D, B, Omega, kappa, tau)/tau;
    dcval = dcvec(end)

end
