function K = Keq_nasa(T)
%KEQ_NASA  WGS equilibrium constant from NASA-7 Gibbs energies (dimensionless).
%   K = exp(-dG0/(R T)),  dG0 = sum(nu_i * g0_i).  Because delta n = 0 the
%   value does not depend on the standard-state pressure.

c = pbmr_const();
[~, ~, ~, g0] = nasa_props(T);
dG0 = c.NU' * g0;
K = exp(-dG0 / (c.R * T));
end
