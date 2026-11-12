function xi = wgs_extent_eq(F0, T)
%WGS_EXTENT_EQ  Equilibrium extent xi [mol/s] of WGS alone (CH4, N2 inert).
%   Solves  K (a - xi)(d - xi) = (b + xi)(c + xi)  with a = CO, d = H2O,
%   b = CO2, c = H2 (inlet flows, mol/s). Exact quadratic; independent of
%   total pressure because delta n = 0 (ideal gas).

k  = pbmr_const();
F0 = F0(:);
a = F0(k.idx.CO);   d = F0(k.idx.H2O);
b = F0(k.idx.CO2);  c = F0(k.idx.H2);

K = Keq_nasa(T);
A = K - 1;
B = -(K*(a + d) + (b + c));
C = K*a*d - b*c;

lo = -min(b, c);
hi =  min(a, d);

if abs(A) < 1e-12
    roots_xi = -C / B;
else
    disc = sqrt(B*B - 4*A*C);
    roots_xi = [(-B - disc)/(2*A), (-B + disc)/(2*A)];
end

tol = 1e-12 * max(abs(F0));
valid = roots_xi(roots_xi >= lo - tol & roots_xi <= hi + tol);
if numel(valid) ~= 1
    error('wgs_extent_eq:root', ...
          'Expected one physical root, got [%g %g] in [%g, %g].', ...
          roots_xi(1), roots_xi(end), lo, hi);
end
xi = valid;
end
