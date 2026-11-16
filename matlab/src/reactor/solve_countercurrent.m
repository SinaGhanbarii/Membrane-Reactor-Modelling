function sol = solve_countercurrent(p, c)
% Counter-current with a sweep gas. The permeate stream leaves at z = 0 and is empty of H2 at z = L,
% so the H2 flow G0 in the permeate stream at z = 0 is found by shooting until G(L) = 0.
Gmax = c.F_in(4) + c.F_in(1);
miss = @(G0) endpoint(p, c, G0);
lo = 1e-12 * Gmax;
hi = 2 * Gmax;   % the root is at most Gmax; the margin keeps the sign change when recovery is almost 100 %
if miss(lo) > 0 || miss(hi) < 0
    error('solve_countercurrent:bracket', 'No sign change in G(L) between %.3g and %.3g mol/s.', lo, hi);
end
G0 = fzero(miss, [lo, hi], optimset('TolX', 1e-15));
sol = integrate_bed(p, c, G0);
if numel(sol.z) < 201
    error('solve_countercurrent:stopped', 'Integration stopped at z = %.4g m of %.4g m.', sol.z(end), c.L);
end
end

function g = endpoint(p, c, G0)
s = integrate_bed(p, c, G0);
g = s.Fp(end);
end
