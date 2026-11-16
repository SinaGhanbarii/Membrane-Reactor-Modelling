function sol = solve_cocurrent(p, c)
% Steady co-current solution along the bed; the permeate stream starts empty at z = 0.
sol = integrate_bed(p, c, 0);
if numel(sol.z) < 201
    error('solve_cocurrent:stopped', 'Integration stopped at z = %.4g m of %.4g m.', sol.z(end), c.L);
end
end
