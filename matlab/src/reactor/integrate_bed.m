function sol = integrate_bed(p, c, G0)
% Integrate the bed from z = 0 to L with G0 = H2 in the permeate stream at z = 0.
y0 = [c.F_in; G0; c.T_in; c.P_in];
z = linspace(0, c.L, 201);
tol = odeset('RelTol', 1e-8, 'AbsTol', [1e-13 * ones(1, 7), 1e-6, 1e-1], 'NonNegative', 1:6, 'Events', @pressure_event);
[z, y] = ode15s(@(z, y) reactor_rhs(z, y, p, c), z, y0, tol);
sol = struct('z', z, 'F', y(:, 1:6), 'Fp', y(:, 7), 'T', y(:, 8), 'P', y(:, 9));
end

function [value, stop, direction] = pressure_event(~, y)
% Stop the run if the pressure drop brings the retentate down to 1 kPa; the solvers then report where it stopped.
value = y(9) - 1e3;  stop = 1;  direction = 0;
end
