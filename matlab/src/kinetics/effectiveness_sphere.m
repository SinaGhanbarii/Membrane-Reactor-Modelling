function eta = effectiveness_sphere(phi)
% Effectiveness factor of an isothermal sphere, phi based on R/3:
%   eta = (1/phi) * (1/tanh(3 phi) - 1/(3 phi))
% Below phi = 1e-3 the formula loses precision, so the series 1 - (3 phi)^2/15 is used.

eta = ones(size(phi));
small = phi < 1e-3;
eta(small) = 1 - (3 * phi(small)).^2 / 15;
phi_big = phi(~small);
eta(~small) = (1 ./ phi_big) .* (1 ./ tanh(3 * phi_big) - 1 ./ (3 * phi_big));
end
