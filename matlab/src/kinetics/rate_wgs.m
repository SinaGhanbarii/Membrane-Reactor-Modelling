function [r, info] = rate_wgs(p, law, T, P, y, opts)
% Net WGS rate in mol per kg catalyst per s.
%   law: 'K1' (Keiski 1996) or 'G2' (Hla 2009 form)     T [K]    P [Pa]
%   y:   mole fractions [CO H2O CO2 H2 CH4 N2]
%   opts.eta = 'on' multiplies the rate by a pore-diffusion factor and needs
%   opts.d_p [m] and opts.D_eff [m2/s]; opts.rho_p [kg/m3] is optional.

if nargin < 6, opts = struct(); end
if ~isfield(opts, 'eta'), opts.eta = 'off'; end
if ~any(strcmp(opts.eta, {'on', 'off'})), error('rate_wgs:eta', 'opts.eta must be ''on'' or ''off''.'); end
R = pbmr_const().R;
y_in = y(:);
y = max(y_in, 1e-6);                  % G2 has negative orders in CO2 and H2, so keep them above zero

beta = y(3) * y(4) / (y(1) * y(2) * Keq_nasa(T));       % equals 1 at WGS equilibrium

switch law
    case 'K1', r0 = rate_K1(p, T, P, y, R);
    case 'G2', r0 = rate_G2(p, T, P, y, R);
    case 'G1', error('rate_wgs:G1', 'G1 is not available: its unit basis is unresolved (scripts/check_rate_units.py). Use K1.');
    case 'G3', error('rate_wgs:G3', 'G3 is not implemented.');
    otherwise, error('rate_wgs:set', 'Unknown rate law "%s". Use K1 or G2.', law);
end
r0 = r0 * (1 - beta);

eta = 1;  phi = 0;
if strcmp(opts.eta, 'on')
    if ~isfield(opts, 'd_p') || ~isfield(opts, 'D_eff')
        error('rate_wgs:eta', 'eta = on needs opts.d_p [m] and opts.D_eff [m2/s].');
    end
    rho_p = pget(p, 'catalyst.Hla2015.rho_p');
    if isfield(opts, 'rho_p'), rho_p = opts.rho_p; end
    n_CO = pget(p, ['rate.' law '.n_CO']);
    c_CO = y(1) * P / (R * T);
    r_volume = max(r0, 0) * rho_p;                       % mol/m3/s inside the particle
    phi = opts.d_p / 6 * sqrt((n_CO + 1) / 2 * r_volume / (opts.D_eff * c_CO));
    eta = effectiveness_sphere(phi);
end
r = eta * r0;

y_dry = y_in(1) / max(1 - y_in(2), eps);                % CO on a dry basis
info = struct('r_intrinsic', r0, 'eta', eta, 'phi', phi, 'beta', beta, ...
              'flags', {check_validity(p.rate.(law).validity, T, P / 1e5, y_dry)});
end

function r = rate_K1(p, T, P, y, R)
% Keiski (1996): r = k0 exp(-E/RT) c_CO^n c_H2O^m, before the (1 - beta) factor.
switch pget(p, 'rate.K1.conc_basis_used')    % inferred, see scripts/check_rate_units.py
    case 'mol cm-3', unit = 1e-6;
    case 'mol dm-3', unit = 1e-3;
    case 'mol m-3',  unit = 1;
    otherwise, error('rate_wgs:basis', 'Unknown K1 concentration basis.');
end
c = y * P / (R * T) * unit;
k = exp(pget(p, 'rate.K1.ln_k0') - pget(p, 'rate.K1.E') / (R * T));
r = 1000 * k * c(1)^pget(p, 'rate.K1.n_CO') * c(2)^pget(p, 'rate.K1.m_H2O');     % mol/g/s -> mol/kg/s
end

function r = rate_G2(p, T, P, y, R)
% Hla (2015) Eq. 2: r = 10^a exp(-E/RT) p_CO^n p_H2O^m p_CO2^q p_H2^s, before the (1 - beta) factor.
switch pget(p, 'rate.G2.pressure_unit')      % inferred, see scripts/check_rate_units.py
    case 'kPa', unit = 1e-3;
    case 'Pa',  unit = 1;
    case 'bar', unit = 1e-5;
    otherwise, error('rate_wgs:punit', 'Unknown G2 pressure unit.');
end
pp = y * P * unit;
k = 10^pget(p, 'rate.G2.log10_k0') * exp(-pget(p, 'rate.G2.E') / (R * T));
r = 1000 * k * pp(1)^pget(p, 'rate.G2.n_CO') * pp(2)^pget(p, 'rate.G2.m_H2O') ...
         * pp(3)^pget(p, 'rate.G2.p_CO2') * pp(4)^pget(p, 'rate.G2.q_H2');          % mol/g/s -> mol/kg/s
end
