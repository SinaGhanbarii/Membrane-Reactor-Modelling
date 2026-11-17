function c = make_case(p, ov)
% Case definition for the reactor model. Defaults come from the register; fields in ov replace them.
if nargin < 2, ov = struct(); end
sp = pbmr_const().species;
c.F_in = zeros(6, 1);
for i = 1:6
    c.F_in(i) = 1e-3 * pget(p, ['cases.target.test_feed_mmol_s.' sp{i}]);   % mol/s
end
c.T_in = pget(p, 'cases.base.T_in');
c.P_in = pget(p, 'cases.base.P_ret');
c.P_perm = pget(p, 'cases.target.P_perm');
c.L = pget(p, 'cases.base.membrane_length');
c.D_m = pget(p, 'cases.base.membrane_OD');
c.D_s = pget(p, 'cases.base.shell_ID');
c.eps = pget(p, 'cases.base.bed_voidage');
c.d_p = pget(p, 'cases.base.d_p');
c.w_cat = pget(p, 'cases.base.catalyst_wt_frac');
c.rho_p = pget(p, 'catalyst.Hla2015.rho_p');
c.membrane = pget(p, 'cases.target.membrane_option');
c.law = 'K1';
c.inhibition = 'off';
c.eta = 'off';   c.D_eff = [];
c.wall = 'isothermal';   c.U = [];   c.T_wall = [];
c.dp = 'off';    c.mu = [];
c.perm_scale = 1;
c.flow = 'cocurrent';   c.sweep = 0;   % sweep = inert sweep gas on the permeate side, mol/s
c.backperm = 'on';   % with a sweep, H2 flows back into the retentate where its pressure is higher on the permeate side

fn = fieldnames(ov);
for k = 1:numel(fn)
    if ~isfield(c, fn{k}), error('make_case:field', 'Unknown case field "%s".', fn{k}); end
    c.(fn{k}) = ov.(fn{k});
end

if ~any(strcmp(c.wall, {'isothermal', 'adiabatic', 'cooled'}))
    error('make_case:wall', 'wall must be isothermal, adiabatic or cooled.');
end
if strcmp(c.wall, 'cooled') && (isempty(c.U) || isempty(c.T_wall))
    error('make_case:wall', 'wall = cooled needs U [W/m2/K] and T_wall [K].');
end
if strcmp(c.eta, 'on') && isempty(c.D_eff)
    error('make_case:eta', 'eta = on needs D_eff [m2/s].');
end
if ~any(strcmp(c.flow, {'cocurrent', 'counter'}))
    error('make_case:flow', 'flow must be cocurrent or counter.');
end
if ~any(strcmp(c.backperm, {'on', 'off'}))
    error('make_case:backperm', 'backperm must be on or off.');
end
if ~any(strcmp(c.eta, {'on', 'off'})) || ~any(strcmp(c.dp, {'on', 'off'}))
    error('make_case:switch', 'eta and dp must each be on or off.');
end
if ~isnumeric(c.sweep) || ~isscalar(c.sweep) || ~isfinite(c.sweep) || c.sweep < 0
    error('make_case:sweep', 'sweep must be a single value >= 0 [mol/s].');
end
if c.sweep == 0, c.flow = 'cocurrent'; end   % without a sweep the direction changes nothing
if strcmp(c.dp, 'on') && isempty(c.mu)
    error('make_case:mu', 'dp = on needs the gas viscosity mu [Pa s].');
end

% Sanity checks: catch wrong units and impossible geometry before they turn into plausible-looking numbers.
c.F_in = c.F_in(:);
if ~isnumeric(c.F_in) || numel(c.F_in) ~= 6 || any(~isfinite(c.F_in)) || any(c.F_in < 0) || c.F_in(1) <= 0
    error('make_case:feed', 'F_in must be 6 finite flows >= 0 in mol/s (CO, H2O, CO2, H2, CH4, N2) with CO > 0.');
end
limits = {'T_in', 250, 1500, 'K (a temperature in C?)'; ...
          'P_in', 1e4, 1e8, 'Pa (a pressure in bar?)'; ...
          'P_perm', 1e3, 1e8, 'Pa (a pressure in bar?)'; ...
          'L', 1e-3, 10, 'm (a length in mm?)'; ...
          'D_m', 1e-4, 1, 'm (a diameter in mm?)'; ...
          'D_s', 1e-4, 1, 'm (a diameter in mm?)'; ...
          'd_p', 1e-5, 5e-2, 'm (a size in mm?)'; ...
          'rho_p', 100, 2e4, 'kg/m3'};
for i = 1:size(limits, 1)
    v = c.(limits{i, 1});
    if ~(isnumeric(v) && isscalar(v) && isfinite(v) && v >= limits{i, 2} && v <= limits{i, 3})
        error('make_case:range', '%s = %s is outside %g to %g %s.', limits{i, 1}, num2str(v), limits{i, 2}, limits{i, 3}, limits{i, 4});
    end
end
if strcmp(c.wall, 'cooled') && ~(isnumeric(c.U) && isscalar(c.U) && isfinite(c.U) && c.U >= 0 && ...
                                 isnumeric(c.T_wall) && isscalar(c.T_wall) && isfinite(c.T_wall) && c.T_wall >= 250 && c.T_wall <= 1500)
    error('make_case:wall', 'U must be a single number >= 0 W/m2/K and T_wall a single temperature of 250 to 1500 K.');
end
if strcmp(c.eta, 'on') && ~(isnumeric(c.D_eff) && isscalar(c.D_eff) && isfinite(c.D_eff) && c.D_eff > 0)
    error('make_case:eta', 'D_eff must be a single number > 0 m2/s.');
end
if strcmp(c.dp, 'on') && ~(isnumeric(c.mu) && isscalar(c.mu) && isfinite(c.mu) && c.mu > 0)
    error('make_case:mu', 'mu must be a single number > 0 Pa s.');
end
if ~(isnumeric(c.eps) && isscalar(c.eps) && c.eps > 0 && c.eps < 1)
    error('make_case:range', 'eps (bed voidage) must be a single number between 0 and 1.');
end
if ~(isnumeric(c.w_cat) && isscalar(c.w_cat) && c.w_cat >= 0 && c.w_cat <= 1)
    error('make_case:range', 'w_cat (catalyst weight fraction) must be a single number between 0 and 1.');
end
if ~(isnumeric(c.perm_scale) && isscalar(c.perm_scale) && isfinite(c.perm_scale) && c.perm_scale >= 0)
    error('make_case:range', 'perm_scale must be a single finite number >= 0.');
end
if c.D_s <= c.D_m
    error('make_case:geometry', 'Shell diameter D_s (%g m) must be larger than membrane diameter D_m (%g m).', c.D_s, c.D_m);
end
if ~any(strcmp(c.law, {'K1', 'G2'}))
    error('make_case:law', 'law must be K1 or G2.');
end
if ~ischar(c.membrane) || ~isfield(p.membrane, c.membrane)
    error('make_case:membrane', 'Unknown membrane option.');
end
if ~ischar(c.inhibition) || ~isfield(p.inhibition, c.inhibition)
    error('make_case:inhibition', 'Unknown inhibition option.');
end
if c.P_perm >= c.P_in
    warning('make_case:pressure', 'Permeate pressure (%.2f bar) is not below retentate pressure (%.2f bar): only a sweep gas can drive permeation.', ...
            c.P_perm / 1e5, c.P_in / 1e5);
end

c.A_bed = pi/4 * (c.D_s^2 - c.D_m^2);          % annulus between shell and membrane, m2
c.rho_cat = (1 - c.eps) * c.rho_p * c.w_cat;   % kg catalyst per m3 bed
end
