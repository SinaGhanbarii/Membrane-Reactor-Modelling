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
if ~isscalar(c.sweep) || c.sweep < 0
    error('make_case:sweep', 'sweep must be a single value >= 0 [mol/s].');
end
if c.sweep == 0, c.flow = 'cocurrent'; end   % without a sweep the direction changes nothing
if strcmp(c.dp, 'on') && isempty(c.mu)
    error('make_case:mu', 'dp = on needs the gas viscosity mu [Pa s].');
end

c.A_bed = pi/4 * (c.D_s^2 - c.D_m^2);          % annulus between shell and membrane, m2
c.rho_cat = (1 - c.eps) * c.rho_p * c.w_cat;   % kg catalyst per m3 bed
end
