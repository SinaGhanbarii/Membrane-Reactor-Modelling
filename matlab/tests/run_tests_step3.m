function run_tests_step3()
%RUN_TESTS_STEP3  Checks for the co-current reactor model (make_case, reactor_rhs, solve_cocurrent, run_case).
%   Run:  >> clear functions; setup_paths; run_tests_step3
%   Reference numbers come from an independent Python model (Cantera thermo, SciPy BDF solver,
%   same register values): scripts/ref_reactor.py.
fprintf('=== pbmr1d MATLAB - Step 3 tests (reactor model) ===\n');
np = 0; nf = 0;
p = load_params();
k = pbmr_const();

% 1. base case against the independent Python model
b = run_case(p);
s = b.summary;
[np, nf] = check(np, nf, relerr(s.X_CO, 0.94271687) < 5e-4, sprintf('isothermal base case: X_CO = 0.94272 (got %.5f)', s.X_CO));
[np, nf] = check(np, nf, relerr(b.sol.Fp(end), 2.93392991e-4) < 5e-4, sprintf('isothermal base case: permeated H2 = 0.29339 mmol/s (got %.5f)', 1e3 * b.sol.Fp(end)));
a = run_case(p, struct('wall', 'adiabatic')).summary;
[np, nf] = check(np, nf, relerr(a.X_CO, 0.82940307) < 5e-4, sprintf('adiabatic base case: X_CO = 0.82940 (got %.5f)', a.X_CO));
[np, nf] = check(np, nf, relerr(a.T_out, 823.804164) < 1e-4, sprintf('adiabatic base case: T_out = 823.80 K (got %.2f)', a.T_out));

% 2. conservation and limits
[np, nf] = check(np, nf, s.atom_error < 1e-6, sprintf('C, H, O, N balance closes (relative error %.1e)', s.atom_error));
[np, nf] = check(np, nf, max(abs(b.sol.T - b.case.T_in)) < 1e-9, 'isothermal wall keeps T at the inlet value');
nm = run_case(p, struct('perm_scale', 0));
[np, nf] = check(np, nf, abs(nm.sol.Fp(end)) < 1e-15, 'no membrane: nothing permeates');
[np, nf] = check(np, nf, abs(nm.summary.X_CO - nm.summary.X_eq) < 1e-3, sprintf('no membrane, long bed: X_CO reaches WGS equilibrium (%.5f vs %.5f)', nm.summary.X_CO, nm.summary.X_eq));
[np, nf] = check(np, nf, s.X_CO > s.X_eq, 'with the membrane, X_CO exceeds the no-membrane equilibrium');
[np, nf] = check(np, nf, s.H2_recovery > 0 && s.H2_recovery < 1, sprintf('H2 recovery between 0 and 1 (got %.3f)', s.H2_recovery));

% 3. adiabatic energy balance: enthalpy in = enthalpy out (no membrane)
na = run_case(p, struct('perm_scale', 0, 'wall', 'adiabatic'));
[~, h_in] = nasa_props(na.case.T_in);
[~, h_out] = nasa_props(na.summary.T_out);
Fo = na.sol.F(end, :)';
extent = na.case.F_in(1) - Fo(1);
dHin = na.case.F_in' * h_in;
dHout = Fo' * h_out;
[np, nf] = check(np, nf, abs(dHout - dHin) < 1e-3 * abs(dH_wgs(na.case.T_in) * extent), 'adiabatic, no membrane: enthalpy in = enthalpy out');
[np, nf] = check(np, nf, relerr(na.summary.T_out, 814.310177) < 1e-4, sprintf('adiabatic, no membrane: T_out = 814.31 K (got %.2f)', na.summary.T_out));

% 4. permeation wiring
[np, nf] = check(np, nf, b.profile.perm(1) == 0, 'inlet H2 partial pressure is below the permeate pressure: no flux at z = 0');
[np, nf] = check(np, nf, all(diff(b.sol.Fp) >= -1e-15), 'permeated H2 never decreases along the bed');
i = 100;
yi = b.sol.F(i, :)' / sum(b.sol.F(i, :));
Ji = membrane_flux(p, 'M7', 'off', b.sol.T(i), yi(4) * b.sol.P(i), b.case.P_perm, yi(1) * b.sol.P(i) / 1e3);
[np, nf] = check(np, nf, b.profile.perm(i) > 0 && relerr(b.profile.perm(i), Ji * pi * b.case.D_m) < 1e-10, 'mid-bed permeation equals membrane_flux times membrane perimeter');
hr = [0.4; 1.2; 0.15; 1.0; 0.05; 0.05] * 1e-3;
zr = run_case(p, struct('w_cat', 0, 'F_in', hr));
Fo = zr.sol.F(end, :)';
[np, nf] = check(np, nf, abs(zr.summary.X_CO) < 1e-9 && zr.sol.Fp(end) > 0, 'no catalyst: no CO converted, H2 still permeates');
[np, nf] = check(np, nf, relerr(Fo(4) + zr.sol.Fp(end), hr(4)) < 1e-8, 'no catalyst: H2 in = H2 out + H2 permeated');
y0 = hr / sum(hr);
J0 = membrane_flux(p, 'M7', 'off', zr.case.T_in, y0(4) * zr.case.P_in, zr.case.P_perm, y0(1) * zr.case.P_in / 1e3);
[np, nf] = check(np, nf, relerr(zr.profile.perm(1), J0 * pi * zr.case.D_m) < 1e-10, 'no catalyst: inlet permeation equals membrane_flux by hand');

% 4b. local regimes along the bed (Python reference: 3.5 / 40.8 / 34.8 / 20.9 % of the bed length)
rf = s.regime_fraction(1:4);   % the fifth label (idle) does not occur without a sweep
ref = [0.035; 0.408; 0.348; 0.209];
[np, nf] = check(np, nf, abs(sum(s.regime_fraction) - 1) < 1e-12 && s.regime_fraction(5) == 0 && all(abs(rf - ref) < 0.01), sprintf('base case regime fractions no flux / reaction / pinch / permeance = %.3f / %.3f / %.3f / %.3f', rf));
[np, nf] = check(np, nf, b.profile.regime(1) == 1 && b.profile.regime(end) == 3, 'base case: bed starts with no flux and ends at the pinch');
[np, nf] = check(np, nf, all(isnan(nm.summary.regime_fraction)) && all(nm.profile.regime == 0), 'no membrane: no regime labels');

% 5. pore diffusion, wall cooling, pressure drop, rate law
sb = struct('perm_scale', 0, 'L', 0.01);
off = run_case(p, sb).summary;
sb.eta = 'on';  sb.d_p = 3.2e-3;  sb.D_eff = 1e-7;
on = run_case(p, sb).summary;
[np, nf] = check(np, nf, relerr(off.X_CO, 0.234649) < 5e-4, sprintf('short bed, eta off: X_CO = 0.23465 (got %.5f)', off.X_CO));
[np, nf] = check(np, nf, on.X_CO > 0 && on.X_CO < off.X_CO - 1e-3, sprintf('short bed, eta on with slow diffusion lowers X_CO (%.4f vs %.4f)', on.X_CO, off.X_CO));
cl = run_case(p, struct('wall', 'cooled', 'U', 1e4, 'T_wall', 650));
[np, nf] = check(np, nf, abs(cl.summary.T_out - 650) < 0.05, sprintf('strongly cooled wall: T_out follows T_wall (got %.3f K)', cl.summary.T_out));
cp = make_case(p, struct('dp', 'on', 'mu', 3e-5));
dy = reactor_rhs(0, [cp.F_in; 0; cp.T_in; cp.P_in], p, cp);
Mi = [28.010 18.015 44.009 2.016 16.043 28.014]' / 1000;
yf = cp.F_in / sum(cp.F_in);
rho = cp.P_in * (yf' * Mi) / (k.R * cp.T_in);
u = sum(cp.F_in) * k.R * cp.T_in / (cp.P_in * cp.A_bed);
dp = -(150 * cp.mu * (1 - cp.eps)^2 / (cp.eps^3 * cp.d_p^2) * u + 1.75 * (1 - cp.eps) / (cp.eps^3 * cp.d_p) * rho * u^2);
[np, nf] = check(np, nf, relerr(dy(9), dp) < 2e-3, sprintf('Ergun gradient at the inlet matches the hand calculation (%.3g Pa/m)', dy(9)));
g2 = run_case(p, struct('law', 'G2')).summary;
[np, nf] = check(np, nf, isfinite(g2.X_CO) && g2.X_CO > 0 && g2.X_CO < 1, sprintf('G2 rate law runs (X_CO = %.4f)', g2.X_CO));
[np, nf] = check(np, nf, any(contains(g2.flags, 'EXTRAPOLATED')) && any(contains(s.flags, 'EXTRAPOLATED')), 'K1 and G2 runs at 40 % CO are flagged EXTRAPOLATED');

% 6. input checks
[np, nf] = check(np, nf, stops(@() make_case(p, struct('foo', 1)), 'make_case:field'), 'unknown case field stops with make_case:field');
[np, nf] = check(np, nf, stops(@() make_case(p, struct('wall', 'cooled')), 'make_case:wall'), 'cooled wall without U and T_wall stops with make_case:wall');
[np, nf] = check(np, nf, stops(@() make_case(p, struct('eta', 'on')), 'make_case:eta'), 'eta on without D_eff stops with make_case:eta');
[np, nf] = check(np, nf, stops(@() make_case(p, struct('dp', 'on')), 'make_case:mu'), 'pressure drop on without mu stops with make_case:mu');

% 7. sweep gas and counter-current (reference values from scripts/ref_reactor.py)
b = run_case(p).summary;
nc = run_case(p, struct('flow', 'counter')).summary;
[np, nf] = check(np, nf, abs(nc.X_CO - b.X_CO) < 1e-9 && abs(nc.H2_recovery - b.H2_recovery) < 1e-9, 'counter flag without a sweep gives the co-current result');
co1 = run_case(p, struct('sweep', 1e-3)).summary;
[np, nf] = check(np, nf, abs(co1.X_CO - 0.981360) < 2e-5 && abs(co1.H2_recovery - 0.809386) < 2e-5, sprintf('co-current, 1 mmol/s sweep: X_CO %.6f, recovery %.6f', co1.X_CO, co1.H2_recovery));
cc1 = run_case(p, struct('sweep', 1e-3, 'flow', 'counter'));
[np, nf] = check(np, nf, abs(cc1.summary.X_CO - 0.999973) < 1e-5 && abs(cc1.summary.H2_recovery - 0.999986) < 1e-5, sprintf('counter-current, 1 mmol/s sweep: X_CO %.6f, recovery %.6f', cc1.summary.X_CO, cc1.summary.H2_recovery));
[np, nf] = check(np, nf, abs(cc1.sol.Fp(end)) < 1e-9 * (cc1.case.F_in(4) + cc1.case.F_in(1)), 'counter-current permeate stream is empty of H2 at z = L');
[np, nf] = check(np, nf, cc1.summary.atom_error < 1e-6, sprintf('counter-current atom balance closes (%.2g)', cc1.summary.atom_error));
cc5 = run_case(p, struct('sweep', 5e-3, 'flow', 'counter')).summary;
co5 = run_case(p, struct('sweep', 5e-3)).summary;
[np, nf] = check(np, nf, abs(co5.H2_recovery - 0.941399) < 2e-5 && cc5.H2_recovery > co5.H2_recovery, sprintf('5 mmol/s sweep: co-current %.4f, counter-current %.4f recovery', co5.H2_recovery, cc5.H2_recovery));
[np, nf] = check(np, nf, abs(cc1.summary.permeate_H2_fraction - 7.499892e-4 / (7.499892e-4 + 1e-3)) < 1e-5, 'H2 share of the permeate stream follows from recovery and sweep');
[np, nf] = check(np, nf, stops(@() make_case(p, struct('flow', 'sideways')), 'make_case:flow'), 'invalid flow stops with make_case:flow');
[np, nf] = check(np, nf, stops(@() make_case(p, struct('sweep', -1)), 'make_case:sweep'), 'negative sweep stops with make_case:sweep');
idl = cc1.summary.regime_fraction(5);
[np, nf] = check(np, nf, abs(idl - 0.438) < 0.03, sprintf('counter-current, 1 mmol/s sweep: %.3f of the bed is idle (Python 0.438)', idl));

fprintf('--- %d passed, %d failed ---\n', np, nf);
if nf == 0, fprintf('ALL STEP 3 TESTS PASSED\n'); else, error('Step 3 tests failed'); end
end

function [np, nf] = check(np, nf, ok, msg)
if ok, np = np + 1; fprintf('  PASS  %s\n', msg);
else,  nf = nf + 1; fprintf('  FAIL  %s\n', msg); end
end

function e = relerr(a, b), e = abs(a - b) ./ abs(b); end

function ok = stops(f, id)
ok = false;
try
    f();
catch err
    ok = strcmp(err.identifier, id);
end
end
