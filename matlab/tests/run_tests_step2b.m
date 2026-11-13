function run_tests_step2b()
%RUN_TESTS_STEP2B  Checks for rate_wgs and the effectiveness factor.
%   Run:  >> clear functions; setup_paths; run_tests_step2b
%   Reference rates were computed independently in Python with the Moe (1962) Keq; they agree with
%   the NASA Keq used here to better than 0.2 % at this composition (beta = 0.018).
fprintf('=== pbmr1d MATLAB - Step 2b tests (rate_wgs, effectiveness factor) ===\n');
np = 0; nf = 0;
p = load_params();
y = [0.10 0.30 0.10 0.15 0.05 0.30]; T = 600; P = 5e5;

% 1. absolute values against the independent Python calculation
r1 = rate_wgs(p, 'K1', T, P, y);
r2 = rate_wgs(p, 'G2', T, P, y);
[np, nf] = check(np, nf, relerr(r1, 0.014492901) < 2e-3, sprintf('K1 rate at 600 K, 5 bar = 1.449e-2 mol/kg/s (got %.4e)', r1));
[np, nf] = check(np, nf, relerr(r2, 0.006836777) < 2e-3, sprintf('G2 rate at 600 K, 5 bar = 6.837e-3 mol/kg/s (got %.4e)', r2));

% 2. reversibility: zero at equilibrium, negative beyond it
yH2 = Keq_nasa(T) * 0.1 * 0.3 / 0.2;
ye = [0.1 0.3 0.2 yH2 0 0]; ye = ye / sum(ye);
re = rate_wgs(p, 'K1', T, P, ye);
[np, nf] = check(np, nf, abs(re) < 1e-9, 'K1 rate is zero at the WGS equilibrium composition');
yb = ye; yb(4) = yb(4) * 1.5; yb = yb / sum(yb);
[np, nf] = check(np, nf, rate_wgs(p, 'K1', T, P, yb) < 0, 'K1 rate is negative when products exceed equilibrium');
[np, nf] = check(np, nf, rate_wgs(p, 'G2', T, P, yb) < 0, 'G2 rate is negative when products exceed equilibrium');

% 3. trends
[np, nf] = check(np, nf, rate_wgs(p, 'K1', 650, P, y) > r1, 'K1 rate rises with temperature');
[np, nf] = check(np, nf, rate_wgs(p, 'G2', 650, P, y) > r2, 'G2 rate rises with temperature');

% 4. G2 with no CO2 or H2 (negative orders) stays finite
y0 = [0.2 0.5 0 0 0.1 0.2];
r0 = rate_wgs(p, 'G2', T, P, y0);
[np, nf] = check(np, nf, isfinite(r0) && r0 > 0, 'G2 rate is finite and positive with zero CO2 and H2 (numerical floor)');

% 5. unavailable sets stop with the right error
ok = false; try, rate_wgs(p, 'G1', T, P, y); catch e, ok = strcmp(e.identifier, 'rate_wgs:G1'); end
[np, nf] = check(np, nf, ok, 'G1 stops with rate_wgs:G1');
ok = false; try, rate_wgs(p, 'G3', T, P, y); catch e, ok = strcmp(e.identifier, 'rate_wgs:G3'); end
[np, nf] = check(np, nf, ok, 'G3 stops with rate_wgs:G3');

% 6. validity flags
[~, i1] = rate_wgs(p, 'K1', 600, 5e5, [0.05 0.45 0.05 0.40 0.03 0.02]);
[np, nf] = check(np, nf, isempty(i1.flags), 'K1 at 600 K, 5 bar, 9 % CO dry: no flag');
yh = [0.16 0.60 0.06 0.14 0.02 0.02];
[~, i2] = rate_wgs(p, 'K1', T, P, yh);
[np, nf] = check(np, nf, any(contains(i2.flags, 'CO dry')), 'K1 at 40 % CO dry is flagged EXTRAPOLATED');

% 7. effectiveness factor, direct
[np, nf] = check(np, nf, abs(effectiveness_sphere(1e-3) - (1 - 6e-7)) < 1e-9, 'eta(1e-3) = 1 - 6e-7');
[np, nf] = check(np, nf, abs(effectiveness_sphere(1e-6) - 1) < 1e-10, 'eta -> 1 for small phi');
[np, nf] = check(np, nf, relerr(effectiveness_sphere(100), 0.0099666667) < 1e-6, 'eta(100) = 0.009967 (about 1/phi)');
ph = logspace(-2, 2, 50); ee = effectiveness_sphere(ph);
[np, nf] = check(np, nf, all(diff(ee) < 0), 'eta decreases monotonically with phi');
[np, nf] = check(np, nf, abs(effectiveness_sphere(1e-3 * (1 - 1e-9)) - effectiveness_sphere(1e-3 * (1 + 1e-9))) < 1e-8, 'eta is continuous at the series switch');

% 8. eta switch inside rate_wgs
o = struct('eta', 'on', 'd_p', 1e-3, 'D_eff', 1e-6, 'rho_p', 2100);
[ron, ion] = rate_wgs(p, 'K1', T, P, y, o);
[np, nf] = check(np, nf, abs(ion.eta - 0.9506) < 2e-3 && abs(ion.phi - 0.2976) < 2e-3, sprintf('K1 phi = 0.298, eta = 0.951 (got %.4f, %.4f)', ion.phi, ion.eta));
[np, nf] = check(np, nf, relerr(ron, r1 * ion.eta) < 1e-12, 'rate with eta on = intrinsic rate times eta');
o2 = o; o2.d_p = 5e-3;
[~, ion2] = rate_wgs(p, 'K1', T, P, y, o2);
[np, nf] = check(np, nf, ion2.eta < ion.eta, 'bigger particle gives smaller eta');
o3 = o; o3.D_eff = 1e-4;
[~, ion3] = rate_wgs(p, 'K1', T, P, y, o3);
[np, nf] = check(np, nf, ion3.eta > ion.eta && ion3.eta > 0.999, 'large D_eff gives eta close to 1');
ok = false; try, rate_wgs(p, 'K1', T, P, y, struct('eta', 'on')); catch e, ok = strcmp(e.identifier, 'rate_wgs:eta'); end
[np, nf] = check(np, nf, ok, 'eta = on without d_p and D_eff stops with rate_wgs:eta');

fprintf('--- %d passed, %d failed ---\n', np, nf);
if nf == 0, fprintf('ALL STEP 2B TESTS PASSED\n'); else, error('Step 2b tests failed'); end
end

function [np, nf] = check(np, nf, ok, msg)
if ok, np = np + 1; fprintf('  PASS  %s\n', msg);
else,  nf = nf + 1; fprintf('  FAIL  %s\n', msg); end
end

function e = relerr(a, b), e = abs(a - b) ./ abs(b); end
