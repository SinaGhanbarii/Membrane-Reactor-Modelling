function run_tests_step1()
%RUN_TESTS_STEP1  Step 1 tests: thermodynamic core (MATLAB vs Python/Cantera).
%   Usage (from the matlab folder):   >> setup_paths
%                                     >> run_tests_step1
%   Reference values are in tests/vectors_step1.json (made by
%   scripts/export_matlab_reference.py).  Each check prints PASS or FAIL.

here = fileparts(mfilename('fullpath'));
vec  = jsondecode(fileread(fullfile(here, 'vectors_step1.json')));
c    = pbmr_const();
nPass = 0;  nFail = 0;

fprintf('\n=== pbmr1d MATLAB - Step 1 tests (thermodynamics) ===\n\n');

% ---------------------------------------------------------------- 1. species data
sp = load_species();
check('species order and data file', isequal(sp.names(:)', c.species) && numel(sp.T_mid) == 6, ...
      sprintf('%d species loaded', numel(sp.T_mid)));

% ---------------------------------------------------------------- 2. NASA properties vs Python
nT = numel(vec.T);
e_cp = 0;  e_h = 0;  e_g = 0;  e_K = 0;  e_Km = 0;  e_dH = 0;  e_xi = 0;
for k = 1:nT
    T = vec.T(k);
    [cp, h, ~, g0] = nasa_props(T);
    e_cp = max(e_cp, relerr(cp, vec.cp(k, :)'));
    e_h  = max(e_h,  relerr(h,  vec.h(k, :)'));
    e_g  = max(e_g,  relerr(g0, vec.g0(k, :)'));
    e_K  = max(e_K,  relerr(Keq_nasa(T), vec.K_nasa(k)));
    e_Km = max(e_Km, relerr(Keq_moe(T),  vec.K_moe(k)));
    e_dH = max(e_dH, relerr(dH_wgs(T),   vec.dH(k)));
    e_xi = max(e_xi, relerr(wgs_extent_eq(vec.F0, T), vec.xi(k)));
end
tol = 1e-9;
check('cp(T) matches Python/Cantera',            e_cp < tol, sprintf('worst rel. diff %.1e', e_cp));
check('h(T) matches Python/Cantera',             e_h  < tol, sprintf('worst rel. diff %.1e', e_h));
check('g0(T) matches Python/Cantera',            e_g  < tol, sprintf('worst rel. diff %.1e', e_g));
check('K_eq (NASA) matches Python/Cantera',      e_K  < tol, sprintf('worst rel. diff %.1e', e_K));
check('K_eq (Moe correlation) matches Python',   e_Km < 1e-12, sprintf('worst rel. diff %.1e', e_Km));
check('dH_wgs matches Python/Cantera',           e_dH < tol, sprintf('worst rel. diff %.1e', e_dH));
check('equilibrium extent matches Python',       e_xi < tol, sprintf('worst rel. diff %.1e', e_xi));

% ---------------------------------------------------------------- 3. physics checks
% NASA K_eq within 10 % of the Moe correlation, 250-500 C
worst = 0;
for TC = [250 300 350 400 450 500]
    T = TC + 273.15;
    worst = max(worst, abs(Keq_nasa(T) / Keq_moe(T) - 1));
end
check('K_eq (NASA) within 10 % of Moe, 250-500 C', worst < 0.10, sprintf('worst deviation %.1f %%', 100*worst));

% K falls with T and reaction is exothermic
Tg = linspace(500, 800, 31);
Kg = arrayfun(@Keq_nasa, Tg);
dHg = arrayfun(@dH_wgs, Tg);
check('K_eq decreases with T; dH < 0 (exothermic)', all(diff(Kg) < 0) && all(dHg < 0), ...
      sprintf('dH(650 K) = %.1f kJ/mol', dH_wgs(650)/1e3));

% van 't Hoff: d lnK / dT = dH / (R T^2)
T = 650;  dT = 1e-2;
dlnK = (log(Keq_nasa(T + dT)) - log(Keq_nasa(T - dT))) / (2*dT);
check('van''t Hoff consistency', relerr(dlnK, dH_wgs(T)/(c.R*T^2)) < 1e-5, ...
      sprintf('rel. diff %.1e', relerr(dlnK, dH_wgs(T)/(c.R*T^2))));

% element balances close and mass action holds at equilibrium
F0 = vec.F0(:);
eE = 0;  eQ = 0;
for TC = [250 400 500]
    T = TC + 273.15;
    F = wgs_equilibrium(F0, T);
    eE = max(eE, max(abs(element_flows(F) - element_flows(F0))) / max(element_flows(F0)));
    Q  = F(c.idx.CO2) * F(c.idx.H2) / (F(c.idx.CO) * F(c.idx.H2O));
    eQ = max(eQ, relerr(Q, Keq_nasa(T)));
end
check('element balances close at equilibrium', eE < 1e-12, sprintf('worst rel. imbalance %.1e', eE));
check('mass action Q = K at equilibrium',      eQ < 1e-9,  sprintf('worst rel. diff %.1e', eQ));

% large excess of steam gives (almost) complete CO conversion
Fs = vec.F0_steam(:);
Xs = co_conversion(Fs, wgs_equilibrium(Fs, 673.15));
check('excess steam gives X_CO > 0.9999', Xs > 0.9999 && relerr(Xs, vec.X_steam_673) < 1e-9, ...
      sprintf('X_CO = %.6f', Xs));

% ---------------------------------------------------------------- summary
fprintf('\n   T [C]    K_eq(NASA)   K_eq(Moe)   dH [kJ/mol]   X_CO,eq (test feed)\n');
for k = 1:nT
    T = vec.T(k);
    fprintf('  %6.1f   %9.3f   %9.3f   %10.2f   %8.4f\n', T - 273.15, Keq_nasa(T), Keq_moe(T), ...
            dH_wgs(T)/1e3, co_conversion(F0, wgs_equilibrium(F0, T)));
end

fprintf('\n=== %d passed, %d failed ===\n', nPass, nFail);
if nFail == 0
    fprintf('ALL STEP 1 TESTS PASSED\n\n');
else
    error('run_tests_step1:failed', '%d test(s) failed - paste this output to Claude.', nFail);
end

% ================================================================ nested helpers
    function check(name, ok, detail)
        if ok
            nPass = nPass + 1;  tag = 'PASS';
        else
            nFail = nFail + 1;  tag = 'FAIL';
        end
        fprintf('[%s] %-48s %s\n', tag, name, detail);
    end
end

function e = relerr(a, b)
%RELERR  Worst relative difference between a and b (true relative error).
e = max(abs(a(:) - b(:)) ./ max(abs(b(:)), realmin));
end
