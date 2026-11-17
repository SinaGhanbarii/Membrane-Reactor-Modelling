function run_verification()
%RUN_VERIFICATION  Runs the reactor model over many different conditions and checks each result.
%   Run:  >> clear functions; setup_paths; run_verification
%   A  20 cases against the independent Python model (Cantera thermo, SciPy BDF, same register values)
%   B  physical rules on 54 runs: 3 membranes x 3 temperatures x isothermal/adiabatic x 3 permeate-side modes
%   C  trends that must go one way (pressure, bed length, permeance, sweep, catalyst, feed)
%   D  unusual but legal conditions: they must run and close the atom balance
%   E  bad inputs: each must stop with the right error
%   Writes results/verification_report.csv (one line per check). Nothing here is a validation against
%   measurements; it shows that the code is self-consistent and agrees with a second implementation.
t0 = tic;
fprintf('=== pbmr1d MATLAB - verification over varied conditions ===\n');
p = load_params();
v = struct('group', {{}}, 'name', {{}}, 'ok', logical([]), 'detail', {{}});

% ---------- A. against the independent Python model ----------
fprintf('\n-- A. agreement with the Python reference --\n');
ref = {
    'base co-current, isothermal', {}, 0.942717, 0.391191, 673.150; ...
    'base, adiabatic', {'wall', 'adiabatic'}, 0.829403, 0.319087, 823.804; ...
    'T_in 623 K, 0.54 m', {'T_in', 623.15}, 0.939372, 0.357384, 623.150; ...
    'T_in 723 K', {'T_in', 723.15}, 0.912418, 0.374373, 723.150; ...
    'short bed 0.27 m', {'L', 0.27}, 0.932428, 0.318750, 673.150; ...
    'long bed 1.08 m', {'L', 1.08}, 0.943823, 0.399432, 673.150; ...
    'permeate 0.2 bar', {'P_perm', 0.2e5}, 0.990591, 0.904819, 673.150; ...
    'permeate 0.5 bar', {'P_perm', 0.5e5}, 0.974986, 0.742865, 673.150; ...
    'permeate 2 bar', {'P_perm', 2e5}, 0.910134, 0.000000, 673.150; ...
    'sweep 1 mmol/s co-current', {'sweep', 1e-3, 'backperm', 'off'}, 0.981360, 0.809386, 673.150; ...
    'sweep 1 mmol/s counter', {'sweep', 1e-3, 'flow', 'counter', 'backperm', 'off'}, 0.999973, 0.999986, 673.150; ...
    'sweep 0.2 mmol/s counter', {'sweep', 2e-4, 'flow', 'counter', 'backperm', 'off'}, 0.986529, 0.922621, 673.150; ...
    'sweep 5 mmol/s counter', {'sweep', 5e-3, 'flow', 'counter', 'backperm', 'off'}, 0.999991, 0.999995, 673.150; ...
    'sweep 1 mmol/s counter, adiabatic', {'sweep', 1e-3, 'flow', 'counter', 'wall', 'adiabatic', 'backperm', 'off'}, 1.000000, 1.000000, 860.107; ...
    'membrane at half permeance', {'perm_scale', 0.5}, 0.937471, 0.339134, 673.150; ...
    'membrane at twice permeance', {'perm_scale', 2}, 0.943769, 0.399246, 673.150; ...
    'catalyst 100 wt % (w_cat 1)', {'w_cat', 1}, 0.943424, 0.394998, 673.150; ...
    'no membrane, adiabatic, long bed', {'perm_scale', 0, 'wall', 'adiabatic', 'L', 1.08}, 0.776605, 0.000000, 814.310; ...
    'P_in 10 bar', {'P_in', 10e5}, 0.977488, 0.768656, 673.150; ...
    'T_in 673 K, P_perm 0.2 bar, 1.08 m, counter sweep 2 mmol/s', {'T_in', 673.15, 'P_perm', 0.2e5, 'L', 1.08, 'sweep', 2e-3, 'flow', 'counter', 'backperm', 'off'}, 1.000000, 1.000000, 673.150; ...
    };
for i = 1:size(ref, 1)
    try
        s = run_case(p, struct(ref{i, 2}{:})).summary;
        tol = 3e-5;  if ref{i, 4} > 0.99, tol = 2e-4; end
        ok = abs(s.X_CO - ref{i, 3}) < tol && abs(s.H2_recovery - ref{i, 4}) < tol && abs(s.T_out - ref{i, 5}) < 0.05;
        v = add(v, 'A', ref{i, 1}, ok, sprintf('X_CO %.6f (ref %.6f), recovery %.6f (ref %.6f), T_out %.2f (ref %.2f)', ...
                s.X_CO, ref{i, 3}, s.H2_recovery, ref{i, 4}, s.T_out, ref{i, 5}));
    catch err
        v = add(v, 'A', ref{i, 1}, false, ['error: ' err.message]);
    end
end

% ---------- B. physical rules over a grid ----------
fprintf('\n-- B. physical rules on a grid of conditions --\n');
membranes = {'M7', 'M1', 'M8'};
temps = [623.15 673.15 723.15];
walls = {'isothermal', 'adiabatic'};
modes = {'no sweep', 'co-current sweep', 'counter-current sweep'};
for m = membranes
    for T0 = temps
        for w = walls
            rec_mode = zeros(1, 3);
            for im = 1:3
                ov = struct('membrane', m{1}, 'T_in', T0, 'wall', w{1});
                if im == 2, ov.sweep = 1e-3; end
                if im == 3, ov.sweep = 1e-3;  ov.flow = 'counter'; end
                lab = sprintf('%s, %.0f K, %s, %s', m{1}, T0, w{1}, modes{im});
                try
                    o = run_case(p, ov);
                    sm = o.summary;  s = o.sol;  c = o.case;
                    rec_mode(im) = sm.H2_recovery;
                    ok = sm.atom_error < 1e-8;
                    v = add(v, 'B', [lab ': atom balance closes'], ok, sprintf('relative error %.1e', sm.atom_error));
                    ok = sm.X_CO >= 0 && sm.X_CO <= 1 + 1e-9 && sm.H2_recovery >= 0 && sm.H2_recovery <= 1 + 1e-6;
                    v = add(v, 'B', [lab ': X_CO and recovery between 0 and 1'], ok, sprintf('X_CO %.5f, recovery %.5f', sm.X_CO, sm.H2_recovery));
                    ok = min(s.F(:)) > -1e-10 && min(s.P) > 0;
                    v = add(v, 'B', [lab ': no negative flows or pressure'], ok, sprintf('min flow %.1e mol/s', min(s.F(:))));
                    if im == 3, dG = diff(s.Fp); else, dG = -diff(s.Fp); end
                    ok = max(dG) < 1e-9;
                    v = add(v, 'B', [lab ': permeate H2 builds up in the flow direction of the permeate stream'], ok, sprintf('largest step against it %.1e mol/s', max(dG)));
                    if strcmp(w{1}, 'isothermal')
                        ok = max(abs(s.T - c.T_in)) < 1e-9;
                    else
                        ok = sm.T_out >= c.T_in - 1e-6 && all(isfinite(s.T));
                    end
                    v = add(v, 'B', [lab ': temperature behaves for the wall mode'], ok, sprintf('T_out %.2f K', sm.T_out));
                catch err
                    v = add(v, 'B', lab, false, ['error: ' err.message]);
                end
            end
            ok = rec_mode(3) >= rec_mode(2) - 1e-9 && rec_mode(2) >= rec_mode(1) - 1e-9;
            v = add(v, 'B', sprintf('%s, %.0f K, %s: recovery no sweep <= co-current sweep <= counter-current sweep', m{1}, T0, w{1}), ok, ...
                    sprintf('%.4f <= %.4f <= %.4f', rec_mode));
            if strcmp(w{1}, 'isothermal')
                nm = run_case(p, struct('membrane', m{1}, 'T_in', T0, 'perm_scale', 0)).summary;
                mm = run_case(p, struct('membrane', m{1}, 'T_in', T0)).summary;
                v = add(v, 'B', sprintf('%s, %.0f K, isothermal: the membrane never lowers CO conversion', m{1}, T0), mm.X_CO >= nm.X_CO - 1e-9, ...
                        sprintf('%.5f with membrane, %.5f without', mm.X_CO, nm.X_CO));
            end
        end
    end
end

% ---------- C. trends ----------
fprintf('\n-- C. trends --\n');
Pp = [0.1 0.2 0.5 1 2 3] * 1e5;
M = scan(p, 'P_perm', Pp, struct(), {'H2_recovery', 'X_CO'});
v = add(v, 'C', 'recovery falls as permeate pressure rises (0.1 to 3 bar)', all(diff(M(:, 1)) <= 1e-9), sprintf('%.3f ', M(:, 1)));
v = add(v, 'C', 'CO conversion falls as permeate pressure rises', all(diff(M(:, 2)) <= 1e-9), sprintf('%.4f ', M(:, 2)));
M = scan(p, 'L', [0.05 0.1 0.27 0.54 1.08 2.16], struct(), {'X_CO', 'H2_recovery'});
v = add(v, 'C', 'CO conversion does not fall with bed length', all(diff(M(:, 1)) >= -1e-9), sprintf('%.4f ', M(:, 1)));
v = add(v, 'C', 'H2 recovery does not fall with bed length', all(diff(M(:, 2)) >= -1e-9), sprintf('%.4f ', M(:, 2)));
M = scan(p, 'perm_scale', [0 0.25 0.5 1 2 4], struct(), {'H2_recovery'});
v = add(v, 'C', 'H2 recovery does not fall with membrane permeance (0 to 4 times)', all(diff(M) >= -1e-9), sprintf('%.4f ', M));
M = scan(p, 'sweep', [0 0.2 0.5 1 2 5] * 1e-3, struct('flow', 'counter'), {'H2_recovery', 'permeate_H2_fraction'});
v = add(v, 'C', 'counter-current recovery does not fall with sweep flow', all(diff(M(:, 1)) >= -1e-9), sprintf('%.4f ', M(:, 1)));
v = add(v, 'C', 'H2 share of the permeate falls as sweep rises', all(diff(M(:, 2)) <= 1e-9), sprintf('%.3f ', M(:, 2)));
M = scan(p, 'P_in', [3 4.5 7 10 15] * 1e5, struct(), {'H2_recovery'});
v = add(v, 'C', 'H2 recovery rises with retentate pressure', all(diff(M) >= -1e-9), sprintf('%.4f ', M));
M = scan(p, 'w_cat', [0.02 0.095 0.3 1], struct('L', 0.1, 'perm_scale', 0), {'X_CO'});
v = add(v, 'C', 'CO conversion rises with catalyst fraction (short bed, no membrane)', all(diff(M) >= -1e-9), sprintf('%.4f ', M));
M = scan(p, 'T_in', [573.15 673.15], struct('L', 0.1, 'perm_scale', 0), {'X_CO'});
v = add(v, 'C', 'a short bed converts more CO at 673 K than at 573 K (kinetics)', M(2) > M(1), sprintf('%.4f at 573 K, %.4f at 673 K', M(1), M(2)));
o = run_case(p, struct('P_perm', 2e5)).summary;
v = add(v, 'C', 'permeate at 2 bar (above the retentate H2 pressure): no permeation, equilibrium conversion', o.H2_recovery == 0 && abs(o.X_CO - o.X_eq) < 1e-3, sprintf('recovery %.2g, X_CO %.4f, X_eq %.4f', o.H2_recovery, o.X_CO, o.X_eq));
a = run_case(p, struct('wall', 'adiabatic', 'perm_scale', 0, 'L', 2)).summary;
i = run_case(p, struct('perm_scale', 0, 'L', 2)).summary;
v = add(v, 'C', 'adiabatic bed ends hotter and below the isothermal conversion (exothermic, equilibrium shifts back)', a.T_out > 673.15 && a.X_CO < i.X_CO, sprintf('X_CO %.4f adiabatic, %.4f isothermal, T_out %.1f K', a.X_CO, i.X_CO, a.T_out));
c0 = run_case(p, struct('wall', 'cooled', 'U', 1e4, 'T_wall', 673.15, 'perm_scale', 0));
v = add(v, 'C', 'a very strongly cooled wall holds the bed at the wall temperature', max(abs(c0.sol.T - 673.15)) < 1.0, sprintf('largest deviation %.2f K', max(abs(c0.sol.T - 673.15))));
c1 = run_case(p, struct('wall', 'cooled', 'U', 50, 'T_wall', 673.15, 'perm_scale', 0));
c2 = run_case(p, struct('wall', 'cooled', 'U', 300, 'T_wall', 673.15, 'perm_scale', 0));
v = add(v, 'C', 'peak temperature falls as wall heat transfer rises (adiabatic >= U 50 >= U 300)', a.T_out >= max(c1.sol.T) && max(c1.sol.T) >= max(c2.sol.T), sprintf('%.1f >= %.1f >= %.1f K', a.T_out, max(c1.sol.T), max(c2.sol.T)));
d0 = run_case(p).summary;
d1 = run_case(p, struct('dp', 'on', 'mu', 2.5e-5)).summary;
v = add(v, 'C', 'pressure drop on: outlet pressure below inlet, recovery not above the no-drop case', d1.P_out < 4.5e5 && d1.P_out > 0 && d1.H2_recovery <= d0.H2_recovery + 1e-9, sprintf('P_out %.4f bar, recovery %.4f vs %.4f', d1.P_out / 1e5, d1.H2_recovery, d0.H2_recovery));
e1 = run_case(p, struct('eta', 'on', 'D_eff', 1e-6)).summary;
v = add(v, 'C', 'pore diffusion limits (eta on) cannot raise CO conversion', e1.X_CO <= d0.X_CO + 1e-9, sprintf('%.4f with eta, %.4f without', e1.X_CO, d0.X_CO));
g1 = run_case(p, struct('law', 'G2')).summary;
v = add(v, 'C', 'G2 rate law runs and keeps X_CO between 0 and 1', g1.X_CO >= 0 && g1.X_CO <= 1, sprintf('X_CO %.4f', g1.X_CO));

% ---------- D. unusual but legal conditions ----------
fprintf('\n-- D. unusual conditions --\n');
base = make_case(p);
dry = {
    'almost steam-free feed (steam-to-CO 0.01)', struct('F_in', 1e-3 * [1; 0.01; 0; 0; 0; 0])
    'very CO-rich feed (80 % CO dry), steam-to-CO 3', struct('F_in', 1e-3 * [0.8; 2.4; 0.05; 0.05; 0.05; 0.05])
    'diluted feed (5 % CO)',                struct('F_in', 1e-3 * [0.05; 0.15; 0.2; 0.2; 0.1; 0.3])
    'pure CO and steam, no products',       struct('F_in', 1e-3 * [1; 3; 0; 0; 0; 0])
    'inlet already at WGS equilibrium composition', struct('F_in', wgs_equilibrium(base.F_in, 673.15))
    'extremely short bed (1 mm)',           struct('L', 1e-3)
    'very long bed (5 m)',                  struct('L', 5)
    'low temperature (500 K)',              struct('T_in', 500)
    'high temperature (800 K)',             struct('T_in', 800)
    'high pressure (30 bar)',               struct('P_in', 30e5)
    'low retentate pressure (1.5 bar), permeate 0.1 bar', struct('P_in', 1.5e5, 'P_perm', 0.1e5)
    'no catalyst (w_cat = 0), membrane on', struct('w_cat', 0)
    'counter-current with a very small sweep (1e-6 mol/s)', struct('sweep', 1e-6, 'flow', 'counter')
    'counter-current with a very large sweep (1 mol/s)',   struct('sweep', 1, 'flow', 'counter')
    'counter-current, adiabatic, G2 law',   struct('sweep', 1e-3, 'flow', 'counter', 'wall', 'adiabatic', 'law', 'G2')
    'counter-current, cooled wall, eta on', struct('sweep', 1e-3, 'flow', 'counter', 'wall', 'cooled', 'U', 50, 'T_wall', 673.15, 'eta', 'on', 'D_eff', 5e-6)
    'counter-current with M1 at 573 K',     struct('sweep', 1e-3, 'flow', 'counter', 'membrane', 'M1', 'T_in', 573.15)
    'counter-current, pressure drop on',    struct('sweep', 1e-3, 'flow', 'counter', 'dp', 'on', 'mu', 2.5e-5)
    'Barbieri inhibition on',               struct('inhibition', 'barbieri2008')
    'Melendez inhibition (alpha > 1) on',   struct('inhibition', 'melendez2017_60H2_15CO')
    'permeate pressure equal to retentate, with sweep and back-permeation', struct('P_perm', 4.5e5, 'sweep', 1e-3, 'flow', 'counter')
    };
warning('off', 'make_case:pressure');
for i = 1:size(dry, 1)
    try
        o = run_case(p, dry{i, 2});  sm = o.summary;
        ok = sm.atom_error < 1e-6 && sm.X_CO <= 1 + 1e-9 && min(o.sol.F(:)) > -1e-10 && all(isfinite([sm.X_CO sm.H2_recovery sm.T_out]));
        v = add(v, 'D', dry{i, 1}, ok, sprintf('X_CO %.4f, recovery %.4f, T_out %.1f K, atom error %.1e', sm.X_CO, sm.H2_recovery, sm.T_out, sm.atom_error));
    catch err
        v = add(v, 'D', dry{i, 1}, false, ['error: ' err.message]);
    end
end
warning('on', 'make_case:pressure');

% ---------- E. bad inputs ----------
fprintf('\n-- E. bad inputs must stop with the right error --\n');
bad = {
    'unknown field',                    struct('foo', 1),                              'make_case:field'
    'wall = cooled without U',          struct('wall', 'cooled', 'T_wall', 673),       'make_case:wall'
    'wall = cooled with negative U',    struct('wall', 'cooled', 'U', -5, 'T_wall', 673), 'make_case:wall'
    'unknown wall',                     struct('wall', 'hot'),                         'make_case:wall'
    'eta on without D_eff',             struct('eta', 'on'),                           'make_case:eta'
    'eta on with D_eff = 0',            struct('eta', 'on', 'D_eff', 0),               'make_case:eta'
    'eta mistyped ("yes")',             struct('eta', 'yes'),                          'make_case:switch'
    'dp mistyped ("true")',             struct('dp', 'true'),                          'make_case:switch'
    'dp on without mu',                 struct('dp', 'on'),                            'make_case:mu'
    'dp on with negative mu',           struct('dp', 'on', 'mu', -1),                  'make_case:mu'
    'unknown flow',                     struct('flow', 'cross'),                       'make_case:flow'
    'negative sweep',                   struct('sweep', -1e-3),                        'make_case:sweep'
    'NaN sweep',                        struct('sweep', NaN),                          'make_case:sweep'
    'vector sweep',                     struct('sweep', [1 2] * 1e-3),                 'make_case:sweep'
    'backperm mistyped',                struct('backperm', 'maybe'),                   'make_case:backperm'
    'inlet temperature in C',           struct('T_in', 25),                            'make_case:range'
    'inlet temperature NaN',            struct('T_in', NaN),                           'make_case:range'
    'pressure in bar',                  struct('P_in', 4.5),                           'make_case:range'
    'bed length in mm',                 struct('L', 540),                              'make_case:range'
    'negative voidage',                 struct('eps', -0.1),                           'make_case:range'
    'catalyst fraction above 1',        struct('w_cat', 1.5),                          'make_case:range'
    'negative permeance scale',         struct('perm_scale', -1),                      'make_case:range'
    'shell smaller than membrane',      struct('D_s', 0.005),                          'make_case:geometry'
    'feed with 5 entries',              struct('F_in', ones(5, 1)),                    'make_case:feed'
    'feed with negative flow',          struct('F_in', [1; 3; -1; 0; 0; 0] * 1e-3),    'make_case:feed'
    'feed without CO',                  struct('F_in', [0; 3; 1; 0; 0; 0] * 1e-3),     'make_case:feed'
    'unknown rate law',                 struct('law', 'K9'),                           'make_case:law'
    'unknown membrane',                 struct('membrane', 'M99'),                     'make_case:membrane'
    'unknown inhibition set',           struct('inhibition', 'nope'),                  'make_case:inhibition'
    'pressure drop that empties the bed', struct('dp', 'on', 'mu', 1e4),               'solve_cocurrent:stopped'
    };
for i = 1:size(bad, 1)
    got = 'no error';
    try
        run_case(p, bad{i, 2});
    catch err
        got = err.identifier;
    end
    v = add(v, 'E', bad{i, 1}, strcmp(got, bad{i, 3}), sprintf('expected %s, got %s', bad{i, 3}, got));
end

% ---------- summary ----------
n = numel(v.ok);  nf = sum(~v.ok);
fprintf('\n=== %d checks, %d passed, %d failed, %.0f s ===\n', n, n - nf, nf, toc(t0));
for g = 'ABCDE'
    j = strcmp(v.group, g);
    fprintf('  group %s: %d of %d passed\n', g, sum(v.ok(j)), sum(j));
end
if nf > 0
    fprintf('Failed checks:\n');
    for j = find(~v.ok)
        fprintf('  [%s] %s -> %s\n', v.group{j}, v.name{j}, v.detail{j});
    end
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
T = table(v.group(:), v.name(:), v.ok(:), v.detail(:), 'VariableNames', {'group', 'check', 'passed', 'detail'});
writetable(T, fullfile(out, 'verification_report.csv'));
fprintf('Wrote %s\n', fullfile(out, 'verification_report.csv'));
if nf == 0, fprintf('ALL VERIFICATION CHECKS PASSED\n'); else, error('Verification failed: %d checks', nf); end
end

function v = add(v, group, name, ok, detail)
v.group{end + 1} = group;  v.name{end + 1} = name;  v.ok(end + 1) = logical(ok);  v.detail{end + 1} = detail;
if ok, fprintf('  PASS  %s\n', name); else, fprintf('  FAIL  %s\n        %s\n', name, detail); end
end

function M = scan(p, key, vals, ov, fields)
% One run per value of ov.(key); columns follow fields.
M = zeros(numel(vals), numel(fields));
for i = 1:numel(vals)
    o = ov;  o.(key) = vals(i);
    sm = run_case(p, o).summary;
    for j = 1:numel(fields), M(i, j) = sm.(fields{j}); end
end
end
