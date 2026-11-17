function run_tests_step2()
%RUN_TESTS_STEP2  Checks for the parameter file, loader, validity flags and membrane flux.
%   Run:  >> setup_paths; run_tests_step2
fprintf('=== pbmr1d MATLAB - Step 2 tests (parameters, membrane flux) ===\n');
np = 0; nf = 0;
p = load_params();

% 1. every leaf has value, unit, source, kind
[bad, n] = audit(p, '');
[np, nf] = check(np, nf, isempty(bad), sprintf('all %d leaves have unit, source and kind  %s', n, strjoin(bad, ', ')));

% 2. missing value stops a run
ok = false; try, pget(p, 'membrane.M4.Q0'); catch e, ok = strcmp(e.identifier, 'pget:missing'); end
[np, nf] = check(np, nf, ok, 'pget stops on a null value (M4 pre-exponential)');
[np, nf] = check(np, nf, strcmp(p.rate.G2.pressure_unit.kind, 'reported') && strcmp(pget(p,'rate.G2.pressure_unit'), 'kPa'), 'G2 pressure unit is kPa and labelled reported (Hla 2009 Tables 4-6)');
[np, nf] = check(np, nf, strcmp(p.rate.K1.conc_basis_used.kind, 'inferred') && strcmp(pget(p,'rate.K1.conc_basis_used'), 'mol cm-3'), 'K1 concentration basis is mol cm-3 and labelled inferred');

% 3. permeance of each option at 400 C, unit driving force (review table 3C)
T = 673.15;
[np, nf] = check(np, nf, relerr(membrane_flux(p,'M7','off',T,1,0,0), 1.10e-3) < 0.015, 'M7 permeance at 400 C = 1.10e-3 (table 3C)');
[np, nf] = check(np, nf, relerr(membrane_flux(p,'M1','off',T,1,0,0), 2.80e-4) < 0.015, 'M1 permeance at 400 C = 2.80e-4 (table 3C)');
[np, nf] = check(np, nf, relerr(membrane_flux(p,'M8','off',T,1,0,0), 1.51e-4) < 0.015, 'M8 permeance at 400 C = 1.51e-4 (table 3C)');

% 4. flux scales with (pr^n - pp^n) and vanishes at zero driving force
J1 = membrane_flux(p,'M7','off',T,4e5,1e5,0);
J2 = membrane_flux(p,'M7','off',T,1e5,1e5,0);
[np, nf] = check(np, nf, J1 > 0 && J2 == 0, 'flux positive with driving force, zero without');

% 5. CO inhibition: Barbieri 2008 reduces the flux, and off equals no inhibition
Jb = membrane_flux(p,'M1','barbieri2008',647,4e5,1e5,100);
Jo = membrane_flux(p,'M1','off',647,4e5,1e5,100);
th = 0.023*100/(1+0.023*100);
[np, nf] = check(np, nf, relerr(Jb/Jo, 1-0.2319*th) < 1e-12, 'Barbieri (1 - alpha*theta) factor reproduced');

% 6. validity flags
r = p.rate.K1.validity;
[np, nf] = check(np, nf, isempty(check_validity(r, 600, 5, 0.05)), 'K1 at 600 K, 5 bar, 5 % CO: no flag');
w = check_validity(r, 600, 5, 0.40);
[np, nf] = check(np, nf, numel(w) == 1 && contains(w{1}, 'CO dry'), 'K1 at 40 % CO is flagged EXTRAPOLATED');
w = check_validity(r, 700, 4.5, 0.05);
[np, nf] = check(np, nf, numel(w) == 1 && contains(w{1}, 'T ='), 'K1 at 700 K is flagged EXTRAPOLATED');

% 7. internal consistency of transcribed numbers
[np, nf] = check(np, nf, relerr(pget(p,'rate.G1.k0'), exp(pget(p,'rate.K1.ln_k0'))) < 0.03, 'G1 k0 equals exp(ln k0 of K1) within 3 %');
c = pbmr_const();
[np, nf] = check(np, nf, relerr(pget(p,'rate.K1.E')/c.R, pget(p,'rate.G1.E_over_R')) < 1e-3, 'K1 E/R equals G1 E/R within 0.1 %');
cs = pget(p,'cases.hla2015.catalyst_mass_no_membrane'); is = pget(p,'cases.hla2015.inert_mass_no_membrane');
[np, nf] = check(np, nf, all(abs(cs./(cs+is) - 0.095) < 0.002), 'Hla Table 3: catalyst is 9.5 wt % in all three loadings');
sp = {'CO','CO2','H2','CH4','N2'}; fd = cellfun(@(s) pget(p, ['cases.target.test_feed_mmol_s.' s]), sp);
[np, nf] = check(np, nf, abs(fd(1)/sum(fd) - 0.40) < 1e-12, 'target test feed is 40 % CO dry');
[np, nf] = check(np, nf, isequal(size(p.rate.K1.validity.T_K), [2 1]), 'validity ranges load as [min; max]');

% 8. relations resolve and sweep bounds are ordered
lv = collect_leaves(p, '');
badrel = {}; nrel = 0; badsw = {};
for k = 1:numel(lv)
    L = lv{k}.leaf;
    if ~isempty(L.related)
        items = strsplit(L.related, ';');
        for q = 1:numel(items)
            ty = strsplit(strtrim(items{q}), ':'); tgt = strjoin(ty(2:end), ':'); nrel = nrel + 1;
            top = strsplit(tgt, '.'); 
            if isfield(p, top{1}) && ~has_path(p, tgt), badrel{end+1} = [lv{k}.id(2:end) ' -> ' tgt]; end %#ok<AGROW>
        end
    end
    if ~isempty(L.sweep) && ~(L.sweep(1) < L.sweep(2)), badsw{end+1} = lv{k}.id; end %#ok<AGROW>
end
[np, nf] = check(np, nf, isempty(badrel), sprintf('all %d relation links resolve  %s', nrel, strjoin(badrel, ', ')));
[np, nf] = check(np, nf, isempty(badsw), 'sweep bounds are ordered');
r = sweep_range(p, 'sweep.U_wall');
[np, nf] = check(np, nf, isequal(r, [0; 300]), 'sweep_range returns [0;300] for sweep.U_wall');
ok = false; try, sweep_range(p, 'rate.K1.E'); catch e, ok = strcmp(e.identifier, 'sweep_range:none'); end
[np, nf] = check(np, nf, ok, 'sweep_range stops where no bounds are recorded');
[np, nf] = check(np, nf, isempty(p.sweep.U_wall.value) && strcmp(p.sweep.U_wall.kind, 'sweep'), 'swept parameters have no single value and kind = sweep');

fprintf('--- %d passed, %d failed ---\n', np, nf);
if nf == 0, fprintf('ALL STEP 2 TESTS PASSED\n'); else, error('Step 2 tests failed'); end
end

function [np, nf] = check(np, nf, ok, msg)
if ok, np = np + 1; fprintf('  PASS  %s\n', msg);
else,  nf = nf + 1; fprintf('  FAIL  %s\n', msg); end
end

function e = relerr(a, b), e = abs(a - b) ./ abs(b); end

function [bad, n] = audit(s, path)
bad = {}; n = 0;
if isstruct(s) && isfield(s, 'value')
    n = 1;
    if ~all(isfield(s, {'unit','source','kind'})), bad = {path}; end
    return
end
if isstruct(s) && isscalar(s)
    fn = fieldnames(s);
    for k = 1:numel(fn)
        if any(strcmp(fn{k}, {'validity','x_meta','label','basis','status','form','scope_note', ...
                'description','built','keq_form','membrane_option','membrane_thick_option'})), continue; end
        [b, m] = audit(s.(fn{k}), [path '.' fn{k}]);
        bad = [bad, b]; n = n + m; %#ok<AGROW>
    end
end
end

function lv = collect_leaves(s, path)
lv = {};
if isstruct(s) && isfield(s, 'value')
    lv = {struct('id', path, 'leaf', s)}; return
end
if isstruct(s) && isscalar(s)
    fn = fieldnames(s);
    for k = 1:numel(fn)
        if any(strcmp(fn{k}, {'validity'})), continue; end
        sub = collect_leaves(s.(fn{k}), [path '.' fn{k}]);
        lv = [lv, sub]; %#ok<AGROW>
    end
end
end

function ok = has_path(p, path)
parts = strsplit(path, '.'); x = p; ok = true;
for k = 1:numel(parts)
    if ~isstruct(x) || ~isfield(x, parts{k}), ok = false; return; end
    x = x.(parts{k});
end
ok = isstruct(x) && isfield(x, 'value');
end
