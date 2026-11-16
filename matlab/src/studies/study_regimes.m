function T = study_regimes()
% Regime study with the K1/G2 rate laws, 9.5 wt % catalyst unless stated, feed H2 held fixed when CO varies.
%   A: membrane (M7, M1) x bed length x permeate pressure x temperature
%   B: pore diffusion (eta on) x d_p x D_eff x temperature x catalyst loading
%   C: wall mode (isothermal, adiabatic, cooled) x temperature x catalyst loading
%   D: feed CO x temperature x rate law x membrane on/off
% Regime label (a proposal, thresholds to be agreed): reaction if beta_out < 0.9, pinch if the outlet H2
% partial pressure is within 5 % of the permeate pressure, otherwise permeance.
% Writes results/study_regimes.csv at the repository root (results/ is ignored by Git).
p = load_params();
b = make_case(p);
oth = b.F_in(3:6);
feed = @(y) [sum(oth) * y / (1 - y); 3 * sum(oth) * y / (1 - y); oth];

cases = struct('part', {}, 'ov', {});
for T0 = [623.15 673.15 723.15]                                             % A
    for mem = {'M7', 'M1'}
        for L = [0.27 0.54 1.08]
            for Pp = [20265 50663 101325]
                cases(end + 1) = struct('part', 'A', 'ov', struct('T_in', T0, 'membrane', mem{1}, 'L', L, 'P_perm', Pp)); %#ok<AGROW>
            end
        end
    end
end
for T0 = [573.15 623.15 673.15]                                             % B
    for w = [0.095 1]
        cases(end + 1) = struct('part', 'B', 'ov', struct('T_in', T0, 'w_cat', w)); %#ok<AGROW>
        for dp = [0.5e-3 0.7e-3 1.5e-3 3.2e-3]
            for De = [1e-7 1e-6 1e-5]
                cases(end + 1) = struct('part', 'B', 'ov', struct('T_in', T0, 'w_cat', w, 'eta', 'on', 'd_p', dp, 'D_eff', De)); %#ok<AGROW>
            end
        end
    end
end
for T0 = [623.15 673.15]                                                    % C
    for w = [0.095 1]
        cases(end + 1) = struct('part', 'C', 'ov', struct('T_in', T0, 'w_cat', w)); %#ok<AGROW>
        cases(end + 1) = struct('part', 'C', 'ov', struct('T_in', T0, 'w_cat', w, 'wall', 'adiabatic')); %#ok<AGROW>
        for U = [10 50 300]
            cases(end + 1) = struct('part', 'C', 'ov', struct('T_in', T0, 'w_cat', w, 'wall', 'cooled', 'U', U, 'T_wall', T0)); %#ok<AGROW>
        end
    end
end
for y = [0.25 0.40 0.50]                                                    % D
    for T0 = [623.15 673.15 723.15]
        for law = {'K1', 'G2'}
            for m = [1 0]
                cases(end + 1) = struct('part', 'D', 'ov', struct('F_in', feed(y), 'T_in', T0, 'law', law{1}, 'perm_scale', m)); %#ok<AGROW>
            end
        end
    end
end

nan = NaN;
tmpl = struct('part', '', 'T_K', nan, 'yCO_dry', nan, 'law', '', 'w_cat', nan, 'membrane', '', 'L', nan, 'P_perm_bar', nan, ...
              'wall', '', 'U', nan, 'd_p', nan, 'D_eff', nan, 'X_CO', nan, 'X_eq', nan, 'beta_out', nan, 'H2_recovery', nan, ...
              'pinch_ratio', nan, 'T_max', nan, 'T_out', nan, 'eta_min', nan, 'regime', '', 'error', "");
rows = tmpl([]);
for i = 1:numel(cases)
    r = tmpl;  r.part = cases(i).part;
    try
        c = make_case(p, cases(i).ov);
        out = run_case(p, cases(i).ov);
        s = out.summary;  Fo = out.sol.F(end, :);
        r.T_K = c.T_in;  r.yCO_dry = c.F_in(1) / (sum(c.F_in) - c.F_in(2));  r.law = c.law;  r.w_cat = c.w_cat;
        r.membrane = c.membrane;  if c.perm_scale == 0, r.membrane = 'none'; end
        r.L = c.L;  r.P_perm_bar = c.P_perm / 1e5;  r.wall = c.wall;
        if strcmp(c.wall, 'cooled'), r.U = c.U; end
        if strcmp(c.eta, 'on'), r.d_p = c.d_p;  r.D_eff = c.D_eff; end
        r.X_CO = s.X_CO;  r.X_eq = s.X_eq;  r.beta_out = s.beta_out;  r.H2_recovery = s.H2_recovery;
        r.pinch_ratio = Fo(4) / sum(Fo) * out.sol.P(end) / c.P_perm;
        r.T_max = max(out.sol.T);  r.T_out = s.T_out;  r.eta_min = min(out.profile.eta);
        if c.perm_scale == 0,         r.regime = 'no membrane';
        elseif s.H2_recovery < 0.01,  r.regime = 'membrane inactive';
        elseif s.beta_out < 0.9,      r.regime = 'reaction';
        elseif r.pinch_ratio < 1.05,  r.regime = 'pinch';
        else,                         r.regime = 'permeance';
        end
    catch err
        r.error = string(err.message);
    end
    rows(end + 1) = r; %#ok<AGROW>
    fprintf('%3d/%d  %s  T=%.0f  X_CO=%.4f  %s\n', i, numel(cases), r.part, r.T_K, r.X_CO, r.regime);
end

T = struct2table(rows);
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
writetable(T, fullfile(out, 'study_regimes.csv'));
fprintf('Wrote %s\n', fullfile(out, 'study_regimes.csv'));
end
