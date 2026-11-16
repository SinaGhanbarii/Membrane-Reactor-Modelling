function T = study_h5_h6()
% H5 and H6 sweep, isothermal wall: temperature x feed CO x rate law x catalyst loading x membrane on/off.
% Feed: 1 mmol/s dry gas, steam-to-CO 3, other dry species in the proportions of the register test feed.
% Writes results/study_h5_h6.csv at the repository root (results/ is ignored by Git).
p = load_params();
base = make_case(p);
others = base.F_in(3:6) / sum(base.F_in(3:6));    % CO2 H2 CH4 N2
dry = 1e-3;
temps = [573.15 623.15 673.15 723.15 743.15];
ycos = [0.25 0.40 0.50];
laws = {'K1', 'G2'};
loads = [0.095 1];

rows = struct('T_K', {}, 'yCO_dry', {}, 'law', {}, 'w_cat', {}, 'membrane', {}, ...
              'X_CO', {}, 'X_eq', {}, 'beta_out', {}, 'H2_recovery', {}, 'extrapolated', {}, 'error', {});
n = 0;
total = numel(temps) * numel(ycos) * numel(laws) * numel(loads) * 2;
for T0 = temps
    for y = ycos
        F = [y * dry; 3 * y * dry; (1 - y) * dry * others];
        for law = laws
            for w = loads
                for m = [1 0]
                    n = n + 1;
                    r = struct('T_K', T0, 'yCO_dry', y, 'law', law{1}, 'w_cat', w, 'membrane', m, ...
                               'X_CO', NaN, 'X_eq', NaN, 'beta_out', NaN, 'H2_recovery', NaN, 'extrapolated', NaN, 'error', "");
                    try
                        s = run_case(p, struct('F_in', F, 'T_in', T0, 'law', law{1}, 'w_cat', w, 'perm_scale', m)).summary;
                        r.X_CO = s.X_CO;  r.X_eq = s.X_eq;  r.beta_out = s.beta_out;
                        r.H2_recovery = s.H2_recovery;  r.extrapolated = ~isempty(s.flags);
                    catch err
                        r.error = string(err.message);
                    end
                    rows(end + 1) = r; %#ok<AGROW>
                    fprintf('%3d/%d  T=%.0f K  yCO=%.2f  %s  w=%.3f  membrane=%d  X_CO=%.4f\n', n, total, T0, y, law{1}, w, m, r.X_CO);
                end
            end
        end
    end
end

T = struct2table(rows);
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
writetable(T, fullfile(out, 'study_h5_h6.csv'));
fprintf('Wrote %s\n', fullfile(out, 'study_h5_h6.csv'));
end
