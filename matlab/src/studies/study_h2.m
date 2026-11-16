function T = study_h2()
% H2 (CO inhibition): membrane x temperature x permeate pressure x inhibition set, 40 % CO base feed, K1, 9.5 wt % catalyst.
% rec_change is the relative change of H2 recovery against the same case without inhibition.
% The Melendez set has an unverified mapping and alpha > 1, so treat it as a stress case, not a prediction.
p = load_params();
sets = {'off', 'barbieri2008', 'melendez2017_60H2_15CO'};
rows = struct('membrane', {}, 'T_K', {}, 'P_perm_bar', {}, 'inhibition', {}, 'X_CO', {}, 'H2_recovery', {}, 'pinch_ratio', {}, 'rec_change', {});
for m = {'M7', 'M1'}
    for T0 = [623.15 673.15 723.15]
        for Pp = [101325 20265]
            for k = 1:numel(sets)
                ov = struct('T_in', T0, 'membrane', m{1}, 'P_perm', Pp, 'inhibition', sets{k});
                o = run_case(p, ov);
                rows(end + 1) = struct('membrane', m{1}, 'T_K', T0, 'P_perm_bar', Pp / 1e5, 'inhibition', sets{k}, ...
                                       'X_CO', o.summary.X_CO, 'H2_recovery', o.summary.H2_recovery, ...
                                       'pinch_ratio', o.profile.pinch_ratio(end), 'rec_change', NaN); %#ok<AGROW>
                fprintf('%s  %.0f K  %.2f bar  %-24s  X_CO=%.4f  recovery=%.4f\n', m{1}, T0, Pp / 1e5, sets{k}, rows(end).X_CO, rows(end).H2_recovery);
            end
        end
    end
end
T = struct2table(rows);
for i = 1:height(T)
    j = strcmp(T.membrane, T.membrane{i}) & T.T_K == T.T_K(i) & T.P_perm_bar == T.P_perm_bar(i) & strcmp(T.inhibition, 'off');
    T.rec_change(i) = T.H2_recovery(i) / T.H2_recovery(j) - 1;
end
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
writetable(T, fullfile(out, 'study_h2.csv'));
fprintf('Wrote %s\n', fullfile(out, 'study_h2.csv'));
end
