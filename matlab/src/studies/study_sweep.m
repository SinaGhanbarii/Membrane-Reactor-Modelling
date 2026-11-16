function T = study_sweep()
% Sweep gas and flow direction: membrane x temperature x sweep flow x co/counter-current, 40 % CO base feed, K1.
% Writes results/study_sweep.csv. permeate_H2_fraction is the H2 share of the permeate stream (H2 + sweep).
p = load_params();
sweeps = [0 0.5 1 2 5] * 1e-3;
rows = struct('membrane', {}, 'T_K', {}, 'sweep_mmol_s', {}, 'flow', {}, 'X_CO', {}, 'H2_recovery', {}, ...
              'permeate_H2_fraction', {}, 'pinch_fraction', {}, 'atom_error', {});
for m = {'M7', 'M1'}
    for T0 = [623.15 673.15 723.15]
        for sw = sweeps
            for fl = {'cocurrent', 'counter'}
                if sw == 0 && strcmp(fl{1}, 'counter'), continue; end
                o = run_case(p, struct('membrane', m{1}, 'T_in', T0, 'sweep', sw, 'flow', fl{1}));
                sm = o.summary;
                rows(end + 1) = struct('membrane', m{1}, 'T_K', T0, 'sweep_mmol_s', sw * 1e3, 'flow', fl{1}, ...
                                       'X_CO', sm.X_CO, 'H2_recovery', sm.H2_recovery, ...
                                       'permeate_H2_fraction', sm.permeate_H2_fraction, ...
                                       'pinch_fraction', sm.regime_fraction(3), 'atom_error', sm.atom_error); %#ok<AGROW>
                fprintf('%s  %.0f K  sweep %.1f mmol/s  %-9s  X_CO=%.4f  recovery=%.4f  H2 in permeate=%.3f\n', ...
                        m{1}, T0, sw * 1e3, fl{1}, sm.X_CO, sm.H2_recovery, sm.permeate_H2_fraction);
            end
        end
    end
end
T = struct2table(rows);
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
writetable(T, fullfile(out, 'study_sweep.csv'));
end
