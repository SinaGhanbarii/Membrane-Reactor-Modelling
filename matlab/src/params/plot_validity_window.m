function coverage = plot_validity_window(p)
%PLOT_VALIDITY_WINDOW  Data range of each literature parameter set versus the target window.
%   plot_validity_window(p)        three panels: temperature (C), CO dry fraction (%), pressure (bar)
%   coverage = plot_validity_window(p)  also returns a table: fraction of the target window each set covers
%   Bars are the ranges recorded in params_register.csv (kind = range).  A diamond is a single
%   data point.  Grey "n/r" means no range is recorded for that quantity (not that it is unlimited).
%   The shaded band is the target window from cases.target.
sets = {'rate.K1','K1 rate (Keiski)',1; 'rate.G1','G1 rate (Keiski constants, Barbieri)',1; ...
        'rate.G2','G2 rate (Hla)',1; 'rate.G3','G3 rate (Cu-ZnO, contrast)',1; ...
        'membrane.M7','M7 membrane (4.5 um, base)',2; 'membrane.M1','M1 membrane (60 um)',2; ...
        'membrane.M8','M8 membrane (100 um)',2; 'membrane.M4','M4 membrane (3.2 um)',2; ...
        'inhibition.barbieri2008','CO inhibition, Barbieri 2008',3; ...
        'inhibition.melendez2017_60H2_15CO','CO inhibition, Melendez (15 % CO)',3; ...
        'inhibition.melendez2017_60H2_5CO','CO inhibition, Melendez (5 % CO)',3};
grp = [0.12 0.47 0.71; 0.17 0.63 0.17; 1.00 0.50 0.05];

fields = {'T_K', 'yCO_dry', 'P_bar'};
ptitle = {'Temperature (C)', 'CO, dry basis (%)', 'Total pressure (bar)'};
scale  = [1 100 1];
offset = [-273.15 0 0];
tgt    = {pget(p,'cases.target.T_range_K'), pget(p,'cases.target.yCO_dry'), pget(p,'cases.target.P_ret_bar')};
xl     = {[150 520], [0 80], [0 20]};
n = size(sets, 1);
cv = nan(n, numel(fields));

figure('Name', 'Validity window', 'Color', 'w', 'Position', [60 60 1400 560]);
tl = tiledlayout(1, numel(fields), 'TileSpacing', 'compact');
for k = 1:numel(fields)
    nexttile; hold on; box on;
    tg = sort(tgt{k}(:)) * scale(k) + offset(k);
    xregion(tg(1), tg(2), 'FaceColor', [0.80 0.88 1.00], 'EdgeColor', 'none');
    xlim(xl{k}); ylim([0.4 n + 0.6]); set(gca, 'YDir', 'reverse');
    for i = 1:n
        v = get_validity(p, sets{i,1});
        col = grp(sets{i,3}, :);
        if ~isfield(v, fields{k})
            text(xl{k}(1) + 0.02*diff(xl{k}), i, 'n/r', 'Color', [0.55 0.55 0.55], 'FontSize', 8, 'VerticalAlignment', 'middle');
            continue
        end
        r = sort(v.(fields{k})(:)) * scale(k) + offset(k);
        if diff(r) < 1e-9
            plot(r(1), i, 'd', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        else
            rectangle('Position', [r(1) i-0.3 diff(r) 0.6], 'FaceColor', col, 'EdgeColor', 'none');
        end
        cv(i, k) = max(0, min(r(2), tg(2)) - max(r(1), tg(1))) / diff(tg);
    end
    xlabel(ptitle{k}); grid on;
    if k == 1
        yticks(1:n); yticklabels(sets(:,2));
    else
        yticks(1:n); yticklabels({});
    end
end
title(tl, 'Where each literature set has data (bars) versus the target window (shaded)');
coverage = array2table(cv, 'VariableNames', {'T_cover', 'yCO_cover', 'P_cover'}, 'RowNames', sets(:,1));
if nargout == 0
    disp('Fraction of the target window covered by each set (NaN = no range recorded):');
    disp(coverage);
    clear coverage
end
end

function v = get_validity(p, path)
v = struct();
parts = strsplit(path, '.');
s = p;
for k = 1:numel(parts)
    if ~isfield(s, parts{k}), return; end
    s = s.(parts{k});
end
if isfield(s, 'validity'), v = s.validity; end
end
