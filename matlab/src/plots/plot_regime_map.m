function plot_regime_map()
% Figures from results/study_regimes.csv: regime maps (A), rate-law comparison (D), hot spots (C).
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
d = readtable(fullfile(root, 'results', 'study_regimes.csv'), 'TextType', 'string');

Ls = [0.27 0.54 1.08];  Ps = [0.20265 0.50663 1.01325];  Ts = [623.15 673.15 723.15];
A = d(d.part == "A", :);
figure('Color', 'w');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact');
title(tl, 'Outlet regime and H2 recovery, 40 % CO, 9.5 wt % catalyst (rxn = reaction, perm = permeance)');
xlabel(tl, 'bed length (m)');  ylabel(tl, 'permeate pressure (bar)');
for m = ["M7" "M1"]
    for T0 = Ts
        nexttile; hold on;
        rows = A(A.membrane == m & near(A.T_K, T0), :);
        for i = 1:height(rows)
            [~, ix] = min(abs(Ls - rows.L(i)));  [~, iy] = min(abs(Ps - rows.P_perm_bar(i)));
            rectangle('Position', [ix - 0.5, iy - 0.5, 1, 1], 'FaceColor', regime_color(rows.regime(i)), 'EdgeColor', 'w');
            text(ix, iy, sprintf('%s\n%.2f', regime_tag(rows.regime(i)), rows.H2_recovery(i)), 'HorizontalAlignment', 'center');
        end
        axis([0.5 3.5 0.5 3.5]);  axis equal tight;
        xticks(1:3); xticklabels(string(Ls));  yticks(1:3); yticklabels(string(round(Ps, 2)));
        title(sprintf('%s, %.0f C', m, T0 - 273.15));
    end
end

D = d(d.part == "D" & d.membrane ~= "none", :);
figure('Color', 'w');
labels = strings(0, 1);  Xk = [];  Xg = [];
for y = [0.25 0.40 0.50]
    for T0 = Ts
        k = D(near(D.yCO_dry, y) & near(D.T_K, T0) & D.law == "K1", :);
        g = D(near(D.yCO_dry, y) & near(D.T_K, T0) & D.law == "G2", :);
        labels(end + 1, 1) = sprintf('%.0f C, %.0f %% CO', T0 - 273.15, 100 * y); %#ok<AGROW>
        Xk(end + 1, 1) = k.X_CO;  Xg(end + 1, 1) = g.X_CO; %#ok<AGROW>
    end
end
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
title(tl, 'Rate-law comparison with membrane (feed H2 held fixed)');
nexttile; bar([Xk Xg]);  grid on;  xticks(1:numel(labels));  xticklabels([]);
ylabel('CO conversion');  legend('K1', 'G2', 'Location', 'southwest');
nexttile; bar(Xk - Xg);  grid on;  xticks(1:numel(labels));  xticklabels(labels);  xtickangle(45);
ylabel('K1 minus G2');

C = d(d.part == "C", :);
figure('Color', 'w');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
title(tl, 'Peak temperature along the bed');
for T0 = [623.15 673.15]
    nexttile; hold on;
    names = ["isothermal" "adiabatic" "U = 10" "U = 50" "U = 300"];
    vals = zeros(5, 2);
    for j = 1:2
        w = [0.095 1];
        r = C(near(C.T_K, T0) & near(C.w_cat, w(j)), :);
        vals(:, j) = [r.T_max(r.wall == "isothermal"); r.T_max(r.wall == "adiabatic"); ...
                      r.T_max(r.wall == "cooled" & near(r.U, 10)); r.T_max(r.wall == "cooled" & near(r.U, 50)); r.T_max(r.wall == "cooled" & near(r.U, 300))];
    end
    plot(1:5, vals(:, 1) - 273.15, 'o', 'MarkerSize', 9, 'LineWidth', 1.5);
    plot(1:5, vals(:, 2) - 273.15, 's', 'MarkerSize', 9, 'LineWidth', 1.5);
    yline(723.15 - 273.15, '--', 'M7 data limit');  yline(743.15 - 273.15, ':', 'window limit');
    xlim([0.5 5.5]);  xticks(1:5);  xticklabels(names);  xtickangle(30);  ylabel('peak T (C)');  grid on;
    title(sprintf('inlet %.0f C', T0 - 273.15));  legend('9.5 wt %', '100 %', 'Location', 'northwest');
end
end

function c = regime_color(r)
switch r
    case "reaction",  c = [0.93 0.69 0.35];
    case "pinch",     c = [0.40 0.62 0.86];
    case "permeance", c = [0.50 0.78 0.55];
    otherwise,        c = [0.80 0.80 0.80];
end
end

function t = regime_tag(r)
switch r
    case "reaction",  t = 'rxn';
    case "pinch",     t = 'pinch';
    case "permeance", t = 'perm';
    otherwise,        t = 'inactive';
end
end

function b = near(a, x)
b = abs(a - x) < 1e-6;
end
