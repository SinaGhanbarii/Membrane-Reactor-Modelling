function plot_flow_compare(p, ov)
% Co-current against counter-current with the same sweep gas: permeate H2 partial pressure, retentate H2
% partial pressure, H2 recovered along the bed, and the regime strips. ov sets the case (default M7, 400 C, 1 mmol/s).
if nargin < 2, ov = struct('sweep', 1e-3); end
flows = {'cocurrent', 'counter'};
tags = {'co', 'counter'};
o = cell(1, 2);
for i = 1:2
    ov.flow = flows{i};
    o{i} = run_case(p, ov);
end
c = o{1}.case;
colours = [0.00 0.45 0.74; 0.85 0.33 0.10];

figure('Color', 'w', 'Position', [100 100 760 900]);
tl = tiledlayout(5, 1, 'TileSpacing', 'compact');
title(tl, sprintf('%s, %.0f K, %.2f m, sweep %.1f mmol/s, permeate %.2f bar', ...
      c.membrane, c.T_in, c.L, c.sweep * 1e3, c.P_perm / 1e5));
xlabel(tl, 'position along the bed (m)');

nexttile; hold on;
for i = 1:2
    s = o{i}.sol;
    plot(s.z, s.F(:, 4) ./ sum(s.F, 2) .* s.P / 1e5, 'Color', colours(i, :), 'LineWidth', 1.5);
end
for i = 1:2
    plot(o{i}.sol.z, o{i}.profile.pperm / 1e5, '--', 'Color', colours(i, :), 'LineWidth', 1.5);
end
grid on; xlim([0 c.L]); ylabel('pH2 (bar)');
legend('retentate, co-current', 'retentate, counter', 'permeate, co-current', 'permeate, counter', 'Location', 'northeast', 'FontSize', 7);

nexttile; hold on;
for i = 1:2
    plot(o{i}.sol.z, 1 - o{i}.sol.F(:, 1) / c.F_in(1), 'Color', colours(i, :), 'LineWidth', 1.5);
end
grid on; xlim([0 c.L]); ylabel('CO conv.'); legend('co-current', 'counter', 'Location', 'southeast', 'FontSize', 7);

nexttile; hold on;
for i = 1:2
    plot(o{i}.sol.z, o{i}.profile.perm, 'Color', colours(i, :), 'LineWidth', 1.5);
end
grid on; xlim([0 c.L]); ylabel('H2 flux');

for i = 1:2
    nexttile;
    imagesc(o{i}.sol.z, 1, o{i}.profile.regime');
    colormap(gca, [0.80 0.80 0.80; 0.93 0.69 0.35; 0.40 0.62 0.86; 0.50 0.78 0.55; 0.82 0.78 0.90]);
    clim([0.5 5.5]); yticks(1); yticklabels(tags(i));
end
cb = colorbar('Location', 'southoutside');
cb.Ticks = 1:5;  cb.TickLabels = o{1}.summary.regime_names;
end
