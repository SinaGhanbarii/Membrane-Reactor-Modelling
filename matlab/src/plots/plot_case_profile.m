function plot_case_profile(out)
% Along-bed view of one run_case result: conversion and recovery, beta, pinch ratio, and the regime strip.
s = out.sol;  pr = out.profile;  c = out.case;
X = 1 - s.F(:, 1) / c.F_in(1);
rec = s.Fp / (c.F_in(4) + c.F_in(1));   % counter-current: H2 in the permeate stream at z, which leaves at z = 0

figure('Color', 'w');
tl = tiledlayout(4, 1, 'TileSpacing', 'compact');
title(tl, sprintf('%s, %.0f K inlet, %s wall, %s, L = %.2f m, permeate %.2f bar', ...
      c.law, c.T_in, c.wall, c.membrane, c.L, c.P_perm / 1e5));
xlabel(tl, 'position along the bed (m)');

nexttile; plot(s.z, X, s.z, rec, 'LineWidth', 1.5); grid on;
ylabel('fraction'); legend('CO conversion', 'H2 in permeate stream', 'Location', 'southeast');
nexttile; plot(s.z, pr.beta, 'LineWidth', 1.5); yline(0.9, '--'); grid on;
ylabel('beta'); ylim([0 1.05]);
nexttile; semilogy(s.z, pr.pinch_ratio, 'LineWidth', 1.5); yline(1.05, '--'); yline(1, ':'); grid on;
ylabel('pH2 / permeate pressure');
nexttile;
if any(pr.regime > 0)
    imagesc(s.z, 1, pr.regime');
    colormap(gca, [0.80 0.80 0.80; 0.93 0.69 0.35; 0.40 0.62 0.86; 0.50 0.78 0.55; 0.82 0.78 0.90]);
    clim([0.5 5.5]);
    cb = colorbar('Location', 'southoutside');
    cb.Ticks = 1:5;  cb.TickLabels = out.summary.regime_names;
    yticks([]);
else
    text(0.5, 0.5, 'no membrane: no regime labels', 'HorizontalAlignment', 'center');
    axis off;
end
end
