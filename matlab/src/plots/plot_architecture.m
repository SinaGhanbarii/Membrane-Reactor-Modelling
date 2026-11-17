function plot_architecture()
% Layer diagram of the model code: each layer uses only the layers above it.
rows = { ...
    'DATA', {'params_register.csv|value, unit, source, kind', 'nasa7_species.json|species thermodynamics'}, [0.92 0.93 0.94]; ...
    'PARAMETERS', {'load_params, pget|read register; stop if value missing', 'pbmr_const|species order, stoichiometry, R', 'check_validity|EXTRAPOLATED warning'}, [0.88 0.91 0.96]; ...
    'PHYSICS', {'thermo|properties, equilibrium, heat', 'kinetics|rate_wgs (K1, G2), eta', 'membrane|membrane_flux (Sieverts)'}, [0.86 0.92 0.98]; ...
    'REACTOR MODEL', {'make_case|defaults + overrides', 'solvers|co- and counter-current', 'reactor_rhs|species, permeate, T, P'}, [0.84 0.94 0.94]; ...
    'ONE CASE', {'run_case|case, profiles, summary', 'regime_labels|no flux, reaction, pinch, permeance, idle'}, [0.93 0.90 0.97]; ...
    'STUDIES', {'studies|loop over cases', 'results/*.csv|one file per study', 'plots|figures'}, [0.93 0.90 0.97]};

figure('Color', 'w', 'Position', [100 100 1000 720]);
ax = axes('Position', [0 0 1 1]);
axis(ax, [0 100 0 100]); axis off; hold on
h = 11;  gap = 4;  top = 97;
for r = 1:size(rows, 1)
    y = top - r * h - (r - 1) * gap;
    text(1, y + h + 1.2, rows{r, 1}, 'FontSize', 8, 'Color', [0.35 0.4 0.45], 'FontWeight', 'bold');
    items = rows{r, 2};
    n = numel(items);
    w = (98 - 2 * (n - 1)) / n;
    for k = 1:n
        x = 1 + (k - 1) * (w + 2);
        rectangle('Position', [x y w h], 'Curvature', 0.08, 'FaceColor', rows{r, 3}, 'EdgeColor', [0.6 0.67 0.72]);
        parts = strsplit(items{k}, '|');
        text(x + 1, y + h - 3.2, parts{1}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        text(x + 1, y + 3.4, parts{2}, 'FontSize', 8, 'Color', [0.3 0.35 0.4], 'Interpreter', 'none');
    end
    if r < size(rows, 1)
        plot([50 50], [y y - gap + 0.6], '-', 'Color', [0.4 0.45 0.5], 'LineWidth', 1.2);
        plot(50, y - gap + 0.6, 'v', 'Color', [0.4 0.45 0.5], 'MarkerFaceColor', [0.4 0.45 0.5], 'MarkerSize', 5);
    end
end
end
