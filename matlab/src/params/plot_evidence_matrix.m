function plot_evidence_matrix(p)
%PLOT_EVIDENCE_MATRIX  Which parameter of which model block is reported, inferred or missing.
%   Rows = model blocks and options.  Columns = the parameters each block needs.
%   Colour = kind in params_register.csv.  White = not applicable to that block.
%   Red = placeholder (value not available) or expected but absent from the register.
%   Teal = swept parameter: no single value, bounds given in the register (sweep_lo, sweep_hi).
rows = {
  'K1 rate (Keiski)',      'rate.K1',   {'k0','ln_k0'; 'E','E'; 'n_CO','n_CO'; 'm_H2O','m_H2O'; 'p_CO2','p_CO2'; 'q_H2','q_H2'; 'unit basis','conc_basis_used'}
  'G1 rate (Barbieri)',    'rate.G1',   {'k0','k0'; 'E','E_over_R'; 'n_CO','n_CO'; 'm_H2O','m_H2O'}
  'G2 rate (Hla)',         'rate.G2',   {'k0','log10_k0'; 'E','E'; 'n_CO','n_CO'; 'm_H2O','m_H2O'; 'p_CO2','p_CO2'; 'q_H2','q_H2'; 'unit basis','pressure_unit'}
  'M7 membrane (base)',    'membrane.M7', {'Pm0 / Q0','Pm0'; 'Ea','Ea'; 'n','n'; 'delta','delta'}
  'M1 membrane (thick)',   'membrane.M1', {'Pm0 / Q0','Q0'; 'Ea','Ea'; 'n','n'; 'delta','delta'}
  'M8 membrane (Hla)',     'membrane.M8', {'Pm0 / Q0','Pm0'; 'Ea','Ea'; 'n','n'; 'delta','delta'}
  'M4 membrane (3.2 um)',  'membrane.M4', {'Pm0 / Q0','Q0'; 'Ea','Ea'; 'n','n'; 'delta','delta'}
  'CO inhibition, Barbieri 2008',  'inhibition.barbieri2008',         {'alpha','alpha'; 'K_CO','K_CO'}
  'CO inhibition, Melendez 15 %',  'inhibition.melendez2017_60H2_15CO', {'alpha','alpha'; 'K_CO','K_CO'}
  'CO inhibition, Melendez 5 %',   'inhibition.melendez2017_60H2_5CO',  {'alpha','alpha'; 'K_CO','K_CO'}
  'Catalyst (Hla 2015)',   'catalyst.Hla2015', {'rho_p','rho_p'; 'cp','cp'; 'k (catalyst)','k'; 'd_p','d_p'; 'identity','catalyst_identity'}
  'Reactor (Hla 2015 case)', 'cases.hla2015', {'shell ID','shell_inner_diameter'; 'voidage','bed_voidage'; 'membrane area','membrane_area'; 'steam/carbon','steam_to_carbon'}
  'Reactor (target case)', 'cases.target', {'geometry','geometry'}
  'Swept parameters',      'sweep',     {'U wall','U_wall'; 'd_p','d_p'; 'delta','delta_membrane'; 'voidage','bed_voidage'}
};
cols = {};
for r = 1:size(rows,1)
    for c = 1:size(rows{r,3},1)
        if ~any(strcmp(cols, rows{r,3}{c,1})), cols{end+1} = rows{r,3}{c,1}; end %#ok<AGROW>
    end
end
names = {'reported','calculated','assumed','inferred','placeholder','sweep'};
code = zeros(size(rows,1), numel(cols));
for r = 1:size(rows,1)
    for c = 1:size(rows{r,3},1)
        j = find(strcmp(cols, rows{r,3}{c,1}), 1);
        code(r,j) = kind_code(p, rows{r,2}, rows{r,3}{c,2}, names);
    end
end
colors = [1 1 1; 0.17 0.63 0.17; 0.12 0.47 0.71; 1.00 0.50 0.05; 0.58 0.40 0.74; 0.84 0.15 0.16; 0.00 0.60 0.60];
figure('Name', 'Evidence matrix', 'Color', 'w', 'Position', [60 60 1250 640]);
imagesc(code); colormap(colors); clim([-0.5 6.5]);
set(gca, 'XTick', 1:numel(cols), 'XTickLabel', cols, 'XTickLabelRotation', 45, ...
         'YTick', 1:size(rows,1), 'YTickLabel', rows(:,1), 'TickLength', [0 0]);
hold on;
for x = 0.5:1:numel(cols)+0.5, plot([x x], [0.5 size(rows,1)+0.5], 'Color', 'w', 'LineWidth', 2, 'HandleVisibility', 'off'); end
for y = 0.5:1:size(rows,1)+0.5, plot([0.5 numel(cols)+0.5], [y y], 'Color', 'w', 'LineWidth', 2, 'HandleVisibility', 'off'); end
for j = 1:numel(names)
    scatter(nan, nan, 70, colors(j+1,:), 's', 'filled', 'DisplayName', names{j});
end
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
title('Evidence status of each parameter (white = not applicable)');
end

function c = kind_code(p, base, leaf, names)
% 0 = not applicable; 1-5 = index into names; an expected leaf that is absent counts as placeholder
if isempty(base) || isempty(leaf), c = 5; return; end
s = p; parts = strsplit([base '.' leaf], '.');
for k = 1:numel(parts)
    if ~isstruct(s) || ~isfield(s, parts{k}), c = 5; return; end
    s = s.(parts{k});
end
if ~(isstruct(s) && isfield(s, 'kind')), c = 5; return; end
j = find(strcmp(names, s.kind), 1);
if isempty(j), c = 5; else, c = j; end
if isempty(s.value) && ~strcmp(s.kind, 'sweep'), c = 5; end
end
