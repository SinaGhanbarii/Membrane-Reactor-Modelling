function plot_param_relations(p)
%PLOT_PARAM_RELATIONS  Only the parameters that are linked to others, and the links between them.
%   Reads the "related" column of params_register.csv.  Nodes are coloured by evidence kind;
%   arrows run from the parameter that carries the link to the one it points at.
%   reuses blue, conflicts_with red, needs black, inferred_from purple, range_from teal,
%   maybe_affected_by orange (a hypothesis, not a finding).  Black nodes are sources outside the register.
edges = collect(p, '', {});
if isempty(edges), disp('No relations recorded in the register.'); return; end
ids = unique([edges(:,1); edges(:,3)], 'stable');
[~, si] = ismember(edges(:,1), ids);
[~, ti] = ismember(edges(:,3), ids);
G = digraph(si, ti);

kn = {'reported','calculated','assumed','inferred','placeholder','sweep','external'};
kc = [0.17 0.63 0.17; 0.12 0.47 0.71; 1.00 0.50 0.05; 0.58 0.40 0.74; 0.84 0.15 0.16; 0.00 0.60 0.60; 0 0 0];
rn = {'reuses','conflicts_with','needs','inferred_from','range_from','maybe_affected_by'};
rc = [0.12 0.47 0.71; 0.84 0.15 0.16; 0 0 0; 0.58 0.40 0.74; 0.00 0.60 0.60; 1.00 0.50 0.05];

NC = zeros(numel(ids), 3);
nodekind = cell(numel(ids), 1);
for i = 1:numel(ids)
    nodekind{i} = leaf_kind(p, ids{i});
    j = find(strcmp(kn, nodekind{i}), 1); if isempty(j), j = numel(kn); end
    NC(i, :) = kc(j, :);
end
EC = zeros(size(edges, 1), 3);
for e = 1:size(edges, 1)
    j = find(strcmp(rn, edges{e,2}), 1); if isempty(j), j = 3; end
    EC(e, :) = rc(j, :);
end

figure('Name', 'Parameter relations', 'Color', 'w', 'Position', [60 60 1300 760]);
h = plot(G, 'Layout', 'force', 'Iterations', 300, 'NodeLabel', ids, 'Interpreter', 'none', ...
         'ArrowSize', 12, 'LineWidth', 2);
h.NodeColor = NC; h.EdgeColor = EC; h.MarkerSize = 9; h.NodeFontSize = 9;
axis off; hold on;
lh = []; ln = {};
for j = 1:numel(kn)
    if any(strcmp(nodekind, kn{j}))
        lh(end+1) = scatter(nan, nan, 60, kc(j,:), 'filled'); ln{end+1} = kn{j}; %#ok<AGROW>
    end
end
for j = 1:numel(rn)
    if any(strcmp(edges(:,2), rn{j}))
        lh(end+1) = plot(nan, nan, '-', 'Color', rc(j,:), 'LineWidth', 2); ln{end+1} = ['link: ' rn{j}]; %#ok<AGROW>
    end
end
legend(lh, ln, 'Location', 'southoutside', 'NumColumns', 6, 'Box', 'off', 'Interpreter', 'none');
title(sprintf('Parameter relations (%d links, %d parameters)', size(edges, 1), numel(ids)));
end

function edges = collect(s, path, edges)
if isstruct(s) && isfield(s, 'value')
    if ~isempty(s.related)
        items = strsplit(s.related, ';');
        for q = 1:numel(items)
            ty = strsplit(strtrim(items{q}), ':');
            edges(end+1, :) = {path(2:end), ty{1}, strjoin(ty(2:end), ':')}; %#ok<AGROW>
        end
    end
    return
end
if isstruct(s) && isscalar(s)
    fn = fieldnames(s);
    for k = 1:numel(fn)
        if strcmp(fn{k}, 'validity'), continue; end
        edges = collect(s.(fn{k}), [path '.' fn{k}], edges);
    end
end
end

function k = leaf_kind(p, id)
k = 'external';
parts = strsplit(id, '.'); x = p;
for i = 1:numel(parts)
    if ~isstruct(x) || ~isfield(x, parts{i}), return; end
    x = x.(parts{i});
end
if isstruct(x) && isfield(x, 'kind'), k = x.kind; end
end
