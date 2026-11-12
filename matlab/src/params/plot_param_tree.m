function plot_param_tree(p, root, showRanges)
%PLOT_PARAM_TREE  Parameter tree with relation links, nodes coloured by evidence kind.
%   plot_param_tree(p)            whole tree (crowded; validity ranges hidden)
%   plot_param_tree(p, '', true)  whole tree including validity-range nodes
%   plot_param_tree(p, 'rate')    one branch, e.g. 'rate', 'membrane', 'cases.hla2015', 'sweep'
%   Nodes:   reported green, calculated blue, assumed orange, inferred purple, placeholder red,
%            sweep teal, validity range grey, group dark grey, external source black.
%   Dashed lines are the "related" column of the register (only drawn when both ends are shown):
%            reuses blue, conflicts_with red, needs black, inferred_from purple,
%            range_from teal, maybe_affected_by orange (a hypothesis, not a finding).
if nargin < 2, root = ''; end
if nargin < 3, showRanges = false; end
s = p; name = 'params'; prefix = '';
if ~isempty(root)
    parts = strsplit(root, '.');
    s = getfield(p, parts{:}); %#ok<GFLD>
    name = root; prefix = [root '.'];
end
T = struct('ids', {{''}}, 'labels', {{name}}, 'kinds', {{'group'}}, 'src', [], 'dst', [], 'rel', {{}});
T = walk(s, 1, prefix, T, showRanges);

% external relation targets (not part of the register) become extra nodes
relsrc = []; reltype = {}; reltgt = [];
for i = 1:numel(T.rel)
    r = T.rel{i};                                  % {source node, type, target id}
    j = find(strcmp(T.ids, r{3}), 1);
    if isempty(j)
        top = strsplit(r{3}, '.'); top = top{1};
        if isfield(p, top), continue; end          % register id that is outside this branch: skip
        T.ids{end+1} = r{3}; T.labels{end+1} = r{3}; T.kinds{end+1} = 'external';
        j = numel(T.ids); T.src(end+1) = r{1}; T.dst(end+1) = j;
    end
    relsrc(end+1) = r{1}; reltype{end+1} = r{2}; reltgt(end+1) = j; %#ok<AGROW>
end

G = digraph(T.src, T.dst);
kn = {'reported','calculated','assumed','inferred','placeholder','sweep','range','group','external'};
kc = [0.17 0.63 0.17; 0.12 0.47 0.71; 1.00 0.50 0.05; 0.58 0.40 0.74; 0.84 0.15 0.16; ...
      0.00 0.60 0.60; 0.65 0.65 0.65; 0.25 0.25 0.25; 0 0 0];
C = zeros(numel(T.kinds), 3);
for i = 1:numel(T.kinds)
    j = find(strcmp(kn, T.kinds{i}), 1); if isempty(j), j = 8; end
    C(i, :) = kc(j, :);
end
figure('Name', 'Parameter tree', 'Color', 'w', 'Position', [60 60 1300 760]);
h = plot(G, 'Layout', 'layered', 'Direction', 'right', 'NodeLabel', T.labels, 'Interpreter', 'none');
h.NodeColor = C; h.EdgeColor = [0.8 0.8 0.8]; h.MarkerSize = 6; h.NodeFontSize = 8;
axis off; hold on;

rn = {'reuses','conflicts_with','needs','inferred_from','range_from','maybe_affected_by'};
rc = [0.12 0.47 0.71; 0.84 0.15 0.16; 0 0 0; 0.58 0.40 0.74; 0.00 0.60 0.60; 1.00 0.50 0.05];
for k = 1:numel(relsrc)
    j = find(strcmp(rn, reltype{k}), 1); if isempty(j), j = 3; end
    a = [h.XData(relsrc(k)) h.YData(relsrc(k))]; b = [h.XData(reltgt(k)) h.YData(reltgt(k))];
    d = b - a; len = hypot(d(1), d(2)); if len == 0, continue; end
    nv = [-d(2) d(1)] / len; if nv(1) < 0, nv = -nv; end        % bulge to the right
    c = (a + b) / 2 + nv * (0.30 * len + 0.01 * k);              % control point; small shift separates parallel arcs
    tt = linspace(0, 1, 40)';
    xy = (1 - tt).^2 * a + 2 * (1 - tt) .* tt * c + tt.^2 * b;   % quadratic Bezier arc
    plot(xy(:,1), xy(:,2), '--', 'Color', rc(j,:), 'LineWidth', 1.6, 'HandleVisibility', 'off');
end
lh = []; ln = {};
for j = 1:numel(kn)
    if any(strcmp(T.kinds, kn{j}))
        lh(end+1) = scatter(nan, nan, 50, kc(j,:), 'filled'); ln{end+1} = kn{j}; %#ok<AGROW>
    end
end
for j = 1:numel(rn)
    if any(strcmp(reltype, rn{j}))
        lh(end+1) = plot(nan, nan, '--', 'Color', rc(j,:), 'LineWidth', 1.6); ln{end+1} = ['link: ' rn{j}]; %#ok<AGROW>
    end
end
legend(lh, ln, 'Location', 'southoutside', 'NumColumns', 6, 'Box', 'off', 'Interpreter', 'none');
title(sprintf('Parameter register: %s', name), 'Interpreter', 'none');
end

function T = walk(s, parent, prefix, T, showRanges)
fn = fieldnames(s);
for i = 1:numel(fn)
    v = s.(fn{i});
    if ischar(v) || isstring(v), continue; end          % plain-text meta rows are not drawn
    if ~(isstruct(v) && isfield(v, 'value')) && ~isstruct(v) && ~showRanges, continue; end   % validity ranges
    if isstruct(v) && ~isfield(v, 'value') && ~has_content(v, showRanges), continue; end     % empty groups
    id = [prefix fn{i}];
    T.ids{end+1} = id; T.labels{end+1} = fn{i}; idx = numel(T.ids);
    T.src(end+1) = parent; T.dst(end+1) = idx;
    if isstruct(v) && isfield(v, 'value')
        T.kinds{end+1} = v.kind;
        if ~isempty(v.related)
            items = strsplit(v.related, ';');
            for q = 1:numel(items)
                ty = strsplit(strtrim(items{q}), ':');
                T.rel{end+1} = {idx, ty{1}, strjoin(ty(2:end), ':')};
            end
        end
    elseif isstruct(v)
        T.kinds{end+1} = 'group';
        T = walk(v, idx, [id '.'], T, showRanges);
    else
        T.kinds{end+1} = 'range';
    end
end
end

function tf = has_content(s, showRanges)
tf = false;
fn = fieldnames(s);
for i = 1:numel(fn)
    v = s.(fn{i});
    if ischar(v) || isstring(v), continue; end
    if isstruct(v) && isfield(v, 'value'), tf = true; return; end
    if isstruct(v), if has_content(v, showRanges), tf = true; return; end
    elseif showRanges, tf = true; return; end
end
end
