function param_browser(p)
%PARAM_BROWSER  Expandable tree of the parameter struct; click a leaf to see its details.
%   param_browser(load_params())
fig = uifigure('Name', 'Parameter browser', 'Position', [100 100 1000 600]);
g = uigridlayout(fig, [1 2]); g.ColumnWidth = {'1x', '1.6x'};
tr = uitree(g);
ta = uitextarea(g, 'Editable', 'off', 'FontName', 'Consolas');
add_nodes(tr, p);
tr.SelectionChangedFcn = @(~, e) set(ta, 'Value', describe(e.SelectedNodes.NodeData));
end

function add_nodes(parent, s)
fn = fieldnames(s);
for i = 1:numel(fn)
    v = s.(fn{i});
    if isstruct(v) && ~isfield(v, 'value')
        n = uitreenode(parent, 'Text', fn{i});
        add_nodes(n, v);
    elseif isstruct(v)
        uitreenode(parent, 'Text', sprintf('%s = %s   [%s]', fn{i}, valstr(v.value), v.kind), 'NodeData', v);
    elseif isnumeric(v)
        uitreenode(parent, 'Text', sprintf('%s = %s   [range]', fn{i}, valstr(v)), 'NodeData', v);
    else
        uitreenode(parent, 'Text', sprintf('%s: %s', fn{i}, char(v)), 'NodeData', char(v));
    end
end
end

function t = valstr(v)
if isempty(v), t = '(missing)';
elseif isnumeric(v) && numel(v) <= 4, t = mat2str(v(:)', 5);
elseif isnumeric(v), t = sprintf('[%d values]', numel(v));
else, t = char(v); end
end

function lines = describe(d)
if isstruct(d)
    lines = {['value:  ' valstr(d.value)], ['unit:   ' d.unit], ['kind:   ' d.kind], ...
             ['source: ' d.source], ['note:   ' d.note]};
elseif isnumeric(d)
    lines = {['range [min; max]: ' mat2str(d(:)', 6)]};
else
    lines = {char(d)};
end
end
