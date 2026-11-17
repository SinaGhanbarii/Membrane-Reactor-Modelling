function T = plot_param_usage(p)
% Matrix of register keys read by each function: rows are functions, columns are groups of parameters.
% Built from the code itself (see param_usage), so it cannot drift from what the functions really read.
T = param_usage(p);
funcs = unique(T.func, 'stable');
groups = unique(T.group, 'stable');
[~, fi] = ismember(T.func, funcs);
[~, gi] = ismember(T.group, groups);
N = accumarray([fi gi], 1, [numel(funcs) numel(groups)]);

figure('Color', 'w', 'Position', [100 100 900 360]);
imagesc(N)
colormap(gca, [linspace(1, 0.45, 64)', linspace(1, 0.75, 64)', linspace(1, 0.78, 64)'])
for r = 1:size(N, 1)
    for c = 1:size(N, 2)
        if N(r, c) > 0, text(c, r, num2str(N(r, c)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    end
end
xticks(1:numel(groups)), xticklabels(groups), xtickangle(30)
yticks(1:numel(funcs)), yticklabels(funcs)
set(gca, 'TickLabelInterpreter', 'none')
title('Register keys read by each function (count)')
end
