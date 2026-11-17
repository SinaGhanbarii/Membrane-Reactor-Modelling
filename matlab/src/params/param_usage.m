function T = param_usage(p)
% Which function reads which register key, found by scanning the source for pget(p, ...) calls.
% A key built from an option, such as ['rate.' law '.n_CO'], is shown as rate.<law>.n_CO.
% kind is looked up in the register for fixed keys; keys that depend on an option say "by option".
root = fileparts(fileparts(mfilename('fullpath')));      % matlab/src
files = dir(fullfile(root, '**', '*.m'));
func = strings(0, 1);  key = strings(0, 1);
for i = 1:numel(files)
    txt = fileread(fullfile(files(i).folder, files(i).name));
    txt = regexprep(txt, '^\s*%[^\n]*', '', 'lineanchors');          % ignore comment lines
    tok = regexp(txt, 'pget\(\s*p\s*,\s*(\[[^\]]*\]|''[^'']*'')', 'tokens');
    for k = 1:numel(tok)
        e = tok{k}{1};
        e = regexprep(e, '''\s+([A-Za-z_]\w*(\{\w+\})?)\s*''', '<$1>');
        e = regexprep(e, '''\s+([A-Za-z_]\w*(\{\w+\})?)\]', '<$1>]');
        e = regexprep(e, '[\[\]'']', '');
        func(end + 1, 1) = string(files(i).name(1:end - 2));   %#ok<AGROW>
        key(end + 1, 1) = string(e);                            %#ok<AGROW>
    end
end
group = strings(size(key));  kind = strings(size(key));
for i = 1:numel(key)
    parts = strsplit(key(i), '.');
    group(i) = strjoin(parts(1:min(2, numel(parts))), '.');
    kind(i) = kind_of(p, parts);
end
T = unique(table(func, key, group, kind));
T = sortrows(T, {'func', 'key'});
end

function k = kind_of(p, parts)
k = "by option";
if any(contains(parts, '<')), return; end
x = p;
for j = 1:numel(parts)
    if ~isstruct(x) || ~isfield(x, parts{j}), k = "not found"; return; end
    x = x.(parts{j});
end
if isstruct(x) && isfield(x, 'kind'), k = string(x.kind); end
end
