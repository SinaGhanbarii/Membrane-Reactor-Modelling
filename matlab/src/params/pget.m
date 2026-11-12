function v = pget(p, path)
%PGET  Value of a parameter leaf by dotted path, e.g. pget(p,'membrane.M7.Pm0').
%   Stops with an error if the leaf is missing or its value is null/empty,
%   so a run can never continue on a number that no source provides.
parts = strsplit(path, '.');
x = p;
for k = 1:numel(parts)
    if ~isstruct(x) || ~isfield(x, parts{k})
        error('pget:nopath', 'No parameter "%s" (failed at "%s").', path, parts{k});
    end
    x = x.(parts{k});
end
if isstruct(x) && isfield(x, 'value')
    if isempty(x.value)
        error('pget:missing', ['Parameter "%s" has no value (%s). ' ...
              'It is not available in the sources; supply it before running.'], path, x.source);
    end
    v = x.value;
else
    v = x;
end
end
