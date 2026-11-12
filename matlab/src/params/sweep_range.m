function r = sweep_range(p, path)
%SWEEP_RANGE  [lo; hi] sweep bounds of a parameter, e.g. sweep_range(p,'sweep.U_wall').
%   Stops with an error if the parameter has no bounds in the register.
parts = strsplit(path, '.');
x = p;
for k = 1:numel(parts)
    if ~isstruct(x) || ~isfield(x, parts{k})
        error('sweep_range:nopath', 'No parameter "%s".', path);
    end
    x = x.(parts{k});
end
if ~isstruct(x) || ~isfield(x, 'sweep') || isempty(x.sweep)
    error('sweep_range:none', 'Parameter "%s" has no sweep bounds in the register.', path);
end
r = x.sweep;
end
