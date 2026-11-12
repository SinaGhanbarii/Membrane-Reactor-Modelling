function w = check_validity(valid, T, P_bar, yCO_dry)
%CHECK_VALIDITY  List the ways a condition falls outside a parameter set's data range.
%   valid   : the 'validity' struct of a rate law, membrane or inhibition set
%   T [K], P_bar [bar], yCO_dry [-]  (pass [] to skip a check)
%   Returns a cell array of warning strings; empty means inside the range.
w = {};
w = [w, rng(valid, 'T_K',      T,       'T',      'K')];
w = [w, rng(valid, 'P_bar',    P_bar,   'P',      'bar')];
w = [w, rng(valid, 'yCO_dry',  yCO_dry, 'CO dry', '-')];
end

function w = rng(valid, field, x, name, unit)
w = {};
if isempty(x) || ~isfield(valid, field), return; end
r = valid.(field);
if any(x < min(r) - 1e-12 | x > max(r) + 1e-12)
    w{1} = sprintf('EXTRAPOLATED: %s = %g %s outside data range [%g, %g]', ...
                   name, x(find(x < min(r) | x > max(r), 1)), unit, min(r), max(r));
end
end
