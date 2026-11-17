function w = membrane_validity(p, opt, inh, T, P_ret, P_perm, yCO)
%MEMBRANE_VALIDITY  Ways a condition falls outside the data range of a membrane set and its inhibition set.
%   opt: membrane option ('M7', ...)   inh: inhibition option ('off', ...)
%   T [K], P_ret and P_perm [Pa], yCO = CO mole fraction on the retentate side.
%   Uses the register ranges: membrane T_K, dP_bar (retentate minus permeate pressure), P_kPa_ret;
%   inhibition T_K, yCO. Returns a cell array of warnings; empty means inside the ranges.
w = {};
w = [w, rng(p.membrane.(opt), 'T_K',       T,                   'T',            'K',   ['membrane ' opt])];
w = [w, rng(p.membrane.(opt), 'dP_bar',    (P_ret - P_perm)/1e5, 'dP',           'bar', ['membrane ' opt])];
w = [w, rng(p.membrane.(opt), 'P_kPa_ret', P_ret/1e3,           'P retentate',  'kPa', ['membrane ' opt])];
if ~strcmp(inh, 'off')
    w = [w, rng(p.inhibition.(inh), 'T_K', T,   'T',   'K', ['membrane ' opt]), ...
            rng(p.inhibition.(inh), 'yCO', yCO, 'CO',  '-', ['membrane ' opt])];
end
end

function w = rng(s, field, x, name, unit, who)
w = {};
if ~isfield(s, 'validity') || ~isfield(s.validity, field), return; end
r = s.validity.(field);
tol = 1e-3 * max(abs(r));                 % the register rounds some ranges (723 for 723.15 K)
if x < min(r) - tol || x > max(r) + tol
    w{1} = sprintf('EXTRAPOLATED (%s): %s = %g %s outside data range [%g, %g]', who, name, x, unit, min(r), max(r));
end
end
