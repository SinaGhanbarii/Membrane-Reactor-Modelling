function [dy, d] = reactor_rhs(~, y, p, c)
% d/dz of y = [F (6 mol/s); H2 in the permeate stream G (mol/s); T (K); P (Pa)] on the retentate side.
% Permeate H2 partial pressure is P_perm * G / (G + sweep); with no sweep it is P_perm.
k = pbmr_const();
F = y(1:6);  G = max(y(7), 0);  T = y(8);  P = y(9);
yi = F / sum(F);

ro = struct('eta', c.eta);
if strcmp(c.eta, 'on'), ro.d_p = c.d_p;  ro.D_eff = c.D_eff;  ro.rho_p = c.rho_p; end
[r, info] = rate_wgs(p, c.law, T, P, yi, ro);
rxn = c.rho_cat * r * c.A_bed;                       % mol/m/s

pperm = c.P_perm;
if c.sweep > 0, pperm = c.P_perm * G / (G + c.sweep); end
J = membrane_flux(p, c.membrane, c.inhibition, T, yi(k.idx.H2) * P, pperm, yi(k.idx.CO) * P / 1e3);
perm = c.perm_scale * max(J, 0) * pi * c.D_m;        % mol/m/s, no back-permeation

dF = k.NU * rxn;
dF(k.idx.H2) = dF(k.idx.H2) - perm;

Qrxn = -dH_wgs(T) * rxn;                             % W/m
Qwall = 0;
if strcmp(c.wall, 'cooled'), Qwall = c.U * pi * c.D_s * (T - c.T_wall); end
dT = 0;
if ~strcmp(c.wall, 'isothermal')
    dT = (Qrxn - Qwall) / (max(F, 0)' * nasa_props(T));
end

dP = 0;
if strcmp(c.dp, 'on')
    M = (k.ATOMS' * [12.011; 1.008; 15.999; 14.007]) / 1000;   % kg/mol per species
    rho = P * (yi' * M) / (k.R * T);
    u = sum(F) * k.R * T / (P * c.A_bed);
    dP = -(150 * c.mu * (1 - c.eps)^2 / (c.eps^3 * c.d_p^2) * u ...
           + 1.75 * (1 - c.eps) / (c.eps^3 * c.d_p) * rho * u^2);
end

sgn = 1;
if strcmp(c.flow, 'counter'), sgn = -1; end          % permeate flows towards z = 0
dy = [dF; sgn * perm; dT; dP];
if nargout > 1
    d = struct('rxn', rxn, 'perm', perm, 'pperm', pperm, 'Qrxn', Qrxn, 'Qwall', Qwall, ...
               'beta', info.beta, 'eta', info.eta, 'flags', {info.flags});
end
end
