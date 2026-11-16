function out = run_case(p, ov)
% One case in, one result out. ov is a struct of overrides (see make_case).
%   H2 recovery = permeated H2 / (H2 fed + CO fed).
%   X_eq = CO conversion at WGS equilibrium of the feed at the outlet temperature, without a membrane.
if nargin < 2, ov = struct(); end
c = make_case(p, ov);
if strcmp(c.flow, 'counter'), s = solve_countercurrent(p, c); else, s = solve_cocurrent(p, c); end
Fpo = s.Fp(end);                                       % H2 leaving in the permeate stream
if strcmp(c.flow, 'counter'), Fpo = s.Fp(1); end

n = numel(s.z);
names = {'rxn', 'perm', 'pperm', 'Qrxn', 'Qwall', 'beta', 'eta'};
prof = struct();
for j = 1:numel(names), prof.(names{j}) = zeros(n, 1); end
flags = {};
for i = 1:n
    [~, d] = reactor_rhs(0, [s.F(i, :)'; s.Fp(i); s.T(i); s.P(i)], p, c);
    for j = 1:numel(names), prof.(names{j})(i) = d.(names{j}); end
    flags = [flags, d.flags];   %#ok<AGROW>
end

pr = s.F(:, 4) ./ sum(s.F, 2) .* s.P ./ max(prof.pperm, 1e-9 * c.P_perm);   % retentate over permeate H2 partial pressure
prof.pinch_ratio = pr;
if c.perm_scale == 0
    prof.regime = zeros(n, 1);
    frac = nan(5, 1);
else
    idle = abs(prof.rxn) < 1e-2 * max(abs(prof.rxn)) & prof.perm < 1e-2 * max(prof.perm);
    prof.regime = regime_labels(prof.beta, pr, prof.perm > 0, idle);
    frac = arrayfun(@(k) mean(prof.regime == k), 1:5)';   % share of bed length (uniform grid)
end

Fout = s.F(end, :)';
Eout = element_flows(Fout) + element_flows([0 0 0 1 0 0]') * Fpo;
sm.X_CO = co_conversion(c.F_in, Fout);
sm.X_eq = co_conversion(c.F_in, wgs_equilibrium(c.F_in, s.T(end)));
sm.H2_recovery = Fpo / (c.F_in(4) + c.F_in(1));
sm.permeate_H2_fraction = 1;                           % H2 share of the permeate stream (H2 + sweep)
if c.sweep > 0, sm.permeate_H2_fraction = Fpo / (Fpo + c.sweep); end
sm.atom_error = max(abs(Eout - element_flows(c.F_in)) ./ element_flows(c.F_in));
sm.T_out = s.T(end);
sm.P_out = s.P(end);
sm.beta_out = prof.beta(end);
sm.regime_names = {'no flux', 'reaction', 'pinch', 'permeance', 'idle'};
sm.regime_fraction = frac;
sm.flags = unique(flags);
out = struct('case', c, 'sol', s, 'profile', prof, 'summary', sm);
end
