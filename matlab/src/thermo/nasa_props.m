function [cp, h, s, g0] = nasa_props(T)
%NASA_PROPS  Ideal-gas standard-state properties of the six species at T [K].
%   cp  [J/(mol K)]   heat capacity
%   h   [J/mol]       absolute (formation-based) enthalpy
%   s   [J/(mol K)]   standard entropy at P0 = 1 atm
%   g0  [J/mol]       standard Gibbs energy at P0 = h - T*s
%   All outputs are 6x1 column vectors in the species order of pbmr_const.
%
%   NASA-7 form (a = a1..a7):
%     cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
%     h/RT = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
%     s/R  = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7

c  = pbmr_const();
sp = load_species();
R  = c.R;
n  = numel(sp.T_mid);

cp = zeros(n, 1);  h = zeros(n, 1);  s = zeros(n, 1);
for i = 1:n
    if T < sp.T_low(i) || T > sp.T_high(i)
        error('nasa_props:range', 'T = %.1f K outside NASA range [%g, %g] K for %s.', ...
              T, sp.T_low(i), sp.T_high(i), c.species{i});
    end
    if T < sp.T_mid(i)
        a = sp.a_low(i, :);
    else
        a = sp.a_high(i, :);
    end
    cp(i) = R * (a(1) + a(2)*T + a(3)*T^2 + a(4)*T^3 + a(5)*T^4);
    h(i)  = R * T * (a(1) + a(2)*T/2 + a(3)*T^2/3 + a(4)*T^3/4 + a(5)*T^4/5 + a(6)/T);
    s(i)  = R * (a(1)*log(T) + a(2)*T + a(3)*T^2/2 + a(4)*T^3/3 + a(5)*T^4/4 + a(7));
end
g0 = h - T .* s;
end
