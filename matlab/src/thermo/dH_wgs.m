function dH = dH_wgs(T)
%DH_WGS  Reaction enthalpy of WGS at T [K], in J/mol (negative = exothermic).

c = pbmr_const();
[~, h] = nasa_props(T);
dH = c.NU' * h;
end
