function J = membrane_flux(p, opt, inh, T, pH2_ret, pH2_perm, pCO_kPa)
%MEMBRANE_FLUX  H2 flux (mol m-2 s-1) through the Pd membrane, Sieverts form.
%   J = Q0 * exp(-Ea/(R T)) * (pr^n - pp^n) * (1 - alpha*theta)
%   theta = K_CO*pCO/(1 + K_CO*pCO), pCO in kPa.
%   opt : membrane option name ('M7','M1','M8','M4')
%   inh : inhibition option name ('off','barbieri2008','melendez2017_60H2_15CO',...)
%   T [K]; pH2_ret, pH2_perm [Pa]; pCO_kPa [kPa] CO partial pressure on the retentate side.
%   Concentration polarisation and H2O/CO2 effects are NOT modelled (see assumptions).
c = pbmr_const();
m = p.membrane.(opt);
delta = pget(p, ['membrane.' opt '.delta']);
n     = pget(p, ['membrane.' opt '.n']);
Ea    = pget(p, ['membrane.' opt '.Ea']);
if isfield(m, 'Pm0')
    Q0 = pget(p, ['membrane.' opt '.Pm0']) / delta;
else
    Q0 = pget(p, ['membrane.' opt '.Q0']);
end
alpha = pget(p, ['inhibition.' inh '.alpha']);
K     = pget(p, ['inhibition.' inh '.K_CO']);
theta = K .* pCO_kPa ./ (1 + K .* pCO_kPa);
drive = max(pH2_ret, 0).^n - max(pH2_perm, 0).^n;
free = max(1 - alpha .* theta, 0);        % a set with alpha > 1 would otherwise give a negative flux at high CO
J = Q0 .* exp(-Ea ./ (c.R .* T)) .* drive .* free;
end
