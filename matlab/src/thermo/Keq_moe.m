function K = Keq_moe(T)
%KEQ_MOE  Empirical WGS equilibrium constant  K = exp(4577.8/T - 4.33).
%   Moe (1962), Chem. Eng. Prog. 58, 33-36, as quoted by Helmi et al. (2016),
%   Molecules 21, 376 (their ref. [3]).  Moe (1962) itself has not been read;
%   use only as an independent cross-check of Keq_nasa, valid for the
%   temperature range of the original correlation (not verified here).

K = exp(4577.8 ./ T - 4.33);
end
