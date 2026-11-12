function F = wgs_equilibrium(F0, T)
%WGS_EQUILIBRIUM  Outlet molar flows [mol/s] (6x1) at WGS-only equilibrium.

c = pbmr_const();
F = F0(:) + c.NU * wgs_extent_eq(F0, T);
end
