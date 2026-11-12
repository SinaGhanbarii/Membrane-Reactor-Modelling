function E = element_flows(F)
%ELEMENT_FLOWS  Element flows [mol/s] (4x1: C, H, O, N) for species flows F.

c = pbmr_const();
E = c.ATOMS * F(:);
end
