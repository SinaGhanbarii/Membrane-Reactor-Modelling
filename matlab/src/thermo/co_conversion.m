function X = co_conversion(F0, F)
%CO_CONVERSION  CO conversion  X = 1 - F_CO / F0_CO  (flows in mol/s).

c = pbmr_const();
X = 1 - F(c.idx.CO) / F0(c.idx.CO);
end
