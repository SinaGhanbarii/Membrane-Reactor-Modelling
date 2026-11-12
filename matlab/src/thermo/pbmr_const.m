function c = pbmr_const()
%PBMR_CONST  Constants and fixed conventions shared by all pbmr1d functions.
%   Species order is fixed:  1 CO, 2 H2O, 3 CO2, 4 H2, 5 CH4, 6 N2.
%   SI units throughout (K, Pa, mol/s, J/mol).
%   WGS:  CO + H2O <=> CO2 + H2   (delta n = 0)

c.R  = 8.31446261815324;          % J/(mol K), same value as Cantera
c.P0 = 101325;                    % Pa, standard-state pressure of the NASA data

c.species = {'CO', 'H2O', 'CO2', 'H2', 'CH4', 'N2'};
c.idx = struct('CO', 1, 'H2O', 2, 'CO2', 3, 'H2', 4, 'CH4', 5, 'N2', 6);

c.NU = [-1; -1; 1; 1; 0; 0];      % WGS stoichiometric coefficients

% Element matrix: rows C, H, O, N; columns = species order above
c.ATOMS = [1 0 1 0 1 0;
           0 2 0 2 4 0;
           1 1 2 0 0 0;
           0 0 0 0 0 2];
end
