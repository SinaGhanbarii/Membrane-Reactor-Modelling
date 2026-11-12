# Reading notes: Keiski (1996), Adrover (2009), Gallucci (2011) Ch.10 §10.2–10.3.1

*Draft for Sina, 12 Nov 2026. Rewrite in your own words before it goes anywhere. Every statement is from the text of the three sources; items marked **[inference]** are mine and not in the papers.*

## 1. Keiski et al. (1996), Appl. Catal. A 137, 349–370: HT-WGS kinetics, Fe3O4–Cr2O3

- **Model used for fitting.** Stationary kinetics from a pseudo-homogeneous 1D plug-flow tube: dF_CO/dV = −r_CO·ρ_B (Eq. 30).
- **Rate law (Eq. 31).** r_CO = k0·exp(−E/RT)·c_CO^n·c_H2O^m·c_CO2^p·c_H2^q·(1−β), with **concentrations in mol dm⁻³** and β the reversibility factor.
- **Best fit.** Four parameters: ln k0 = 26.1, E = 95.0 kJ/mol (E/R ≈ 11 430 K), n = 1.1, m = 0.53. The CO2 and H2 orders were not retained (Tables 4–5, p and q ≈ 0).
- **Catalyst age matters.** Table 4 (200–280 h) and Table 5 (725–762 h) give different parameter sets; the paper says some deactivation occurs over the first 725 h. The best-fit set comes from the *younger* catalyst.
- **Validity box (my reading of Tables 1–3 and the abstract).** About 573–633 K, CO 3–10 %, tube runs at 3 and 9 bar. Above that is extrapolation. The abstract's "3–5 bar" disagrees with the tables.
- **Unit problem.** The notation gives c in mol dm⁻³, ρ_B in g cm⁻³, F_CO in mol h⁻¹ and V in dm³. It gives **no unit for k0**. **[inference]** The rate is per kg of catalyst. A plausibility estimate points to seconds, not hours, but this must be tested by reproducing one of Keiski's own experiments.
- **Reuse elsewhere.** Barbieri (2011) and Caravella (2016) reuse these numbers with partial pressures and no stated conversion, so they are not independent evidence.

## 2. Adrover et al. (2009), AIChE J. 55, 3206–3213: non-isothermal WGS membrane reactor

- **Model.** 1D pseudo-homogeneous, **steady state**, concentric-tube geometry; catalyst inside the Pd tubes, sweep gas in an adiabatic shell, co- and counter-current.
- **Stated assumptions.**
  - (a) Mass and heat dispersion neglected.
  - (b) Negligible catalyst–gas mass and heat resistances.
  - (c) Isobaric on the permeate side.
  - (d) Infinite H2 selectivity.
  - (e) Adiabatic shell.
- **Kinetics.** Langmuir–Hinshelwood (Podolski & Kim); the review found its constants are not stated in the paper.
- **Heat transfer.** Overall U from Abo-Ghander et al.; wall coefficient and radial conductivity from Dixon & Cresswell; shell side from Kern.
- **Result for my work.** With counter-current sweep, heat feedback gives **steady-state multiplicity** (S-shaped curves; up to three steady states at the same inlet temperature); co-current does not. The authors say a 2D model would be needed to examine radial gradients.
- **Limit for my work.** It is a steady-state study. It says nothing about dynamics, and its low-diameter tubes (0.3 in) support the 1D assumption.

## 3. Gallucci, van Sint Annaland, Kuipers (2011), Ch.10 §10.2 and §10.3.1.1

- **Limit conversion (§10.2).**
  - The membrane reactor's pseudo-equilibrium needs two conditions at once: chemical equilibrium in the reaction zone, and equality of the permeating species' partial pressures across the membrane (Eq. 10.3).
  - With zero permeate pressure, or an infinite inert sweep, the limit is 100 %. With a finite permeate pressure it can be well below 100 %.
  - Beyond the limit, more membrane area or thinner membranes do not help. Non-isothermal cases also need the energy balance.
- **1D pseudo-homogeneous model (§10.3.1.1).**
  - It is the most used model for lab-scale packed-bed membrane reactors, and describes radially averaged axial profiles.
  - Eqs. 10.4–10.12 contain time-derivative terms: continuity, momentum with an Ergun-type friction term, component balances with axial dispersion and a membrane flux term, and an energy balance with gas and solid heat capacity.
  - The steady model drops the time derivatives. Danckwerts boundary conditions are the usual choice.
  - The flux is Richardson/Sieverts (Eq. 10.9) or Sieverts–Langmuir for CO inhibition (Eq. 10.10).
- **When 1D holds.** The chapter says 1D works when flux through the membrane limits performance (thick membranes of 40–100 µm). With ultra-thin membranes, bed-to-wall mass transfer (concentration polarisation) becomes prominent, the hypotheses behind the 1D model break down, and a 2D model is needed. The same applies to the flux: a 1D model overestimates it if the wall partial pressure is lower than the bulk.

## 4. What this means for question A

- **Kinetic input.** Keiski is the only independent primary rate law, and its unit basis must be resolved first.
- **Thermal behaviour.** Adrover shows it can change the steady state, so the energy balance belongs in the 1D model.
- **Concentration polarisation.** For thin membranes the 1D assumption itself is at risk. The model needs a bed-to-wall mass-transfer term, and a lumped film coefficient is the simplest option.
- **Limit conversion.** Computing it first (§10.2) is the check for every later result.
