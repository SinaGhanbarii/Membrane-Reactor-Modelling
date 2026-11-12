# Research log - pbmr1d

## 2026-10-23 (Fri) - Phase 0 started
- Week 1 of 1D modelling plan. Environment and repository set up.

## 2026-10-27 (Tue) - Phase 1 complete (thermodynamic bounds)
**Done**
- `thermo.py`: Cantera properties, exact WGS equilibrium, membrane-reactor (MR) limit,
  adiabatic lines, full Gibbs equilibrium, carbon activity. 22 tests pass.
- Figures 1-5 and `results/phase1/summary.md`.

**Findings (to discuss with supervisor)**
1. Nominal feed: S/CO = 3.73. Conventional equilibrium X_CO = 94 % at 400 C, 91 % at 450 C.
   The membrane can add at most ~5-8 percentage points of CO conversion here.
2. MR limit at 0.15-0.2 bara permeate: H2 recovery ceiling ~89-93 %. H2Site's nominal
   permeate (0.02 kg/h) is ~33 % recovery, i.e. ~35-37 % of the thermodynamic ceiling.
   => at nominal conditions the rig is limited by membrane area/kinetics, not thermodynamics.
3. Excess steam lowers the recovery ceiling (retentate dilution). Recovery ceiling peaks
   near S/CO ~1.3-1.5; trade-off against carbon and catalyst over-reduction risk.
4. Methanation is strongly favoured thermodynamically (CH4 x5 at 250 C, x4 at 450 C, full
   Gibbs). "CH4 inert" is a catalyst-selectivity assumption -> check CH4 in = CH4 out by GC.
5. Graphite activity >> 1 below ~550 C even at WGS equilibrium. Absence of carbon depends on
   kinetics, not thermodynamics.
6. Adiabatic: T_in 250 C -> adiabatic equilibrium at ~416 C, 93 %. Delta T_ad(full) ~178 K.
7. Document inconsistencies: (a) 3.6 SLPM != 0.02 kg/h at any common reference T
   (3.71 SLPM at 0 C); (b) for the "original reactor" feed, 0.02 kg/h implies 102 % H2
   recovery - impossible -> Table 1 permeate is a nameplate value, not a prediction;
   (c) original composition sums to 99.85 %.

**Open questions**
- Which permeate pressure (0.15 or 0.2 bara) and inlet T (202/250/350 C) are real?
- Catalyst identity (selectivity against methanation is now a key question).
- K_eq: Cantera vs Moe-form differ 3-8 %; check against Twigg/other correlation when read.

**Next**: Phase 2 - isothermal reactor model (needs placeholder kinetics + membrane params).

## 6 Nov 2026 - Review of "H2Site - More Details About Membrane" folder
- 14 public papers/patents reviewed; note in docs/membrane_folder_review.md.
- K_eq correlation traced to Moe (1962) via Helmi et al. (2016); Karagoz attribution to be dropped.
- Membrane parameters extracted (Melendez 2017, Helmi thesis, Barbieri 2008); CO inhibition parameter sets disagree strongly at rig-like p_CO -> treat as sensitivity.
- Open: is the 555.003 reactor fixed or fluidised bed? Which catalyst? (decides what transfers).
