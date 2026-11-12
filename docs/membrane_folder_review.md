# Review: "H2Site - More Details About Membrane" folder (6 Nov 2026)

Scope: 14 PDFs, all public literature or patents (TU/e - Tecnalia - Hydrogen Onsite group).
None is an H2Site operating document, so they are citable. Status: read for model-relevant
content; figures were not digitised; PDFs not re-verified page by page beyond what is quoted.

## 1. Inventory

| File | What it is | Use for the model |
|---|---|---|
| 1-s2.0-S036031991402374X | Fernandez et al. (2015) Int. J. Hydrogen Energy 40, 3506-3519. Thin Pd-Ag membranes, WGS-related gases, fluidised-bed tests | Catalyst-membrane interaction (Ti-containing supports cut flux; alumina/zirconia did not); external mass-transfer limits |
| 1-s2.0-S0360319915025914 | Medrano et al. (2016) Int. J. Hydrogen Energy 41, 8706-8718. Metallic-supported Pd-Ag, >800 h, reforming | Membrane stability; not WGS |
| 1-s2.0-S0376738816318191 | Melendez et al. (2017) J. Membr. Sci. 528, 12-23. Ultra-thin (0.46-1.29 um) Pd-Ag on ceramic | Q, Ea, n, CO inhibition, N2/CH4/CO2 effects |
| 1-s2.0-S0378382021000515 | Cechetto et al. (2021) Fuel Process. Technol. 216, 106772. NH3 decomposition MR | Membrane reactor operation with vacuum permeate; not WGS |
| 1-s2.0-S1383586607004790 | Barbieri et al. (2008) Sep. Purif. Technol. 61, 217-224. Sieverts-Langmuir model for CO | CO inhibition model and parameters |
| 20170327_CO_Helmi_Siasi_Farimani | Helmi (2017) PhD thesis, TU/e. FBMR for WGS | Concentration polarisation, Pm0/Ea, FBMR hydrodynamics |
| 387047507446037 | Helmi et al. (2013) ECCE9 abstract/poster. Packed vs fluidised MR for WGS | Qualitative only |
| arradeve2018 | Arratibel et al. (2018) J. Membr. Sci. 550, 536-544. Double-skin membranes | Membrane durability; mass-transfer limits at high permeance |
| catalysts-13-00920 | Cechetto et al. (2023) Catalysts 13, 920. Metallic-supported Pd-Ag, NH3 | Not WGS |
| Modelling_of_membrane_reactors... | Gallucci, van Sint Annaland, Kuipers, Chapter 10 (book details not on the pages read) | Model equations (see section 4) |
| molecules-21-00376 | Helmi et al. (2016) Molecules 21, 376. FBMR WGS demonstration, ~900 h | Closest WGS MR dataset; K_eq source; methanation note |
| processes-12-00040 | Agnolin et al. (2024) Processes 12, 40. Pd-Ag on Hastelloy-X, SMR | Not WGS |
| US11395988 | US 11,395,988 B2 (granted 26 Jul 2022), assignee Hydrogen Onsite S.L. | Claims a double-skin membrane (Pd95-Ag5, YSZ/gamma-Al2O3 protective layer). Membrane, not reactor design |
| WO2018104455A3 | PCT/EP2017/081835, applicants TU/e and Tecnalia; this copy is the front page and search report only (scanned) | Same family as the US patent |

## 2. Does this describe the 555.003 reactor?  NOT ESTABLISHED
- The TU/e/Tecnalia WGS work is a FLUIDISED-bed membrane reactor (bubbling bed, 2 wt% Pt/Al2O3
  catalyst, 1-2.5 bar, CO 5-10 %, 400-450 C, permeate ~30 mbar). Operating pressure and CO
  content are far below the 555.003 design window.
- The H2Site manual speaks of a "catalyst bed" and "wet catalyst"; nothing read says fluidised,
  and the photograph shows only the insulated vessel. Treat the reactor type as UNKNOWN until
  H2Site/supervisor confirm.
- The link between H2Site and Hydrogen Onsite S.L. (patent assignee, authors' affiliation) is my
  inference only (both Bizkaia). Not confirmed.
- Consequence: membrane parameters may transfer; fluidised-bed hydrodynamics (bubble/emulsion,
  Kbe) do NOT transfer to a fixed bed.

## 3. Parameters extracted (membrane)
All at 400 C unless stated. Flux J = (Pe/delta)(p_ret^n - p_perm^n), Pe = Pe0 exp(-Ea/RT).

| Source | Membrane | Values | Notes |
|---|---|---|---|
| Melendez et al. (2017), Tables 6-7, Fig. 7 | Pd-Ag, 0.46 / 0.78 / 1.17 / 1.29 um | Permeability 5.52 / 6.31 / 8.19 / 9.03 e-9 mol s-1 m-1 Pa-0.5 (n=0.5 basis, 100 kPa). Ea 5.31 / 6.40 / 7.58 / 7.44 kJ/mol. Best-fit n 0.58-0.81 | n > 0.5, closer to 1 for thinnest layer; authors attribute to support resistance/composition; thinnest membrane H2/N2 only ~50 |
| Helmi (2017) thesis, Ch. 4 | Pd-Ag, ~4.5 um | Pm0 = 1.76e-8 mol m-1 s-1 Pa-0.5, Ea = 7.1 kJ/mol (n=0.5) | Pm0 is the pre-exponential; my calc gives ~4.95e-9 at 400 C |
| Helmi et al. (2016) | Pd-Ag, five-membrane module | Permeance 3.89e-6 mol m-2 s-1 Pa-1 at 400 C, 1 atm difference (abstract) | Thickness not taken from the text read |
| Barbieri et al. (2008) | 60 um Pd-Ag, 374 C, <= 700 kPa | alpha = 0.2319, K_CO = 0.023 kPa-1 (assumed T-independent alpha) | Lumped Sieverts-Langmuir: J x [1 - alpha K p_CO/(1+K p_CO)] |
| Melendez et al. (2017), Table 9 | Pd-Ag 0.78 um | K_CO 0.009-0.051 kPa-1; alpha 1.11-3.47 | alpha > 1 is physically odd for a coverage-scaling factor; authors say model incomplete |

Cross-check (my calculation, p_H2,ret = 1 bar, permeate 0, clean H2, thickness-scaled to 1.29 um):
Melendez 2.21, Helmi thesis 1.21, Boon (2015) 1.48 mol m-2 s-1. Same order of magnitude; n differences
make this indicative only.

CO inhibition is the weakest link. At p_CO ~ 84 kPa (my earlier rig-type estimate) the flux factor
is 0.85 with the Barbieri pair but 0.26 with the Melendez 60 % H2 / 15 % CO pair (illustrative
arithmetic, not a model result). Melendez's own measured loss was 13-15 % at 15 % CO, so their (alpha, K)
only make sense at much lower p_CO than those tested at rig pressure. Carry this as a sensitivity
parameter, not a fixed value.

Also from Melendez: H2/N2 and H2/CH4 fluxes were virtually identical; H2/CO2 slightly lower
(attributed to reverse WGS). External mass transfer prevailed for mixtures because the membranes
are so permeable (also Arratibel 2018; Fernandez 2015) -> concentration polarisation must be modelled.

## 4. Provenance fix and model structure
- K_eq = exp(4577.8/T - 4.33) is Moe (1962), Chem. Eng. Prog. 58, 33-36, quoted as ref. [3] in Helmi
  et al. (2016). Cite Moe (1962) via Helmi et al. (2016); drop the Karagoz et al. (2020) attribution.
  (Moe 1962 itself was not read.)
- Gallucci et al. Chapter 10: Sec. 10.2 defines limit conversion (pseudo-equilibrium) exactly as in
  Phase 1; Sec. 10.3.1 gives the 1D pseudo-homogeneous PBMR model WITH time-derivative terms
  (eqs 10.4-10.12: continuity, momentum with Ergun-type friction, species with axial dispersion,
  energy with gas+solid heat capacity), Sieverts and Sieverts-Langmuir flux (eqs 10.9-10.10).
  This is the template for the dynamic distributed model; Sec. 10.3.2 covers 2D models.
- Methanation: Helmi et al. (2016) report the (Pt) catalyst active for methanation only at >= ~450 C
  and recommend ~400 C. Catalyst-specific; says nothing about the 555.003 catalyst.

## 5. Still missing (not in this folder)
- WGS kinetics for the catalyst H2Site uses (Helmi's reactor ran near equilibrium, so no rate law).
- Packed-bed concentration-polarisation / wall heat-transfer correlations beyond Chapter 10.
- Numerical WGS membrane-reactor data: Helmi et al. (2016) results appear to be in figures only
  (Tables 1-5 are properties/operating window); digitising needed, not yet done.
- CO inhibition at 25-45 % CO and 4-5 bar.
