# Membrane reactor modelling (pbmr1d)

1D model of a packed-bed water-gas shift (WGS) reactor with a Pd membrane, for the PhD project
*Membrane and Fuel Cell Integration for Low Carbon Hydrogen Utilisation* (EBRI, Aston University).
Last updated 18 November 2026.

**Status: the model is implemented and verified against itself. It is not yet validated against measurements.**

## Start here

1. **Open [`matlab/reactor_walkthrough.mlx`](matlab/reactor_walkthrough.mlx) in MATLAB R2024b or later.**
   It is a live script: 14 short sections with text, code and saved figures, so the results can be read without running anything.
   GitHub cannot display `.mlx` files; download the file or clone the repository. The same content as plain text is in
   [`matlab/reactor_walkthrough.m`](matlab/reactor_walkthrough.m), which GitHub does display.
2. **Read the design document** [`docs/1D_Pd-membrane_WGS_reactor_model_design_document_v3.docx`](docs/1D_Pd-membrane_WGS_reactor_model_design_document_v3.docx):
   question, hypotheses, model equations, parameter sources, verification, risks, and the decisions I need from you (section 8, D1 to D12).

## What the walkthrough covers

| Section | Content |
|---|---|
| 1-3 | Setup, how the model was designed, code architecture |
| 4-7 | Where every parameter came from: reported, inferred, assumed or missing, with sources and validity windows |
| 8-9 | Base case (40 % CO, 400 C, 4.5 bar, membrane M7) and an editable case of your own |
| 10 | Which mechanism limits the bed, and where (five regimes: no/reverse flux, reaction, pinch, permeance, idle) |
| 11-12 | Co-current against counter-current, sweep gas, and how recovery falls as the membrane flux is reduced |
| 13-14 | Source of any single number, and the main limits of the model |

## Main early results (all unvalidated)

- **Temperature.** Without wall cooling the bed peaks at about 510 C (350 C inlet) and 550 C (400 C inlet). Both exceed the 450 C upper limit of the M7 flux data. Only strong cooling keeps the peak below it, and the wall coefficient is assumed.
- **Rate laws.** K1 (Keiski 1996) and G2 (Hla 2009) agree at 400-450 C and differ at 350 C, which is outside the fitted range of G2.
- **Counter-current sweep.** A hydrogen recovery near 100 % is an ideal-membrane upper bound (pure-H2 permeate, inert sweep, no polarisation). In the base counter-current case the reaction, not the membrane, is the limit, and about half of the bed is idle. Recovery falls to 0.95 at one third of the flux, and to 0.30 at one twentieth.

## Main limits and open issues

- No comparison with measured data yet. MATLAB and an independent Python implementation agree to about six digits, but both use the same equations and parameter register, so this confirms the implementation, not the physics.
- The M7 flux data cover 1.2-1.6 bar; the model is used at 4-5 bar.
- **K1 rate units are unresolved.** Keiski prints concentrations in mol dm-3 and rates in mol g-1 s-1, but read literally these give a rate far too high, and the paper's own equation mixes hours, seconds, dm3 and cm3. The plausible readings differ by one constant factor on the whole rate (about 1 to 22 relative to the current setting). A sensitivity test is planned; until then, rate-limited results depend on this choice.
- Wall heat transfer, bed voidage, effective diffusivity and shell diameter are assumed.
- Steady state only. Concentration polarisation and any H2O or CO2 effect on flux are not included.

## Running the checks

```matlab
>> cd('<repo>\matlab')
>> setup_paths
>> run_tests_step3          % 68 checks
>> run_verification         % 385 checks over varied and unusual conditions; writes results/verification_report.csv
```

Each prints a pass/fail summary. See [`matlab/README.md`](matlab/README.md) for the code layout.

## Repository map

| Folder | Content |
|---|---|
| `matlab/` | The MATLAB model, the walkthrough, tests and the parameter register (`matlab/data/params_register.csv`) |
| `docs/` | Design document, research log, review of the membrane literature |
| `pbmr1d/`, `tests/`, `scripts/` | Python part: thermodynamic bounds (Phase 1, 22 tests) and the reference implementation used to cross-check MATLAB |
| `Phase 2 Papers/` | Open-literature PDFs behind the parameters, for reading only (publisher copyright) |
| `params/` | Example rig parameters only |

The Python part needs `conda env create -f environment.yml`, then `python -m pytest -q`.

## Data and confidentiality

The MATLAB model uses open-literature parameters only. No H2Site rig data enter the model. Files named `params/*_private.yaml`
and the `results/` folder (generated output) are excluded from Git. Do not share the repository beyond the supervisory team.

## Conventions

SI units internally (K, Pa, mol/s, J/mol). Species order: CO, H2O, CO2, H2, CH4, N2. Reaction: CO + H2O <=> CO2 + H2.
