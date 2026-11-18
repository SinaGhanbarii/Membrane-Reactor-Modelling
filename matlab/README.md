# MATLAB model of the 1D Pd-membrane WGS reactor

Base MATLAB (R2024b), no toolboxes. This folder holds the whole model. The Python code in `../pbmr1d/` and `../scripts/` is the independent reference used to cross-check it.

## Where to start

Open **`reactor_walkthrough.mlx`** (Live Editor). It runs the model section by section and shows the saved results.
`reactor_walkthrough.m` is the same script as plain text.

```matlab
>> cd('<repo>\matlab')
>> setup_paths
>> open reactor_walkthrough.mlx
```

Sections 10 to 12 use the case defined in section 9, so run section 9 first. Section 12 makes about 22 model runs and takes over a minute.

## Tests

| Script | Checks | What it covers |
|---|---|---|
| `tests/run_tests_step1.m` | - | Thermodynamics against Python reference values |
| `tests/run_tests_step2.m`, `run_tests_step2b.m` | 22 each | Parameter register, rate laws, membrane flux |
| `tests/run_tests_step3.m` | 68 | Reactor: conservation, limits worked out by hand, regimes, input checks |
| `tests/run_verification.m` | 385 | 20 reference cases against Python, 54 varied runs, 17 trends, 21 unusual conditions, 30 bad inputs |

All tests print a pass/fail summary. The reference cases agree with the Python model to about six digits. Both use the same equations and parameter register, so agreement shows a correct implementation, not correct physics.

## Layout

```
matlab/
  reactor_walkthrough.mlx / .m   guided tour of the model and its results
  setup_paths.m                  adds src/ and tests/ to the path
  data/                          params_register.csv (every parameter with value, unit, source, kind), nasa7_species.json
  src/thermo/                    heat capacity, equilibrium constant, WGS equilibrium, atom balances
  src/kinetics/                  rate laws K1 (Keiski 1996) and G2 (Hla 2009), pore-diffusion factor
  src/membrane/                  Sieverts flux with CO inhibition, validity checks
  src/params/                    register loader, validity windows, evidence and relation plots
  src/reactor/                   make_case, run_case, reactor_rhs, co-current and counter-current solvers, regime labels
  src/plots/                     profile, flow comparison, regime map, architecture figure
  src/studies/                   regime, sweep and hypothesis studies (write to ../results/, which Git ignores)
  tests/                         the tests above
```

## Model in brief

Steady plug flow in an annulus between a Pd tube and the shell, with catalyst in the annulus. Balances for six species, H2 in the permeate stream, temperature and pressure. Flux follows Sieverts' law, with an optional inert sweep gas on the permeate side. Counter-current flow is solved by shooting. Reaction rate is the literature rate law multiplied by (1 - beta), the distance from equilibrium.

## Conventions

SI units (K, Pa, mol/s, J/mol). Species order: 1 CO, 2 H2O, 3 CO2, 4 H2, 5 CH4, 6 N2. Flow vectors are 6x1 columns.

## Provenance

Thermodynamic data: GRI-Mech 3.0 NASA-7 polynomials via Cantera. Equilibrium constant cross-check: Moe (1962) as quoted by Helmi et al. (2016). Every other parameter is in `data/params_register.csv` with its source and a kind: reported, inferred, assumed, placeholder, sweep or range.

## Known limits

Nothing here is validated against rig data. See section 14 of the walkthrough and the design document.
