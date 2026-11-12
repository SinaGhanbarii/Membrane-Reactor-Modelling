# pbmr1d - MATLAB port (built step by step)

MATLAB implementation of the 1D Pd-membrane WGS reactor model. Python (`pbmr1d/`) is
the verified reference; every MATLAB step is checked against reference values exported
from Python (`scripts/export_matlab_reference.py`). Base MATLAB only - no toolboxes.

## Run
```matlab
>> cd('<repo>\matlab')
>> setup_paths
>> run_tests_step1
```
Expected last line: `ALL STEP 1 TESTS PASSED`. If anything prints `[FAIL]`, copy the whole
output to Claude.

## Layout
```
matlab/
  setup_paths.m             adds src/ and tests/ to the path
  data/nasa7_species.json   NASA-7 coefficients, GRI-Mech 3.0 (from Cantera)
  src/thermo/               Step 1: cp, h, s, g0, K_eq, dH, equilibrium extent
  tests/                    vectors_step1.json + run_tests_step1.m
```

## Conventions
SI units (K, Pa, mol/s, J/mol). Species order: 1 CO, 2 H2O, 3 CO2, 4 H2, 5 CH4, 6 N2.
WGS: CO + H2O <=> CO2 + H2. Flow vectors are 6x1 columns.

## Provenance
- Thermodynamic data: GRI-Mech 3.0 NASA-7 polynomials via Cantera (version in the JSON note).
- `Keq_moe`: Moe (1962) as quoted by Helmi et al. (2016); cross-check only.

## Steps
1. Thermodynamic core (this step).
2. Parameters (JSON with `source:` fields), rate law, membrane flux (planned).
3. Steady 1D reactor, verification tests, benchmarks (planned).
