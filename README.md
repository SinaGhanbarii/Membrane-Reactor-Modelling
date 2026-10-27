# pbmr1d

1D packed-bed Pd-membrane reactor model for water-gas shift (PhD, Aston University / EBRI).

## Status
- Phase 0 (setup): done
- **Phase 1 (thermodynamic bounds): done** - `pbmr1d/thermo.py`, `scripts/phase1_equilibrium.py`
- Phase 2 (isothermal reactor): next

## Quick start
    conda env create -f environment.yml
    conda activate pbmr1d
    python -m pytest -q                         # 22 tests
    python scripts/phase1_equilibrium.py         # figures + summary in results/phase1/

The scripts read `params/rig_private.yaml` (confidential H2Site data, git-ignored) and fall
back to `params/rig_example.yaml` if it is missing. Keep all rig values in the private file.

## Conventions
SI units internally (mol/s, K, Pa, J/mol). Species order: CO, H2O, CO2, H2, CH4, N2.
Ideal-gas thermodynamics from GRI-Mech 3.0 via Cantera; graphite from nasa_condensed.yaml.
