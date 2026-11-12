"""Export thermodynamic data and Python reference values for the MATLAB port (Step 1).

Writes (default: ./matlab_export/):
    data/nasa7_species.json     NASA-7 coefficients (GRI-Mech 3.0) for the six species
    tests/vectors_step1.json    reference values computed with pbmr1d.thermo (Cantera)

Before writing anything it checks that an independent NASA-7 evaluator (the same
formulas the MATLAB code uses) reproduces Cantera to near machine precision.
Run:  python scripts/export_matlab_reference.py [outdir]
"""
import json
import sys
from pathlib import Path

import numpy as np
import cantera as ct

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repo root
from pbmr1d import thermo as t

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("matlab_export")
R = ct.gas_constant / 1000.0                       # J/(mol K), same constant as Cantera

# ---------------------------------------------------------------- NASA-7 data
all_sp = {s.name: s for s in ct.Species.list_from_file("gri30.yaml")}
rows = []
for name in t.SPECIES:
    th = all_sp[name].thermo
    c = np.asarray(th.coeffs, dtype=float)         # [T_mid, high(7), low(7)] (checked below)
    rows.append(dict(name=name, T_low=th.min_temp, T_mid=c[0], T_high=th.max_temp,
                     a_high=c[1:8].tolist(), a_low=c[8:15].tolist()))


def nasa_eval(row, T):
    a = np.array(row["a_low"] if T < row["T_mid"] else row["a_high"])
    cp = R * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)
    h = R * T * (a[0] + a[1]*T/2 + a[2]*T**2/3 + a[3]*T**3/4 + a[4]*T**4/5 + a[5]/T)
    s = R * (a[0]*np.log(T) + a[1]*T + a[2]*T**2/2 + a[3]*T**3/3 + a[4]*T**4/4 + a[6])
    return cp, h, s, h - T*s


T_LIST = [523.15, 573.15, 623.15, 673.15, 723.15, 773.15]
worst = 0.0
for T in list(T_LIST) + [300.0, 999.0, 1001.0, 1500.0]:       # incl. both sides of T_mid
    cp_c, h_c, g_c = t.species_props(T)
    for i, row in enumerate(rows):
        cp, h, s, g = nasa_eval(row, T)
        for mine, ref in ((cp, cp_c[i]), (h, h_c[i]), (g, g_c[i])):
            worst = max(worst, abs(mine - ref) / max(abs(ref), 1.0))
assert worst < 1e-10, f"NASA evaluator disagrees with Cantera (worst rel err {worst:.2e})"
print(f"NASA-7 evaluator vs Cantera: worst relative difference {worst:.1e}  (OK)")

# ----------------------------------------------------------------- write data
(OUT / "data").mkdir(parents=True, exist_ok=True)
(OUT / "tests").mkdir(parents=True, exist_ok=True)
species = dict(
    note=("NASA-7 polynomials from GRI-Mech 3.0 (gri30.yaml via Cantera "
          f"{ct.__version__}); standard state 1 atm; SI units after multiplying by R."),
    R=R,
    P0=ct.one_atm,
    names=t.SPECIES,
    mw=t.MW.tolist(),                              # kg/mol
    T_low=[r["T_low"] for r in rows],
    T_mid=[r["T_mid"] for r in rows],
    T_high=[r["T_high"] for r in rows],
    a_low=[r["a_low"] for r in rows],
    a_high=[r["a_high"] for r in rows],
)
(OUT / "data" / "nasa7_species.json").write_text(json.dumps(species, indent=1))

# ------------------------------------------------------------ reference values
F0 = np.array([0.40, 1.20, 0.15, 0.35, 0.05, 0.05]) * 1e-3     # mol/s, same as test_thermo.py
F_steam = F0.copy()
F_steam[t.IDX["H2O"]] = 1e4 * F_steam[t.IDX["CO"]]

vec = dict(
    note="Reference values from pbmr1d.thermo (Python/Cantera) for MATLAB Step 1 tests.",
    T=T_LIST,
    cp=[t.species_props(T)[0].tolist() for T in T_LIST],        # rows: T, cols: species
    h=[t.species_props(T)[1].tolist() for T in T_LIST],
    g0=[t.species_props(T)[2].tolist() for T in T_LIST],
    K_nasa=[t.K_wgs(T) for T in T_LIST],
    K_moe=[t.K_wgs_moe(T) for T in T_LIST],
    dH=[t.dH_wgs(T) for T in T_LIST],
    F0=F0.tolist(),
    xi=[t.wgs_extent_eq(F0, T) for T in T_LIST],
    F0_steam=F_steam.tolist(),
    X_steam_673=t.co_conversion(F_steam, t.wgs_equilibrium(F_steam, 673.15)),
)
(OUT / "tests" / "vectors_step1.json").write_text(json.dumps(vec, indent=1))
print("Wrote", OUT / "data" / "nasa7_species.json", "and", OUT / "tests" / "vectors_step1.json")
print("K_nasa:", [f"{k:.5f}" for k in vec["K_nasa"]])
print("xi/CO_in:", [f"{x / F0[0]:.5f}" for x in vec["xi"]])
