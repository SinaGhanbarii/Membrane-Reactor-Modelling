"""Thermodynamics for the WGS packed-bed membrane reactor (Phase 1).

Conventions
-----------
* SI units throughout: T in K, P in Pa, molar flows in mol/s, energies in J/mol.
* Species order is fixed by SPECIES; flow vectors are numpy arrays in that order.
* All gas-phase properties are ideal-gas properties from GRI-Mech 3.0 (NASA-7
  polynomials) via Cantera. Graphite comes from Cantera's nasa_condensed.yaml.

Reaction considered:  CO + H2O <=> CO2 + H2   (WGS, delta_n = 0)
"""
from __future__ import annotations

import numpy as np
import cantera as ct
from scipy.optimize import brentq

R = ct.gas_constant / 1000.0          # J/(mol K)
P_REF = ct.one_atm                    # standard-state pressure of the GRI data (Pa)

SPECIES = ["CO", "H2O", "CO2", "H2", "CH4", "N2"]
IDX = {s: i for i, s in enumerate(SPECIES)}
NU_WGS = np.array([-1.0, -1.0, +1.0, +1.0, 0.0, 0.0])

# Element matrix (rows: C, H, O, N; columns: SPECIES)
ELEMENTS = ["C", "H", "O", "N"]
ATOMS = np.array([
    # CO H2O CO2 H2 CH4 N2
    [1,  0,  1,  0,  1,  0],   # C
    [0,  2,  0,  2,  4,  0],   # H
    [1,  1,  2,  0,  0,  0],   # O
    [0,  0,  0,  0,  0,  2],   # N
], dtype=float)


def _gas_subset(names):
    """Ideal-gas Cantera phase containing only the named GRI-Mech species."""
    all_sp = {s.name: s for s in ct.Species.list_from_file("gri30.yaml")}
    return ct.Solution(thermo="ideal-gas", species=[all_sp[n] for n in names])


_GAS = _gas_subset(SPECIES)
MW = _GAS.molecular_weights / 1000.0  # kg/mol, in SPECIES order


# ---------------------------------------------------------------------------
# Pure-species ideal-gas properties
# ---------------------------------------------------------------------------
def species_props(T):
    """Return (cp, h, g0) arrays [J/mol/K, J/mol, J/mol] at T for SPECIES.

    h is the absolute (formation-based) enthalpy; g0 is the standard Gibbs
    energy at P_REF.
    """
    _GAS.TP = T, P_REF
    RT = R * T
    cp = _GAS.standard_cp_R * R
    h = _GAS.standard_enthalpies_RT * RT
    g0 = _GAS.standard_gibbs_RT * RT
    return cp, h, g0


def dH_wgs(T):
    """Reaction enthalpy of WGS at T [J/mol] (negative = exothermic)."""
    return float(NU_WGS @ species_props(T)[1])


def dG0_wgs(T):
    """Standard Gibbs energy of WGS at T [J/mol]."""
    return float(NU_WGS @ species_props(T)[2])


def K_wgs(T):
    """Equilibrium constant of WGS from Cantera data (dimensionless; delta_n = 0)."""
    return float(np.exp(-dG0_wgs(T) / (R * T)))


def K_wgs_moe(T):
    """Empirical K_eq = exp(4577.8/T - 4.33), as printed in Karagoz et al. (2020).

    Used only as an independent cross-check of the Cantera value.
    """
    return float(np.exp(4577.8 / T - 4.33))


def mixture_enthalpy_flow(F, T):
    """Total enthalpy flow sum(F_i h_i(T)) [W] for flows F [mol/s]."""
    return float(np.asarray(F) @ species_props(T)[1])


# ---------------------------------------------------------------------------
# Feed construction
# ---------------------------------------------------------------------------
def feed_flows(syngas_kg_h, composition_mol_pct, water_kg_h):
    """Molar feed flows [mol/s] in SPECIES order from rig-style inputs.

    composition_mol_pct is normalised to 100 % (the sum is returned so callers
    can report documents whose composition does not close).
    """
    x = np.zeros(len(SPECIES))
    for name, pct in composition_mol_pct.items():
        x[IDX[name]] = pct
    total_pct = x.sum()
    x = x / total_pct
    mw_mix = float(x @ MW)                        # kg/mol
    n_syngas = syngas_kg_h / 3600.0 / mw_mix      # mol/s
    F = x * n_syngas
    F[IDX["H2O"]] += water_kg_h / 3600.0 / MW[IDX["H2O"]]
    return F, {"sum_pct": total_pct, "mw_syngas_g_mol": mw_mix * 1000.0}


def element_flows(F):
    """Element flows [mol/s] (C, H, O, N) for species flows F."""
    return ATOMS @ np.asarray(F)


# ---------------------------------------------------------------------------
# Conventional WGS equilibrium (exact, extent-based)
# ---------------------------------------------------------------------------
def wgs_extent_eq(F0, T):
    """Equilibrium extent xi [mol/s] of WGS alone, CH4 and N2 inert.

    Solves K (a - xi)(d - xi) = (b + xi)(c + xi), with a=CO, d=H2O, b=CO2,
    c=H2. Independent of total pressure because delta_n = 0 (ideal gas).
    """
    a, d = F0[IDX["CO"]], F0[IDX["H2O"]]
    b, c = F0[IDX["CO2"]], F0[IDX["H2"]]
    K = K_wgs(T)
    A = K - 1.0
    B = -(K * (a + d) + (b + c))
    C = K * a * d - b * c
    lo, hi = -min(b, c), min(a, d)
    if abs(A) < 1e-12:
        roots = [-C / B]
    else:
        disc = np.sqrt(B * B - 4 * A * C)
        roots = [(-B - disc) / (2 * A), (-B + disc) / (2 * A)]
    valid = [r for r in roots if lo - 1e-14 <= r <= hi + 1e-14]
    if len(valid) != 1:
        raise RuntimeError(f"Expected one physical root, got {roots} in [{lo}, {hi}]")
    return float(valid[0])


def wgs_equilibrium(F0, T):
    """Outlet flows [mol/s] at WGS-only equilibrium."""
    return np.asarray(F0) + NU_WGS * wgs_extent_eq(F0, T)


def co_conversion(F0, F):
    return 1.0 - F[IDX["CO"]] / F0[IDX["CO"]]


# ---------------------------------------------------------------------------
# Membrane-reactor thermodynamic limit
# ---------------------------------------------------------------------------
def membrane_reactor_limit(F0, T, P_ret, P_perm):
    """Ideal-membrane, infinite-length limit for an isothermal MR with pure-H2 permeate.

    At the retentate outlet two conditions hold simultaneously:
      (1) WGS equilibrium:          K (CO)(H2O) = (CO2)(H2)
      (2) zero H2 driving force:    y_H2 * P_ret = P_perm
    Unknowns: extent xi and permeated H2 flow Fp.

    Returns dict with xi, Fp, retentate flows, X_CO and H2 recovery
    R = Fp / (H2_in + CO_in).
    If the conventional equilibrium already has y_H2 * P_ret <= P_perm, no H2
    can permeate and the conventional result is returned.
    """
    F0 = np.asarray(F0, dtype=float)
    a, d = F0[IDX["CO"]], F0[IDX["H2O"]]
    b, c = F0[IDX["CO2"]], F0[IDX["H2"]]
    N0 = F0.sum()
    r = P_perm / P_ret
    K = K_wgs(T)

    def h2_out(xi):                    # from condition (2) with N = N0 - Fp
        return r * (N0 - c - xi) / (1.0 - r)

    def resid(xi):
        return K * (a - xi) * (d - xi) - (b + xi) * h2_out(xi)

    xi_conv = wgs_extent_eq(F0, T)
    F_conv = F0 + NU_WGS * xi_conv
    if F_conv[IDX["H2"]] / F_conv.sum() * P_ret <= P_perm:
        Fp, xi, F_ret = 0.0, xi_conv, F_conv
    else:
        hi = min(a, d) * (1 - 1e-12)
        xi = brentq(resid, xi_conv, hi, xtol=1e-16, rtol=1e-14, maxiter=500)
        H2o = h2_out(xi)
        Fp = c + xi - H2o
        F_ret = F0 + NU_WGS * xi
        F_ret[IDX["H2"]] -= Fp
    return {
        "xi": xi,
        "Fp": Fp,
        "F_ret": F_ret,
        "X_CO": co_conversion(F0, F_ret),
        "R_H2": Fp / (c + a),
    }


# ---------------------------------------------------------------------------
# Adiabatic temperature
# ---------------------------------------------------------------------------
def adiabatic_T(F0, T_in, X_CO, T_bracket=(250.0, 2000.0)):
    """Outlet T [K] of an adiabatic, membrane-free reactor at CO conversion X_CO."""
    F0 = np.asarray(F0, dtype=float)
    F = F0 + NU_WGS * X_CO * F0[IDX["CO"]]
    H_in = mixture_enthalpy_flow(F0, T_in)
    return brentq(lambda T: mixture_enthalpy_flow(F, T) - H_in, *T_bracket, xtol=1e-9)


def adiabatic_equilibrium(F0, T_in):
    """Intersection of the adiabatic operating line with WGS equilibrium.

    Returns (X_CO, T_out). This is the most an adiabatic, membrane-free bed
    can achieve from inlet temperature T_in.
    """
    def g(X):
        T = adiabatic_T(F0, T_in, X)
        return X - co_conversion(F0, wgs_equilibrium(F0, T))
    X = brentq(g, 0.0, co_conversion(F0, wgs_equilibrium(F0, T_in)), xtol=1e-12)
    return X, adiabatic_T(F0, T_in, X)


# ---------------------------------------------------------------------------
# Full Gibbs equilibrium (side reactions allowed) and carbon activity
# ---------------------------------------------------------------------------
def gibbs_equilibrium(F0, T, P, species=("CO", "H2O", "CO2", "H2", "CH4", "N2"),
                      with_graphite=False):
    """Full Gibbs-energy minimisation at fixed T, P with Cantera.

    species: gas species allowed to form. Including CH4 allows methanation.
    with_graphite: also allow solid carbon C(gr) (nasa_condensed.yaml).
    Returns (gas mole flows [mol/s] in SPECIES order, graphite flow [mol/s]).
    """
    gas = _gas_subset(list(species))
    x = {s: F0[IDX[s]] for s in species if F0[IDX[s]] > 0}
    N0 = float(sum(F0[IDX[s]] for s in species))
    gas.TPX = T, P, x
    if not with_graphite:
        gas.equilibrate("TP")
        out = np.zeros(len(SPECIES))
        for s in species:
            out[IDX[s]] = gas[s].X[0] * N0
        return out, 0.0

    cs = [s for s in ct.Species.list_from_file("nasa_condensed.yaml") if s.name == "C(gr)"]
    graphite = ct.Solution(thermo="fixed-stoichiometry", species=cs)
    graphite.TP = T, P
    mix = ct.Mixture([(gas, N0), (graphite, 0.0)])
    mix.T, mix.P = T, P
    mix.equilibrate("TP", solver="vcs", max_steps=2000)
    out = np.zeros(len(SPECIES))
    n_gas = mix.phase_moles(0)
    for s in species:
        out[IDX[s]] = gas[s].X[0] * n_gas
    return out, float(mix.phase_moles(1))


def _g_graphite(T):
    cs = [s for s in ct.Species.list_from_file("nasa_condensed.yaml") if s.name == "C(gr)"][0]
    return (cs.thermo.h(T) - T * cs.thermo.s(T)) / 1000.0   # J/mol


def carbon_activity(F, T, P):
    """Carbon activities from Boudouard and CO reduction for gas flows F.

    a_C > 1 means graphite formation is thermodynamically favoured.
    Boudouard:     2 CO  <=> C + CO2       a_C = K_B  p_CO^2 / p_CO2
    CO reduction:  CO + H2 <=> C + H2O     a_C = K_R  p_CO p_H2 / p_H2O
    Partial pressures relative to P_REF.
    """
    F = np.asarray(F, dtype=float)
    p = F / F.sum() * P / P_REF
    _, _, g0 = species_props(T)
    gC = _g_graphite(T)
    RT = R * T
    K_B = np.exp(-(gC + g0[IDX["CO2"]] - 2 * g0[IDX["CO"]]) / RT)
    K_R = np.exp(-(gC + g0[IDX["H2O"]] - g0[IDX["CO"]] - g0[IDX["H2"]]) / RT)
    with np.errstate(divide="ignore"):
        aB = K_B * p[IDX["CO"]] ** 2 / p[IDX["CO2"]]
        aR = K_R * p[IDX["CO"]] * p[IDX["H2"]] / p[IDX["H2O"]]
    return {"boudouard": float(aB), "co_reduction": float(aR)}
