"""Correctness tests for pbmr1d.thermo (Phase 1).

Run with:  python -m pytest -q
"""
import numpy as np
import pytest

from pbmr1d import thermo as t

P_RET = 6.01325e5   # 5 barg
P_PERM = 0.2e5      # 0.2 bara

# Generic test feed (mol/s); not the rig composition
F0 = np.array([0.40, 1.20, 0.15, 0.35, 0.05, 0.05]) * 1e-3


@pytest.mark.parametrize("TC", [250, 300, 350, 400, 450, 500])
def test_keq_matches_empirical_correlation(TC):
    """Cantera K_eq within 10 % of the Moe-type correlation (Karagoz et al., 2020)."""
    T = TC + 273.15
    assert t.K_wgs(T) == pytest.approx(t.K_wgs_moe(T), rel=0.10)


def test_keq_decreases_with_temperature():
    """Exothermic reaction: K falls as T rises (van 't Hoff)."""
    T = np.linspace(500, 800, 31)
    K = np.array([t.K_wgs(x) for x in T])
    assert np.all(np.diff(K) < 0)
    assert all(t.dH_wgs(x) < 0 for x in T)


def test_vant_hoff_consistency():
    """d ln K / dT = dH / (R T^2), checked numerically."""
    T, h = 650.0, 1e-2
    dlnK = (np.log(t.K_wgs(T + h)) - np.log(t.K_wgs(T - h))) / (2 * h)
    assert dlnK == pytest.approx(t.dH_wgs(T) / (t.R * T**2), rel=1e-5)


@pytest.mark.parametrize("TC", [250, 400, 500])
def test_equilibrium_closes_element_balances(TC):
    F = t.wgs_equilibrium(F0, TC + 273.15)
    np.testing.assert_allclose(t.element_flows(F), t.element_flows(F0), rtol=1e-12, atol=1e-18)


@pytest.mark.parametrize("TC", [250, 400, 500])
def test_equilibrium_satisfies_mass_action(TC):
    T = TC + 273.15
    F = t.wgs_equilibrium(F0, T)
    i = t.IDX
    Q = F[i["CO2"]] * F[i["H2"]] / (F[i["CO"]] * F[i["H2O"]])
    assert Q == pytest.approx(t.K_wgs(T), rel=1e-9)


def test_large_excess_steam_gives_full_conversion():
    F = F0.copy()
    F[t.IDX["H2O"]] = 1e4 * F[t.IDX["CO"]]
    X = t.co_conversion(F, t.wgs_equilibrium(F, 673.15))
    assert X > 0.9999


def test_cantera_gibbs_reproduces_extent_solution():
    """With only CO, H2O, CO2, H2, N2 allowed (no CH4), Gibbs minimisation has a
    single independent reaction, so it must equal the extent-based solution."""
    F = F0.copy()
    F[t.IDX["CH4"]] = 0.0
    for TC in (300, 400, 500):
        T = TC + 273.15
        g, _ = t.gibbs_equilibrium(F, T, P_RET, species=("CO", "H2O", "CO2", "H2", "N2"))
        np.testing.assert_allclose(g, t.wgs_equilibrium(F, T), rtol=1e-6, atol=1e-12)


def test_pressure_does_not_change_conventional_equilibrium():
    F = F0.copy()
    F[t.IDX["CH4"]] = 0.0
    sp = ("CO", "H2O", "CO2", "H2", "N2")
    g1, _ = t.gibbs_equilibrium(F, 673.15, 1e5, species=sp)
    g2, _ = t.gibbs_equilibrium(F, 673.15, 20e5, species=sp)
    np.testing.assert_allclose(g1, g2, rtol=1e-6)


def test_membrane_limit_conditions_and_balances():
    T = 673.15
    m = t.membrane_reactor_limit(F0, T, P_RET, P_PERM)
    Fr, i = m["F_ret"], t.IDX
    # H2 driving force exhausted at outlet
    assert Fr[i["H2"]] / Fr.sum() * P_RET == pytest.approx(P_PERM, rel=1e-8)
    # Retentate in WGS equilibrium
    Q = Fr[i["CO2"]] * Fr[i["H2"]] / (Fr[i["CO"]] * Fr[i["H2O"]])
    assert Q == pytest.approx(t.K_wgs(T), rel=1e-6)
    # Elements: retentate + permeate (pure H2) = feed
    perm = np.zeros(6); perm[i["H2"]] = m["Fp"]
    np.testing.assert_allclose(t.element_flows(Fr + perm), t.element_flows(F0), rtol=1e-10)


def test_membrane_limit_exceeds_conventional_and_falls_with_permeate_pressure():
    T = 673.15
    X_conv = t.co_conversion(F0, t.wgs_equilibrium(F0, T))
    X = [t.membrane_reactor_limit(F0, T, P_RET, p)["X_CO"] for p in (0.1e5, 0.5e5, 1.0e5)]
    R = [t.membrane_reactor_limit(F0, T, P_RET, p)["R_H2"] for p in (0.1e5, 0.5e5, 1.0e5)]
    assert all(x > X_conv for x in X)
    assert X[0] > X[1] > X[2]
    assert R[0] > R[1] > R[2]


def test_membrane_limit_reduces_to_conventional_when_no_driving_force():
    m = t.membrane_reactor_limit(F0, 673.15, P_RET, P_RET * 0.999)
    assert m["Fp"] == 0.0


def test_adiabatic_energy_balance():
    T_in = 523.15
    X = 0.8
    T_out = t.adiabatic_T(F0, T_in, X)
    F = F0 + t.NU_WGS * X * F0[t.IDX["CO"]]
    assert t.mixture_enthalpy_flow(F, T_out) == pytest.approx(
        t.mixture_enthalpy_flow(F0, T_in), rel=1e-10)
    assert T_out > T_in


def test_carbon_activity_routes_agree_at_wgs_equilibrium():
    """At WGS equilibrium the Boudouard and CO-reduction carbon activities are
    linked by K_WGS, so they must be identical."""
    F = t.wgs_equilibrium(F0, 673.15)
    a = t.carbon_activity(F, 673.15, P_RET)
    assert a["boudouard"] == pytest.approx(a["co_reduction"], rel=1e-8)
