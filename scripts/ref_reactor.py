"""Independent Python reference for the Step 3 reactor model (Cantera thermo, SciPy solver).
Reads the same register as MATLAB. Usage: python ref_reactor.py [path/to/params_register.csv]"""
import csv, math, sys
from pathlib import Path
import numpy as np
import cantera as ct
from scipy.integrate import solve_ivp

default = Path(__file__).resolve().parents[1] / 'matlab' / 'data' / 'params_register.csv'
REG = sys.argv[1] if len(sys.argv) > 1 else default
reg = {r['id']: r['value'] for r in csv.DictReader(open(REG, newline='', encoding='utf-8-sig'))}
g = lambda k: float(reg[k])

R = 8.31446261815324
gas = ct.Solution('gri30.yaml')
names = ['CO', 'H2O', 'CO2', 'H2', 'CH4', 'N2']
idx = [gas.species_index(n) for n in names]
NU = np.array([-1, -1, 1, 1, 0, 0.])


def thermo(T):
    gas.TP = T, 101325
    cp = gas.standard_cp_R[idx] * R
    h = gas.standard_enthalpies_RT[idx] * R * T
    return cp, h, gas.standard_gibbs_RT[idx]


def Keq(T):
    return math.exp(-float(NU @ thermo(T)[2]))


def rate_K1(T, P, y):
    y = np.maximum(y, 1e-6)
    beta = y[2] * y[3] / (y[0] * y[1] * Keq(T))
    c = y * P / (R * T) * 1e-6                      # mol cm-3, inferred basis (register)
    k = math.exp(g('rate.K1.ln_k0') - g('rate.K1.E') / (R * T))
    return 1000 * k * c[0] ** g('rate.K1.n_CO') * c[1] ** g('rate.K1.m_H2O') * (1 - beta)


def flux(T, pH2, pPerm):
    Q0 = g('membrane.M7.Pm0') / g('membrane.M7.delta')
    n = g('membrane.M7.n')
    J = Q0 * math.exp(-g('membrane.M7.Ea') / (R * T)) * (max(pH2, 0) ** n - max(pPerm, 0) ** n)
    return max(J, 0.0)


def base():
    return dict(
        F_in=np.array([g('cases.target.test_feed_mmol_s.' + n) for n in names]) * 1e-3,
        T_in=g('cases.base.T_in'), P_in=g('cases.base.P_ret'), P_perm=g('cases.target.P_perm'),
        L=g('cases.base.membrane_length'), D_m=g('cases.base.membrane_OD'), D_s=g('cases.base.shell_ID'),
        eps=g('cases.base.bed_voidage'), w_cat=g('cases.base.catalyst_wt_frac'))


def solve(F_in, T_in, P_in, P_perm, L, D_m, D_s, eps, w_cat, wall='isothermal', perm_scale=1.0, sweep=0.0, flow='co'):
    """flow 'co' or 'counter'. Without a sweep the direction does not matter. G = H2 in the permeate stream at z."""
    from scipy.optimize import brentq
    rho_cat = (1 - eps) * g('catalyst.Hla2015.rho_p') * w_cat
    A = math.pi / 4 * (D_s ** 2 - D_m ** 2)
    sgn = -1.0 if (flow == 'counter' and sweep > 0) else 1.0

    def rhs(z, y):
        F, G, T, P = y[:6], max(y[6], 0.0), y[7], y[8]
        yi = F / F.sum()
        rxn = rho_cat * rate_K1(T, P, yi) * A
        pperm = P_perm if sweep == 0 else P_perm * G / (G + sweep)
        perm = perm_scale * flux(T, yi[3] * P, pperm) * math.pi * D_m
        dF = NU * rxn
        dF[3] -= perm
        dT = 0.0
        if wall == 'adiabatic':
            cp, h, _ = thermo(T)
            dT = (-float(NU @ h) * rxn) / (F @ cp)
        return np.concatenate([dF, [sgn * perm, dT, 0.0]])

    def integrate(G0):
        y0 = np.concatenate([F_in, [G0, T_in, P_in]])
        return solve_ivp(rhs, (0, L), y0, method='BDF', t_eval=np.linspace(0, L, 201), rtol=1e-10, atol=1e-16)

    if sgn > 0:
        return integrate(0.0)
    hi = F_in[3] + F_in[0]
    G0 = brentq(lambda g0: integrate(g0).y[6, -1], 1e-12 * hi, hi, xtol=1e-16, rtol=1e-12)
    return integrate(G0)


if __name__ == '__main__':
    b = base()
    for label, kw in [('isothermal', {}), ('adiabatic', dict(wall='adiabatic')),
                      ('no membrane, isothermal', dict(perm_scale=0.0)),
                      ('no membrane, adiabatic', dict(perm_scale=0.0, wall='adiabatic'))]:
        y = solve(**{**b, **kw}).y[:, -1]
        print(f'{label:26s} X_CO={1 - y[0] / b["F_in"][0]:.8f}  Fp={y[6]:.8e} mol/s  T_out={y[7]:.6f} K')
    y = solve(**{**b, 'L': 0.01, 'perm_scale': 0.0}).y[:, -1]
    print(f'{"short bed 0.01 m, no mem":26s} X_CO={1 - y[0] / b["F_in"][0]:.6f}')
