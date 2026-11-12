"""check_rate_units.py - which unit basis makes the literature rate laws physically consistent?

Method: integrate an isothermal plug-flow bed at the conditions the papers actually ran, under each unit
hypothesis, and compare X_CO with equilibrium X_CO.  A hypothesis that gives X = Xeq (or ~0) at every
condition cannot be the law that produced the paper's measurements.
Assumptions NOT from the papers: Keiski bulk density 1.2 g/cm3 (not printed), space velocity S referred to
dry gas at 273.15 K and 1 atm per dm3 of bed, isothermal operation, Keq from Moe (1962).
This is a consistency check, not a reproduction of measured conversions.
Run:  python check_rate_units.py        (needs scipy)
"""
import math, itertools, statistics
from scipy.integrate import solve_ivp
R = 8.314462618

def Keq(T): return math.exp(4577.8 / T - 4.33)                      # Moe 1962 (within ~1 % of NASA at 300-470 C)
def beta(T, y):
    yCO, yH2O, yCO2, yH2 = y
    return (yCO2 * yH2) / (yCO * yH2O * Keq(T)) if yCO * yH2O > 0 else 1.0

def rate_K1(T, P, y, scale, per):    # c in (1/scale) mol per unit volume; rate divided by `per`
    c = P / (R * T)
    return (math.exp(26.1) * math.exp(-95000 / (R * T)) * (c * y[0] * scale) ** 1.1
            * (c * y[1] * scale) ** 0.53 * (1 - beta(T, y)) / per)

def rate_G1(T, P, y, f):             # partial pressures in Pa*f
    return 2.2e11 * math.exp(-11426.5 / T) * (P * y[0] * f) ** 1.1 * (P * y[1] * f) ** 0.53 * (1 - beta(T, y))

def rate_G2(T, P, y, f):
    p = [P * v * f for v in y]
    return (10 ** 0.659 * math.exp(-88000 / (R * T)) * p[0] ** 0.9 * p[1] ** 0.31
            * max(p[2], 1e-30) ** -0.156 * max(p[3], 1e-30) ** -0.05 * (1 - beta(T, y)))

def pfr(rate, T, P, feed, W):        # feed: mol/s of CO, H2O, CO2, H2, inert ; W in g
    n0 = [feed['CO'], feed['H2O'], feed['CO2'], feed['H2'], feed['inert']]
    def f(_, n):
        nt = sum(n); y = [max(n[i], 0) / nt for i in range(4)]
        r = rate(T, P, y)
        return [-r, -r, r, r, 0]
    s = solve_ivp(f, [0, W], n0, method='LSODA', rtol=1e-8, atol=1e-14)
    return 1 - s.y[0, -1] / n0[0]

def eqX(T, feed):
    lo, hi = 0.0, feed['CO']
    for _ in range(100):
        x = (lo + hi) / 2
        K = (feed['CO2'] + x) * (feed['H2'] + x) / ((feed['CO'] - x) * (feed['H2O'] - x))
        lo, hi = (lo, x) if K > Keq(T) else (x, hi)
    return x / feed['CO']

def keiski_corners(rate_of_hyp, rho=1.2):
    """32 corners of Keiski Tables 2 (set 1): CO 5/10 %, H2O/dry 0.83/1.24, T 573/623 K, S 2000/7000 1/h, H2 55/80 %, p 3/9 bar."""
    out = []
    for yco, w, T, S, h2, p in itertools.product((5, 10), (0.83, 1.24), (573, 623), (2000, 7000), (55, 80), (3, 9)):
        nd = S / 22.414 / 3600
        feed = {'CO': nd * yco / 100, 'H2': nd * h2 / 100, 'CO2': 0.0, 'H2O': nd * w, 'inert': nd * (1 - yco / 100 - h2 / 100)}
        out.append(pfr(lambda T_, P_, y_: rate_of_hyp(T_, P_, y_), T, p * 1e5, feed, rho * 1000) / eqX(T, feed))
    return min(out), statistics.median(out), max(out)

def main():
    print("== K1 (Keiski) : X_CO / X_eq over the 32 corners of Keiski's own design (min / median / max) ==")
    for name, scale, per in [("c in mol/dm3, rate mol/g/s   (as printed)", 1e-3, 1),
                             ("c in mol/dm3, rate mol/g/h", 1e-3, 3600),
                             ("c in mol/m3 , rate mol/g/s", 1.0, 1),
                             ("c in mol/cm3, rate mol/g/s", 1e-6, 1),
                             ("c in mol/cm3, rate mol/g/h", 1e-6, 3600)]:
        print("  %-42s %.3g / %.3g / %.3g" % ((name,) + keiski_corners(lambda T, P, y: rate_K1(T, P, y, scale, per))))
    print("\n== G1 (Barbieri 2011 form, 2.2e11 exp(-11426.5/T) P^1.1 P^0.53) on the same corners ==")
    for lab, f in (("P in Pa", 1.0), ("P in kPa", 1e-3), ("P in bar", 1e-5)):
        print("  %-42s %.3g / %.3g / %.3g" % ((lab + ", mol/g/s",) + keiski_corners(lambda T, P, y: rate_G1(T, P, y, f))))
        print("  %-42s %.3g / %.3g / %.3g" % ((lab + ", mol/m3/s per bed vol",) + keiski_corners(lambda T, P, y: rate_G1(T, P, y, f) / 1.2e6)))
    print("\n== G2 (Hla) vs Hla 2015 WGS-only runs: 350 C inlet, 16.01 bar abs, dry CO/H2/CO2 = 64.5/33/2.5, S:C 3 (carbon = CO+CO2) ==")
    print("   Paper: CO conversion < 10 % for 200/450/780 mm beds at 4-10 L/min; Run 5 (W/F = 11.6 kg/(mol CO/s)) approaches equilibrium")
    for Ldry, W in ((4, 3.22), (10, 3.22), (4, 12.18), (10, 12.18)):
        nd = Ldry / 22.414 / 60; nco, nco2 = nd * .645, nd * .025
        feed = {'CO': nco, 'H2O': 3 * (nco + nco2), 'CO2': nco2, 'H2': nd * .33, 'inert': 0}
        row = "  ".join("%s %.3g" % (u, pfr(lambda T, P, y: rate_G2(T, P, y, f), 623.15, 16.01e5, feed, W)) for u, f in
                        (("Pa", 1.0), ("kPa", 1e-3), ("bar", 1e-5)))
        print("  %2d L/min, %5.2f g (W/F %.2f)  Xeq %.2f | %s" % (Ldry, W, W / 1000 / nco, eqX(623.15, feed), row))
    nd = 4 / 22.414 / 60; nco, nco2 = nd * .645, nd * .025
    feed = {'CO': nco, 'H2O': 3 * (nco + nco2), 'CO2': nco2, 'H2': nd * .33, 'inert': 0}
    print("  Run-5-like (W/F 11.6): " + "  ".join("%s %.3g" % (u, pfr(lambda T, P, y: rate_G2(T, P, y, f), 623.15, 16.01e5, feed, 11.6 * nco * 1000))
          for u, f in (("Pa", 1.0), ("kPa", 1e-3), ("bar", 1e-5))) + "   Xeq %.3f" % eqX(623.15, feed))

if __name__ == "__main__":
    main()
