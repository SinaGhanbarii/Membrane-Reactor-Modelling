"""Phase 1: thermodynamic bounds for the 555.003 WGS membrane reactor.

Usage:  python scripts/phase1_equilibrium.py [params/rig_private.yaml]
Writes figures and summary.md to results/phase1/.
"""
import sys
from pathlib import Path

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pbmr1d import thermo as t  # noqa: E402

OUT = ROOT / "results" / "phase1"
OUT.mkdir(parents=True, exist_ok=True)

# ---- style (reference palette slots 1-3, recessive axes) -------------------
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "legend.frameon": False,
    "savefig.dpi": 200, "figure.facecolor": "white",
})

# ---- case ----------------------------------------------------------------
cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "params" / "rig_private.yaml"
if not cfg_path.exists():
    cfg_path = ROOT / "params" / "rig_example.yaml"
cfg = yaml.safe_load(open(cfg_path))
fd, op = cfg["feed"], cfg["operating"]
F0, finfo = t.feed_flows(fd["syngas_mass_flow_kg_h"], fd["composition_mol_pct"],
                         fd["water_mass_flow_kg_h"])
i = t.IDX
SCO_nom = F0[i["H2O"]] / F0[i["CO"]]
P_ret_list = [p * 1e5 + t.P_REF for p in op["P_ret_barg"]]
P_perm_list = [p * 1e5 for p in op["P_perm_bara"]]
P_ret_lo = P_ret_list[0]

lines = [f"# Phase 1 results: thermodynamic bounds\n",
         f"Case file: `{cfg_path.relative_to(ROOT)}` ({cfg['source']})\n"]

def kgh(x):  return x * 3600.0
def molh(x): return x * 3600.0

# ---- 0. Feed and Table 1 mass-balance check -------------------------------
lines.append("## 1. Feed and document consistency\n")
lines.append("| Species | Feed [mol/h] | mol % (wet) |\n|---|---|---|")
for s in t.SPECIES:
    lines.append(f"| {s} | {molh(F0[i[s]]):.2f} | {100*F0[i[s]]/F0.sum():.1f} |")
lines.append(f"\nSyngas mean molar mass {finfo['mw_syngas_g_mol']:.2f} g/mol; "
             f"total wet feed {molh(F0.sum()):.1f} mol/h; **S/CO = {SCO_nom:.2f}**; "
             f"steam/dry gas = {F0[i['H2O']]/(F0.sum()-F0[i['H2O']]):.2f}.\n")
H2_max = F0[i["H2"]] + F0[i["CO"]]
perm_nom = cfg["nominal_permeate"]["H2_kg_h"] / 3600 / t.MW[i["H2"]]
R_nom = perm_nom / H2_max
lines.append(f"- Nominal permeate {cfg['nominal_permeate']['H2_kg_h']} kg/h = "
             f"{molh(perm_nom):.2f} mol/h H2. Feed H2 = {molh(F0[i['H2']]):.2f} mol/h; "
             f"maximum after full shift = {molh(H2_max):.2f} mol/h. "
             f"**Implied H2 recovery = {100*R_nom:.1f} %.**")
for Tref, lab in ((273.15, "0 C"), (288.15, "15 C"), (293.15, "20 C")):
    slpm = perm_nom * 60 * t.R * Tref / t.P_REF * 1000
    lines.append(f"- 0.02 kg/h expressed at 1 atm, {lab}: {slpm:.2f} SLPM (document: "
                 f"{cfg['nominal_permeate']['H2_SLPM']} SLPM)")
if "original_reactor" in cfg:
    o = cfg["original_reactor"]
    Fo, oinfo = t.feed_flows(o["syngas_mass_flow_kg_h"], o["composition_mol_pct"],
                             o["water_mass_flow_kg_h"])
    lines.append(f"- 'Original reactor' composition sums to {oinfo['sum_pct']:.2f} % "
                 f"(normalised here). Its S/CO = {Fo[i['H2O']]/Fo[i['CO']]:.2f}; "
                 f"max H2 = {molh(Fo[i['H2']]+Fo[i['CO']]):.2f} mol/h, so the same "
                 f"0.02 kg/h implies {100*perm_nom/(Fo[i['H2']]+Fo[i['CO']]):.1f} % recovery.")
lines.append("")

# ---- 1. Conventional equilibrium map --------------------------------------
TC = np.linspace(250, 500, 126)
T = TC + 273.15
SCO = np.linspace(1.0, 4.5, 71)
X = np.zeros((len(SCO), len(T)))
for a, sco in enumerate(SCO):
    F = F0.copy(); F[i["H2O"]] = sco * F0[i["CO"]]
    for b, Tk in enumerate(T):
        X[a, b] = t.co_conversion(F, t.wgs_equilibrium(F, Tk))

fig, ax = plt.subplots(figsize=(6.4, 4.4))
lev = [0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99]
cs = ax.contour(TC, SCO, X, levels=lev, colors=C1, linewidths=1.6)
ax.clabel(cs, fmt=lambda v: f"{100*v:.0f} %", fontsize=8, colors=INK)
ax.axhline(SCO_nom, color=INK2, lw=1, ls="--")
ax.axvspan(op["T_inlet_C"], op["T_outlet_max_C"], color=C2, alpha=0.08, lw=0)
ax.text(405, SCO_nom + 0.08, f"nominal S/CO = {SCO_nom:.2f}", color=INK2, fontsize=8)
ax.text(op["T_inlet_C"] + 3, 1.08, "documented reactor T window", color=C2, fontsize=8)
ax.set_xlabel("Temperature [°C]"); ax.set_ylabel("Steam / CO [mol/mol]")
ax.set_title("Equilibrium CO conversion, WGS only (independent of pressure)", fontsize=10, loc="left")
fig.tight_layout(); fig.savefig(OUT / "fig1_Xeq_map.png"); plt.close(fig)

lines.append("## 2. Conventional WGS equilibrium (no membrane)\n")
lines.append("| T [°C] | K_eq (Cantera) | K_eq (Moe form) | ΔH_WGS [kJ/mol] | X_CO,eq at nominal S/CO |\n|---|---|---|---|---|")
for c in (250, 300, 350, 400, 450, 470, 500):
    Tk = c + 273.15
    lines.append(f"| {c} | {t.K_wgs(Tk):.2f} | {t.K_wgs_moe(Tk):.2f} | {t.dH_wgs(Tk)/1e3:.1f} | "
                 f"{100*t.co_conversion(F0, t.wgs_equilibrium(F0, Tk)):.1f} % |")
lines.append("")

# ---- 2. Membrane-reactor limit vs permeate pressure -----------------------
Pp = np.linspace(0.05e5, 1.5e5, 59)
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
cols = {350: C3, 400: C1, 450: C2}
for c, col in cols.items():
    Tk = c + 273.15
    Xc = t.co_conversion(F0, t.wgs_equilibrium(F0, Tk))
    mr = [t.membrane_reactor_limit(F0, Tk, P_ret_lo, p) for p in Pp]
    axs[0].plot(Pp / 1e5, [100 * m["X_CO"] for m in mr], color=col, label=f"{c} °C, MR limit")
    axs[0].axhline(100 * Xc, color=col, lw=1, ls=":")
    axs[1].plot(Pp / 1e5, [100 * m["R_H2"] for m in mr], color=col, label=f"{c} °C")
axs[0].plot([], [], color=INK2, lw=1, ls=":", label="conventional equilibrium")
for ax in axs:
    for p in P_perm_list:
        ax.axvline(p / 1e5, color=INK2, lw=0.8, ls="--")
axs[1].axhline(100 * R_nom, color=INK, lw=1.2, ls="-.")
axs[1].text(0.45, 100 * R_nom - 6, f"H2Site nominal permeate ({100*R_nom:.0f} %)", fontsize=8, color=INK)
axs[0].set_ylabel("CO conversion [%]"); axs[1].set_ylabel("Max. H2 recovery [%]")
axs[0].set_ylim(85, 100.5); axs[1].set_ylim(0, 100)
for ax in axs: ax.set_xlabel("Permeate pressure [bara]")
axs[0].set_title(f"(a) CO conversion, P_ret = {P_ret_lo/1e5:.2f} bara", fontsize=10, loc="left")
axs[1].set_title("(b) H2 recovery = F_perm / (H2,in + CO,in)", fontsize=10, loc="left")
axs[0].legend(fontsize=8, loc="lower left"); axs[1].legend(fontsize=8, loc="lower left")
fig.tight_layout(); fig.savefig(OUT / "fig2_membrane_limit.png"); plt.close(fig)

lines.append("## 3. Membrane-reactor thermodynamic limit (ideal membrane, infinite length, isothermal)\n")
lines.append("| T [°C] | P_ret [bara] | P_perm [bara] | X_CO conv. eq. | X_CO MR limit | R_H2 max | Max permeate [mol/h] | Nominal / max |\n|---|---|---|---|---|---|---|---|")
for c in (350, 400, 450):
    Tk = c + 273.15
    Xc = t.co_conversion(F0, t.wgs_equilibrium(F0, Tk))
    for Pr in P_ret_list:
        for pp in P_perm_list:
            m = t.membrane_reactor_limit(F0, Tk, Pr, pp)
            lines.append(f"| {c} | {Pr/1e5:.2f} | {pp/1e5:.2f} | {100*Xc:.1f} % | {100*m['X_CO']:.2f} % | "
                         f"{100*m['R_H2']:.1f} % | {molh(m['Fp']):.1f} | {100*perm_nom/m['Fp']:.0f} % |")
lines.append("")

# ---- 3. Steam penalty on recovery ----------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
SCOs = np.linspace(1.0, 5.0, 81)
Tk = 400 + 273.15
for pp, col in zip([0.1e5, 0.2e5, 0.5e5, 1.0e5], [C3, C1, C2, INK2]):
    Rs, Xs = [], []
    for sco in SCOs:
        F = F0.copy(); F[i["H2O"]] = sco * F0[i["CO"]]
        m = t.membrane_reactor_limit(F, Tk, P_ret_lo, pp)
        Rs.append(100 * m["R_H2"]); Xs.append(100 * m["X_CO"])
    axs[0].plot(SCOs, Xs, color=col, label=f"P_perm = {pp/1e5:.2f} bara")
    axs[1].plot(SCOs, Rs, color=col, label=f"P_perm = {pp/1e5:.2f} bara")
Xconv = []
for sco in SCOs:
    F = F0.copy(); F[i["H2O"]] = sco * F0[i["CO"]]
    Xconv.append(100 * t.co_conversion(F, t.wgs_equilibrium(F, Tk)))
axs[0].plot(SCOs, Xconv, color=INK, lw=1, ls=":", label="conventional equilibrium")
for ax in axs:
    ax.axvline(SCO_nom, color=INK2, lw=0.8, ls="--"); ax.set_xlabel("Steam / CO [mol/mol]")
axs[0].set_ylabel("CO conversion [%]"); axs[1].set_ylabel("Max. H2 recovery [%]")
axs[0].set_title(f"(a) CO conversion, 400 °C, P_ret = {P_ret_lo/1e5:.2f} bara", fontsize=10, loc="left")
axs[1].set_title("(b) Steam dilutes the retentate and lowers the recovery ceiling", fontsize=10, loc="left")
axs[0].legend(fontsize=8); axs[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "fig3_steam_tradeoff.png"); plt.close(fig)

lines.append("## 4. Steam trade-off at 400 °C (membrane-reactor limit)\n")
lines.append("| S/CO | X_CO conv. | X_CO MR (0.2 bara) | R_H2 max (0.2 bara) | R_H2 max (1.0 bara) |\n|---|---|---|---|---|")
for sco in (1.5, 2.0, 2.5, 3.0, round(SCO_nom, 2), 4.5):
    F = F0.copy(); F[i["H2O"]] = sco * F0[i["CO"]]
    m2 = t.membrane_reactor_limit(F, Tk, P_ret_lo, 0.2e5)
    m10 = t.membrane_reactor_limit(F, Tk, P_ret_lo, 1.0e5)
    lines.append(f"| {sco} | {100*t.co_conversion(F, t.wgs_equilibrium(F, Tk)):.1f} % | {100*m2['X_CO']:.2f} % | "
                 f"{100*m2['R_H2']:.1f} % | {100*m10['R_H2']:.1f} % |")
lines.append("")

# ---- 4. Adiabatic operating lines -----------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 4.4))
Teq = np.linspace(250, 650, 161)
ax.plot(Teq, [100 * t.co_conversion(F0, t.wgs_equilibrium(F0, x + 273.15)) for x in Teq],
        color=INK, lw=1.5, label="equilibrium (no membrane)")
lines.append("## 5. Heat release and adiabatic temperature rise (no membrane)\n")
lines.append("| T_in [°C] | ΔT_ad at X_CO = 50 % | ΔT_ad at 100 % | Adiabatic-equilibrium X_CO | T_out [°C] |\n|---|---|---|---|---|")
for Tin, col in ((op["T_inlet_C"], C1), (350, C2)):
    Xs = np.linspace(0, 1, 51)
    Ts = [t.adiabatic_T(F0, Tin + 273.15, x) - 273.15 for x in Xs]
    Xae, Tae = t.adiabatic_equilibrium(F0, Tin + 273.15)
    ax.plot(Ts, 100 * Xs, color=col, label=f"adiabatic line, T_in = {Tin} °C")
    ax.plot(Tae - 273.15, 100 * Xae, "o", ms=8, color=col, mec="white", mew=2)
    lines.append(f"| {Tin} | {t.adiabatic_T(F0, Tin+273.15, 0.5)-273.15-Tin:.0f} K | "
                 f"{t.adiabatic_T(F0, Tin+273.15, 1.0)-273.15-Tin:.0f} K | {100*Xae:.1f} % | {Tae-273.15:.0f} |")
ax.axvline(op["T_outlet_max_C"], color=C2, lw=1, ls="--")
ax.text(op["T_outlet_max_C"] + 4, 5, "documented max. 470 °C", color=C2, fontsize=8, rotation=90)
ax.set_xlabel("Temperature [°C]"); ax.set_ylabel("CO conversion [%]"); ax.set_ylim(0, 102)
ax.set_title("Adiabatic operating lines vs WGS equilibrium", fontsize=10, loc="left")
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout(); fig.savefig(OUT / "fig4_adiabatic.png"); plt.close(fig)
lines.append("")

# ---- 5. Side reactions: methanation and carbon ----------------------------
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
TCs = np.linspace(250, 550, 61)
for Pr, col in zip(P_ret_list, [C1, C2]):
    ch4 = []
    for c in TCs:
        g, _ = t.gibbs_equilibrium(F0, c + 273.15, Pr)
        ch4.append(g[i["CH4"]] / F0[i["CH4"]])
    axs[0].plot(TCs, ch4, color=col, label=f"P = {Pr/1e5:.2f} bara")
axs[0].axhline(1, color=INK2, lw=1, ls=":")
axs[0].text(255, 1.15, "inert CH4 (model assumption)", fontsize=8, color=INK2)
axs[0].set_ylabel("CH4 out / CH4 in at full Gibbs equilibrium")
axs[0].set_title("(a) Methanation potential", fontsize=10, loc="left")
aB_eq, aB_feed = [], []
for c in TCs:
    Tk = c + 273.15
    aB_eq.append(t.carbon_activity(t.wgs_equilibrium(F0, Tk), Tk, P_ret_lo)["boudouard"])
    aB_feed.append(t.carbon_activity(F0, Tk, P_ret_lo)["boudouard"])
axs[1].semilogy(TCs, aB_feed, color=C2, label="feed gas (before shift)")
axs[1].semilogy(TCs, aB_eq, color=C1, label="WGS-equilibrium gas")
axs[1].axhline(1, color=INK2, lw=1, ls=":")
axs[1].text(470, 1.4, "a_C = 1", fontsize=8, color=INK2)
axs[1].set_ylabel("Carbon activity a_C (Boudouard)")
axs[1].set_title("(b) Graphite formation potential", fontsize=10, loc="left")
for ax in axs: ax.set_xlabel("Temperature [°C]"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "fig5_side_reactions.png"); plt.close(fig)

lines.append("## 6. Side reactions\n")
lines.append(f"| T [°C] | CH4 out/in (full Gibbs, {P_ret_lo/1e5:.2f} bara) | X_CO full Gibbs | H2 at full Gibbs / H2 at WGS eq. | a_C feed | a_C WGS eq. | Graphite at full eq. |\n|---|---|---|---|---|---|---|")
for c in (250, 300, 350, 400, 450, 500):
    Tk = c + 273.15
    g, _ = t.gibbs_equilibrium(F0, Tk, P_ret_lo)
    gw = t.wgs_equilibrium(F0, Tk)
    _, gr = t.gibbs_equilibrium(F0, Tk, P_ret_lo, with_graphite=True)
    lines.append(f"| {c} | {g[i['CH4']]/F0[i['CH4']]:.2f} | {100*t.co_conversion(F0, g):.1f} % | "
                 f"{g[i['H2']]/gw[i['H2']]:.2f} | {t.carbon_activity(F0, Tk, P_ret_lo)['boudouard']:.3g} | "
                 f"{t.carbon_activity(gw, Tk, P_ret_lo)['boudouard']:.3g} | {'yes' if gr > 1e-12 else 'none'} |")
lines.append("")

(OUT / "summary.md").write_text("\n".join(lines))
print("\n".join(lines))
