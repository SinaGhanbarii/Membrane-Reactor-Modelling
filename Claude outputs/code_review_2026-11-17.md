# Code review of the 1D Pd-membrane WGS reactor model (MATLAB), 17 Nov 2026

Scope: every file under `matlab/src` (kinetics, membrane, params, reactor, thermo, studies, plots) plus the walkthrough. Read line by line.
Limits: MATLAB was not available to the reviewer. Findings come from reading the code, plus a few numerical checks in the Python reference. Items marked "reasoned" were not run.

## High

| # | Where | Problem | Evidence | Suggested fix |
|---|---|---|---|---|
| H1 | `load_params.m`, `load_species.m` | Both cache their result (`persistent cache`). After the register or the JSON is edited, MATLAB keeps the old values until `clear functions` or `clear all`. Today's register edits (G2 pressure unit now `reported`, G2 validity window) are not seen by an open session, and `run_tests_step2` line 15 will fail on the stale `inferred`. | Read from code | Run `clear functions` before tests; or reload when the file's modification time changes |
| H2 | `solve_countercurrent.m` | (a) With no net permeation possible (`perm_scale = 0` with a sweep and `flow = 'counter'`), G(L) stays at the lower bound and `miss(lo) > 0` raises `solve_countercurrent:bracket`. (b) `endpoint` uses `s.Fp(end)` even if `ode15s` stopped early, so `fzero` can converge to a root of a truncated run; only the final run is checked. | Reasoned from code | If `perm_scale == 0`, set G0 = 0 and integrate once. Make `endpoint` error when `numel(s.z) < 201` |
| H3 | `membrane_flux.m`, `check_validity.m`, `run_case.m` | Membrane validity is never checked. The register gives M7 a range (T 623-723 K, dP 1.2-1.6 bar) and M1 a range (`P_kPa_ret`), but only `rate_wgs` calls `check_validity`, and `check_validity` knows only `T_K`, `P_bar`, `yCO_dry`. The run's EXTRAPOLATED flags therefore cover the rate law only. The membrane flux is the largest extrapolation in the model (1.2-1.6 bar data used at 4-5 bar). | Read from code and register | Add membrane flags: T against `T_K`, and P_ret - P_perm against `dP_bar`. Pass them into `sm.flags` |
| H4 | `make_case.m` | Overrides are not sanity-checked. `D_s <= D_m` gives a negative bed area and a negative reaction rate with no error. A value in mm instead of m, a temperature in C instead of K, `P_in <= P_perm`, `eps` outside (0,1), `L <= 0`, a row-vector `F_in` or a NaN `perm_scale` all pass or fail obscurely. This matters most when H2Site values are entered by hand. | Read from code | Add range and unit checks with clear error ids |

## Medium

| # | Where | Problem | Suggested fix |
|---|---|---|---|
| M1 | `reactor_rhs.m` (dp = on) | The pressure can fall to zero or below on a long bed or small particles. A negative base to a fractional power in the rate laws then gives complex numbers without an error. | Error or clip when P <= 0 |
| M2 | `membrane_flux.m` | The inhibition factor (1 - alpha*theta) is unclipped. Currently positive (Melendez set: 0.27 at 82 kPa CO, 0.02 at 300 kPa), but any alpha > 1 set can go negative, and with `backperm = on` and a sweep a negative factor would reverse the flux. | `max(1 - alpha*theta, 0)` |
| M3 | `study_regimes.m` vs `run_case.m` | Two regime definitions. The study classifies by outlet beta and outlet pinch ratio. `run_case` labels every position along the bed and has the idle and reverse-flux classes. The regime map and the walkthrough bar chart can disagree. The study also ignores sweep, counter-current and `backperm`. | Use `run_case`'s regime fractions in the study |
| M4 | `plot_regime_map.m` | Assumes exactly one row per condition. A failed case (row with an error message and NaN fields) breaks parts C and D with a dimension error. | Skip rows with NaN or an error |
| M5 | `reactor_rhs.m` (energy) | The energy balance has no heat capacity for the sweep stream and no permeate-retentate heat exchange. For adiabatic or cooled counter-current runs with 1-5 mmol/s sweep this matters (similar order to the retentate heat capacity). The isothermal default hides it. | State as a limitation; add later |

## Low

- L1 `rate_wgs.m`: mole fractions are floored at 1e-6, which caps X_CO near 0.99999 and keeps CO from reaching equilibrium below that.
- L2 `load_params.m`: a numeric typo such as `0,4` becomes text, and arithmetic on text gives ASCII codes with no error. Check numeric-looking strings.
- L3 `nasa_props.m`: a trial temperature outside 200-3500 K inside `ode15s` raises an error, so the run aborts instead of rejecting the step.
- L4 `rate_wgs.m` / `reactor_rhs.m`: flags are rebuilt (with `sprintf`) at every right-hand-side call, which slows counter-current shooting (about 40 integrations).
- L5 `study_h2.m`: `rec_change` divides by the baseline recovery, giving Inf or NaN when it is zero.
- L6 `reactor_rhs.m`: without a sweep, the permeate is treated as pure H2 at P_perm, so back-permeation is clamped. This is an assumption, not a bug.
- L7 `param_usage.m`: a `pget` call inside a trailing comment is counted.

## Checked and found sound

- Units and mass balance in `reactor_rhs` (rate per kg catalyst times bed density times area; flux times membrane perimeter); the energy balance sign and the cancellation of the permeate enthalpy.
- `wgs_extent_eq` quadratic and root selection over the model's temperature range.
- Element matrix and `dp` density term dimensions (`ATOMS'` is 6x4).
- `study_*` file paths (results go to `results/` at the repository root, ignored by Git).
- The `ode15s` setup (tolerances, 201 output points, truncation check in the final run).
