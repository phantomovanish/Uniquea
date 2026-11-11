"""Pi precision test (Uniquea, November 11, 2026).
Option A: recompute selected Uniquea quantities with pi to 15 decimals.
Option B: expand pi one decimal at a time (0 to 20), truncated and rounded, and watch the trend.
Reference ("truest"): pi to 50 decimals, Python Decimal at 60 significant digits.
Each sample also carries its data uncertainty: once pi's error is far below it, more decimals add no clarity.
Usage: python3 pi_precision.py   (writes results.csv, summary.json, pi_precision.png)
"""
from decimal import Decimal as D, getcontext, ROUND_DOWN, ROUND_HALF_EVEN
import csv, json
getcontext().prec = 60
PI = D("3.14159265358979323846264338327950288419716939937510")
c  = D(299792458)                      # m/s, exact
G  = D("6.67430e-11")                  # m^3 kg^-1 s^-2, CODATA 2018, rel. unc. 2.2e-5
MPC = D("3.0856775814913673e22")       # m per Mpc
H0 = D("67.4") * 1000 / MPC            # Planck 2018: 67.4 +/- 0.5 km/s/Mpc
R_EARTH = D("6371.0")                  # km, mean radius

samples = {
  "cH0/2pi (the a0 lead, m/s^2)":       (lambda p: c*H0/(2*p),            D("0.5")/D("67.4")),
  "Upper tap c^4/(80 pi G) (N)":         (lambda p: c**4/(80*p*G),          D("2.2e-5")),
  "e^(-pi) (the withdrawn 0.0432)":      (lambda p: (-p).exp(),             D(0)),
  "Earth circumference 2 pi R (km)":     (lambda p: 2*p*R_EARTH,            D("0.003")),
  "KFA = 8 pi^2 G / c^8":                (lambda p: 8*p*p*G/c**8,           D("2.2e-5")),
}

def pi_at(d, mode):
    q = D(1).scaleb(-d)
    return PI.quantize(q, rounding=ROUND_DOWN if mode == "truncated" else ROUND_HALF_EVEN)

rows, summary = [], {}
for name, (f, unc) in samples.items():
    truest = f(PI)
    summary[name] = {"truest": f"{truest:.15E}", "data_rel_uncertainty": float(unc),
                     "option_A_15_decimals": f"{f(pi_at(15,'truncated')):.15E}",
                     "option_A_rel_error": float(abs(f(pi_at(15,'truncated'))-truest)/abs(truest))}
    enough = None
    for d in range(0, 21):
        for mode in ("truncated", "rounded"):
            v = f(pi_at(d, mode)); rel = (v - truest) / truest
            rows.append([name, d, mode, str(pi_at(d, mode)), f"{v:.20E}", f"{float(rel):.3e}"])
        rel_t = abs(f(pi_at(d, "truncated")) - truest) / abs(truest)
        if enough is None and unc > 0 and rel_t < unc / 10:
            enough = d
    summary[name]["decimals_where_pi_error_drops_below_one_tenth_of_data_uncertainty"] = enough

with open("results.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample", "decimals", "mode", "pi_used", "value", "relative_error"]); w.writerows(rows)
json.dump(summary, open("summary.json", "w"), indent=1)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5.2), dpi=150)
name0 = "cH0/2pi (the a0 lead, m/s^2)"
for mode, style in (("truncated", "-o"), ("rounded", "--s")):
    xs, ys, signs = [], [], []
    for r in rows:
        if r[0] == name0 and r[2] == mode:
            e = float(r[5]); xs.append(r[1]); ys.append(max(abs(e), 1e-22)); signs.append(e)
    ax.semilogy(xs, ys, style, ms=4, label=f"pi {mode}")
ax.axhline(float(D('0.5')/D('67.4')), color="#b0452a", ls=":", label="H0 measurement uncertainty (0.74%)")
ax.axhline(2.2e-5, color="#6a6a6a", ls=":", label="G measurement uncertainty (0.0022%)")
ax.axhline(1.1e-16, color="#2a7a3a", ls="-.", lw=1, label="ordinary computer arithmetic (~16 digits)")
ax.set_xlabel("decimal places of pi"); ax.set_ylabel("relative error in the result")
ax.set_title("Pi precision: each decimal adds a ring of clarity, until the data's own error wins", fontsize=10)
ax.set_xticks(range(0, 21)); ax.grid(alpha=0.25); ax.legend(fontsize=7.5, frameon=False)
fig.tight_layout(); fig.savefig("pi_precision.png")
print(json.dumps(summary, indent=1))
for r in rows:
    if r[0] == name0 and r[1] in (0,1,2,3,4,5,8,15,20): print(r[1], r[2], r[3], r[5])
