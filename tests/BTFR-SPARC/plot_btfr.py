"""Plot for Test 2. Usage: python3 plot_btfr.py <SPARC_Lelli2016c.mrt> <results.json> <out.png>"""
import sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.argv, (tab, rj, out) = sys.argv[:1] + [sys.argv[1], "/dev/null"], sys.argv[1:4]
import importlib.util
src = open(__file__.replace("plot_btfr.py", "btfr_sparc.py")).read().split("res = dict(")[0]
ns = {"__name__": "x"}; exec(compile(src, "btfr", "exec"), ns)
x, y, sx, sy = ns["data"](0.5); r = json.load(open(rj))["orthogonal Y=0.5"]
G, MSUN = 6.674e-11, 1.989e30
fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=150)
ax.errorbar(x, y, xerr=sx, yerr=sy, fmt="o", ms=3.5, color="#2a6fb0", ecolor="#9bbbd9", elinewidth=0.8, label="SPARC, 123 galaxies (official Table1, Υ = 0.5)")
xx = np.linspace(1.3, 2.6, 50)
ax.plot(xx, r["slope"] * xx + r["intercept"], color="#1F3864", lw=1.6, label=f"Fit: slope {r['slope']:.2f} ± {r['slope_err']:.2f}")
a_fit = r["slope4_a0_weighted"]
for a0, ls, tag in ((1.2e-10, "--", "standard value"), (a_fit, ":", "fitted here")):
    I = -np.log10(G * a0 * MSUN / 1e12)
    ax.plot(xx, 4 * xx + I, ls, color="#b0452a", lw=1.4, label=f"Rung zero, slope 4, a₀ = {a0/1e-10:.2f}×10⁻¹⁰ m/s² ({tag})")
ax.set_xlabel("log V_flat (km/s)"); ax.set_ylabel("log M_b (solar masses)")
ax.set_title("Test 2: baryonic Tully–Fisher relation, M_b = V⁴/(G a₀)", fontsize=11)
ax.legend(fontsize=8, loc="upper left", frameon=False); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(out)
