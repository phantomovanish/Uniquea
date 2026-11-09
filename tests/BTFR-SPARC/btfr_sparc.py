"""Test 2: the baryonic Tully-Fisher relation (BTFR) against SPARC.

Rung zero (deep MOND) predicts, for the flat outer speed of a disk galaxy,
    M_b = V_flat^4 / (G a0),   i.e.  log M_b = 4 log V_flat - log(G a0)  (slope exactly 4).
Data: SPARC Table1 (Lelli, McGaugh & Schombert 2016, AJ 152, 157), official file SPARC_Lelli2016c.mrt.
Benchmark: Lelli et al. 2019, MNRAS 484, 3267 (V_flat: N = 123, s = 3.85 +/- 0.09, I = 1.99 +/- 0.18,
sigma_perp = 0.026 +/- 0.007 dex, observed vertical scatter 0.24 dex).
Method follows the benchmark: Q < 3, i >= 30 deg, V_flat > 0 (5% flatness, already applied in the table);
M_b = 1.33 M_HI + Y L[3.6], Y = 0.5 (varied 0.4-0.6); errors on M_b from flux, distance and 0.11 dex scatter in Y;
errors on V_flat from the table plus the inclination term V di / tan i; orthogonal maximum likelihood with
constant orthogonal intrinsic scatter; posterior mapped with a plain Metropolis MCMC (numpy only).
Both likelihoods are reported: the benchmark's orthogonal one and the vertical one; they bracket the slope.
Usage: python3 btfr_sparc.py <SPARC_Lelli2016c.mrt> <out.json>
"""
import sys, re, json, numpy as np
tab, out = sys.argv[1:3]
G, MSUN = 6.674e-11, 1.989e30
rows = []
for l in open(tab, errors="ignore"):
    q = l.split()
    if len(q) >= 18 and re.match(r"^\d+$", q[1]):
        rows.append(dict(name=q[0], D=float(q[2]), eD=float(q[3]), inc=float(q[5]), einc=float(q[6]),
                         L=float(q[7]), eL=float(q[8]), MHI=float(q[13]), Vf=float(q[15]), eVf=float(q[16]), Q=int(q[17])))
sel = [r for r in rows if r["Q"] < 3 and r["inc"] >= 30 and r["Vf"] > 0]

def data(Y):
    L = np.array([r["L"] for r in sel]) * 1e9; eL = np.array([r["eL"] for r in sel]) * 1e9
    H = np.array([r["MHI"] for r in sel]) * 1e9
    D = np.array([r["D"] for r in sel]); eD = np.array([r["eD"] for r in sel])
    V = np.array([r["Vf"] for r in sel]); eV = np.array([r["eVf"] for r in sel])
    i = np.radians([r["inc"] for r in sel]); ei = np.radians([r["einc"] for r in sel])
    Mb = 1.33 * H + Y * L
    sY = Y * L * np.log(10) * 0.11                          # 0.11 dex scatter in Y
    eMb = np.sqrt((Y * eL) ** 2 + (0.1 * 1.33 * H) ** 2 + sY ** 2 + (2 * eD / D * Mb) ** 2)  # 10% HI flux error
    eVt = np.sqrt(eV ** 2 + (V * ei / np.tan(i)) ** 2)
    return np.log10(V), np.log10(Mb), eVt / (V * np.log(10)), eMb / (Mb * np.log(10))

ORTH = True   # set by run(); True = benchmark's orthogonal likelihood, False = vertical likelihood
def lnlike(p, x, y, sx, sy):
    s, I, sp = p
    if sp <= 0 or sp > 1: return -np.inf
    d = y - s * x - I
    if ORTH:   # Gaussian in the orthogonal distance (as in Lelli et al. 2019)
        var = (sy ** 2 + (s * sx) ** 2) / (1 + s ** 2) + sp ** 2
        return -0.5 * np.sum(d ** 2 / (1 + s ** 2) / var + np.log(var))
    var = sy ** 2 + (s * sx) ** 2 + sp ** 2 * (1 + s ** 2)   # vertical residual, orthogonal intrinsic scatter projected
    return -0.5 * np.sum(d ** 2 / var + np.log(var))

def mcmc(x, y, sx, sy, n=60000, seed=1):
    rng = np.random.default_rng(seed); p = np.array([3.8, 2.0, 0.03]); lp = lnlike(p, x, y, sx, sy)
    step = np.array([0.04, 0.09, 0.004]); ch = []; best = (lp, p.copy())
    for k in range(n):
        q = p + step * rng.standard_normal(3); lq = lnlike(q, x, y, sx, sy)
        if np.log(rng.random()) < lq - lp: p, lp = q, lq
        if lp > best[0]: best = (lp, p.copy())
        if k >= n // 4: ch.append(p.copy())
    ch = np.array(ch); return best[1], ch.std(0), ch

def fixed4(x, y, sx, sy):
    # slope held at 4: a0 = 10^(4x - y) / G, in SI; weighted mean of log a0 with intrinsic scatter ignored
    la = 4 * x - y + np.log10((1e3) ** 4 / (G * MSUN))
    w = 1 / (sy ** 2 + 16 * sx ** 2)
    m = np.sum(w * la) / np.sum(w)
    return 10 ** m, 10 ** np.median(la), la.std(), np.sqrt(1 / np.sum(w))

res = dict(N=len(sel), benchmark=dict(N=123, slope=3.85, slope_err=0.09, intercept=1.99, intercept_err=0.18,
                                       sigma_perp=0.026, sigma_perp_err=0.007, obs_scatter=0.24))
for lik in ("orthogonal", "vertical"):
  ORTH = lik == "orthogonal"
  for Y in (0.4, 0.5, 0.6):
    x, y, sx, sy = data(Y)
    best, err, ch = mcmc(x, y, sx, sy)
    obs = float(np.sqrt(np.mean((y - best[0] * x - best[1]) ** 2)))
    a_w, a_med, a_sc, a_err = fixed4(x, y, sx, sy)
    res[f"{lik} Y={Y}"] = dict(slope=round(float(best[0]), 3), slope_err=round(float(err[0]), 3),
                         intercept=round(float(best[1]), 3), intercept_err=round(float(err[1]), 3),
                         sigma_perp=round(float(best[2]), 4), sigma_perp_err=round(float(err[2]), 4),
                         obs_vertical_scatter=round(obs, 3),
                         slope4_a0_weighted=float(f"{a_w:.3e}"), slope4_a0_median=float(f"{a_med:.3e}"),
                         slope4_log_a0_scatter_dex=round(float(a_sc), 3), slope4_a0_stat_err_dex=round(float(a_err), 4),
                         slope_ge_4_fraction=round(float(np.mean(ch[:, 0] >= 4)), 3))
res["galaxies"] = [r["name"] for r in sel]
json.dump(res, open(out, "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "galaxies"}, indent=1))
