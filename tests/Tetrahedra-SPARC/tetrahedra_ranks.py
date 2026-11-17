"""Test 3, follow-up (November 17, 2026): how unusual is the real sky, counted directly?

The first run (tetrahedra_sparc.py, 300 shuffles) reported z = (real - shuffle mean) / shuffle spread,
which assumes the shuffles form a bell curve. This script also COUNTS: how many shuffles came out as
extreme as the real sky (one-sided, in the predicted direction), and repeats with many more shuffles
so very small chances can be resolved.
Usage: python3 tetrahedra_ranks.py <SPARC_Lelli2016c.mrt> <sparc_sky_coordinates_ned.csv> <n_shuffles> <seed>
"""
import sys, json, numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.stats import norm, skew
tab, coordfile, N, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
co = {}
for l in open(coordfile):
    p = l.strip().split(',')
    if len(p) == 4 and p[1] and p[0] != 'name': co[p[0]] = (float(p[1]), float(p[2]))
D, fD = {}, {}
for l in open(tab, errors="ignore"):
    f = l.split()
    if len(f) >= 18 and f[1].isdigit(): D[f[0]] = float(f[2]); fD[f[0]] = int(f[4])
names = [n for n in co if n in D]
ra = np.radians([co[n][0] for n in names]); de = np.radians([co[n][1] for n in names])
d = np.array([D[n] for n in names]); fd = np.array([fD[n] for n in names])
def xyz(ra, de, d): return np.c_[d*np.cos(de)*np.cos(ra), d*np.cos(de)*np.sin(ra), d*np.sin(de)]
def stats(P):
    t = Delaunay(P); s = t.simplices
    a, b, c, e = P[s[:,0]], P[s[:,1]], P[s[:,2]], P[s[:,3]]
    V = np.sort(np.abs(np.einsum('ij,ij->i', b-a, np.cross(c-a, e-a)))/6); n = len(V)
    nn = cKDTree(P).query(P, k=2)[0][:, 1]
    return dict(median_nn=np.median(nn), small_frac=np.mean(V < np.median(V)*0.01),
                logV_std=np.std(np.log10(V)), top10_share=V[int(0.9*n):].sum()/V.sum())
# predicted direction for each ruler: +1 = real should be LARGER than chance, -1 = SMALLER
direction = dict(median_nn=-1, small_frac=+1, logV_std=+1, top10_share=+1)
rng = np.random.default_rng(seed); out = {}
cuts = [("all 175", np.ones(len(names), bool)), ("within 30 Mpc", d <= 30),
        ("excluding Ursa Major", fd != 4)]
for label, mask in cuts:
    P = xyz(ra[mask], de[mask], d[mask]); real = stats(P)
    nulls = [stats(xyz(ra[mask], de[mask], rng.permutation(d[mask]))) for _ in range(N)]
    out[label] = {}
    print(f"\n== {label}: {mask.sum()} galaxies, {N} shuffles")
    for k, sgn in direction.items():
        arr = np.array([x[k] for x in nulls]); z = (real[k]-arr.mean())/arr.std()
        m = int(np.sum(sgn*arr >= sgn*real[k]))          # shuffles at least as extreme
        p = (m + 1)/(N + 1)                               # standard empirical p (never zero)
        zp = norm.isf(p)                                  # the z this p corresponds to
        out[label][k] = dict(real=real[k], null_mean=arr.mean(), null_sd=arr.std(), z_gauss=z,
                             n_as_extreme=m, p_empirical=p, z_from_p=zp, null_skew=float(skew(arr)))
        print(f"{k:12s} z(bell-curve)={z:+.2f}  shuffles as extreme={m:5d}/{N}  p={p:.2g}  z(from count)={zp:.2f}  skew={skew(arr):+.2f}")
json.dump(out, open(f"ranks_{N}_seed{seed}.json", "w"), indent=1)
