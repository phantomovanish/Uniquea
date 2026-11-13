"""Test 3: Phantom's tetrahedra in the real neighborhood (Uniquea, November 12, 2026).

Prediction (stated before the run): connect the SPARC galaxies in 3D and the tetrahedra fit snugly:
small ones crowd into groups and filaments, big ones span the voids between.

Data: SPARC Table1 (Lelli, McGaugh & Schombert 2016, AJ 152, 157), official file SPARC_Lelli2016c.mrt,
for distances; sky positions for all 175 galaxies from NASA/IPAC NED (sparc_sky_coordinates_ned.csv).
Method: 3D Delaunay tetrahedralization. Null: 300 shuffles keeping every sky position and the set of
distances, scrambling which distance goes with which galaxy (destroys 3D clustering, keeps survey shape).
Usage: python3 tetrahedra_sparc.py <SPARC_Lelli2016c.mrt> <sparc_sky_coordinates_ned.csv>
"""
import sys, numpy as np
from scipy.spatial import Delaunay, cKDTree
tab, coordfile = sys.argv[1:3]
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
    t = Delaunay(P); V = []
    for s in t.simplices:
        a, b, c, e = P[s]; V.append(abs(np.dot(b-a, np.cross(c-a, e-a)))/6)
    V = np.sort(np.array(V)); n = len(V)
    nn = cKDTree(P).query(P, k=2)[0][:, 1]
    return dict(n_tet=n, median_nn=np.median(nn), small_frac=np.mean(V < np.median(V)*0.01),
                logV_std=np.std(np.log10(V)), top10_share=V[int(0.9*n):].sum()/V.sum())
rng = np.random.default_rng(42)
def run(mask, label):
    P = xyz(ra[mask], de[mask], d[mask]); real = stats(P)
    nulls = [stats(xyz(ra[mask], de[mask], rng.permutation(d[mask]))) for _ in range(300)]
    print(f"\n== {label}: {mask.sum()} galaxies")
    for k in real:
        arr = np.array([x[k] for x in nulls]); z = (real[k]-arr.mean())/arr.std()
        print(f"{k:12s} real={real[k]:.4g}  null={arr.mean():.4g}±{arr.std():.3g}  z={z:+.2f}")
run(np.ones(len(names), bool), "all 175")
run(d <= 30, "within 30 Mpc")
run(fd != 4, "excluding Ursa Major (assigned 18 Mpc)")
