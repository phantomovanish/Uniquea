"""Uniquea test 1: does one a0 fit galaxy rotation curves (rung zero = deep-MOND, F = m a^2/a0)?
Data: SPARC (Lelli, McGaugh & Schombert 2016), 175 galaxies, Rotmod_LTG files + SPARC_Lelli2016c table.
Baryonic acceleration g_bar = (Vgas|Vgas| + Yd Vdisk^2 + Yb Vbul^2)/R, Yd = 0.5, Yb = 0.7 (3.6 micron).
Cuts as McGaugh, Lelli & Schombert (2016): quality Q < 3, inclination >= 30 deg.
Tests:
 (1) Deep regime only (g_bar < 0.1 * 1.2e-10): rung zero predicts g_obs = sqrt(a0 g_bar); fit one a0.
 (2) All points with the standard RAR function g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0))) (context; MOND's own).
 (3) Universality: best a0 per galaxy (deep points) and its spread vs. the measurement errors.
Usage: python3 rung_zero_sparc.py <Rotmod_LTG dir> <SPARC_Lelli2016c.mrt> <out.json> [published]
  Add "published" to also drop points with velocity uncertainty >= 10% (McGaugh et al. 2016 cut)."""
import sys, os, re, json, numpy as np
from scipy.optimize import minimize_scalar
rot, tab, out = sys.argv[1:4]; PUB = len(sys.argv)>4 and sys.argv[4]=="published"
KPC=3.0857e19; Yd,Yb=0.5,0.7
meta={}
for l in open(tab,encoding="utf-8",errors="ignore"):
    p=l.split()
    if len(p)>=18 and re.match(r"^\d+$",p[1]) and p[0] not in ("Byte",):
        try: meta[p[0]]=dict(D=float(p[2]),inc=float(p[5]),Q=int(p[17]))
        except: pass
pts=[]; per={}
for f in sorted(os.listdir(rot)):
    if not f.endswith("_rotmod.dat"): continue
    g=f.replace("_rotmod.dat","")
    m=meta.get(g)
    if m is None or m["Q"]>=3 or m["inc"]<30: continue
    a=np.loadtxt(os.path.join(rot,f),comments="#",usecols=range(6),ndmin=2)
    R,Vo,eV,Vg,Vd,Vb=a.T
    ok=(R>0)&(Vo>0)&(eV>0)
    if PUB: ok&=(eV/np.where(Vo>0,Vo,1))<0.1
    R,Vo,eV,Vg,Vd,Vb=[x[ok] for x in (R,Vo,eV,Vg,Vd,Vb)]
    gb=(Vg*np.abs(Vg)+Yd*Vd**2+Yb*Vb**2)*1e6/(R*KPC); go=Vo**2*1e6/(R*KPC)
    eg=2*go*eV/Vo
    keep=gb>0
    for x in zip(gb[keep],go[keep],eg[keep]): pts.append(x)
    per[g]=np.c_[gb[keep],go[keep],eg[keep]]
P=np.array(pts); gb,go,eg=P.T
# log-space residuals with fractional errors (plus 0.1 dex floor for distance/inclination/M/L systematics)
def chi(model,mask): 
    s=np.sqrt((eg[mask]/go[mask]/np.log(10))**2+0.1**2)
    return np.sum(((np.log10(go[mask])-np.log10(model[mask]))/s)**2)
deep=gb<0.1*1.2e-10
f1=lambda la: chi(np.sqrt(10**la*gb),deep)
r1=minimize_scalar(f1,bounds=(-11,-9),method="bounded")
rar=lambda la: gb/(1-np.exp(-np.sqrt(gb/10**la)))
f2=lambda la: chi(rar(la),np.ones_like(gb,bool))
r2=minimize_scalar(f2,bounds=(-11,-9),method="bounded")
# 1-sigma (delta chi2 = 1) on log a0
def err(f,x0):
    xs=np.linspace(x0-0.2,x0+0.2,4001); c=np.array([f(x) for x in xs]); inside=xs[c<=f(x0)+1]; return (inside.max()-inside.min())/2
res_deep=np.log10(go[deep])-np.log10(np.sqrt(10**r1.x*gb[deep]))
# per-galaxy a0 from deep points (galaxies with >= 5 deep points)
pg=[]
for g,a in per.items():
    d=a[:,0]<0.1*1.2e-10
    if d.sum()<5: continue
    b,o,e=a[d].T; s=np.sqrt((e/o/np.log(10))**2+0.1**2)
    f=lambda la: np.sum(((np.log10(o)-0.5*(la+np.log10(b)))/s)**2)
    r=minimize_scalar(f,bounds=(-12,-8),method="bounded"); pg.append((g,r.x))
la=np.array([x[1] for x in pg])
outd=dict(galaxies=len(per),points=len(gb),deep_points=int(deep.sum()),
  a0_deep=10**r1.x, a0_deep_logerr=err(f1,r1.x), chi2_deep=float(r1.fun), dof_deep=int(deep.sum()-1),
  deep_resid_rms_dex=float(np.std(res_deep)), deep_resid_mean_dex=float(np.mean(res_deep)),
  a0_rar=10**r2.x, a0_rar_logerr=err(f2,r2.x), chi2_rar=float(r2.fun), dof_rar=int(len(gb)-1),
  per_galaxy_n=len(pg), per_galaxy_log_a0_mean=float(la.mean()), per_galaxy_log_a0_std=float(la.std()),
  cH0_over_2pi=2.998e8*(67.4e3/3.0857e22)/(2*np.pi))
json.dump(outd,open(out,"w"),indent=1)
for k,v in outd.items(): print(f"{k:28s} {v:.4g}" if isinstance(v,float) else f"{k:28s} {v}")
