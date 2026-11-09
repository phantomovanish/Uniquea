"""Plot g_obs vs g_bar for the SPARC sample with the rung-zero deep law and the RAR fit.
Usage: python3 plot_rar.py <Rotmod_LTG dir> <SPARC_Lelli2016c.mrt> <results.json> <out.png>"""
import sys, os, re, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
rot,tab,res,out=sys.argv[1:5]; r=json.load(open(res)); KPC=3.0857e19
meta={}
for l in open(tab,encoding="utf-8",errors="ignore"):
    p=l.split()
    if len(p)>=18 and re.match(r"^\d+$",p[1]):
        try: meta[p[0]]=dict(inc=float(p[5]),Q=int(p[17]))
        except: pass
B=[];O=[]
for f in os.listdir(rot):
    g=f.replace("_rotmod.dat",""); m=meta.get(g)
    if not f.endswith("_rotmod.dat") or m is None or m["Q"]>=3 or m["inc"]<30: continue
    R,Vo,eV,Vg,Vd,Vb=np.loadtxt(os.path.join(rot,f),comments="#",usecols=range(6),ndmin=2).T
    ok=(R>0)&(Vo>0)&(eV>0); R,Vo,Vg,Vd,Vb=[x[ok] for x in (R,Vo,Vg,Vd,Vb)]
    gb=(Vg*np.abs(Vg)+0.5*Vd**2+0.7*Vb**2)*1e6/(R*KPC); go=Vo**2*1e6/(R*KPC); k=gb>0; B+=list(gb[k]); O+=list(go[k])
B=np.array(B);O=np.array(O); x=np.logspace(-13,-8,300)
fig,ax=plt.subplots(figsize=(7.2,6.2),dpi=170)
ax.loglog(B,O,".",ms=2.2,color="#8AA0C8",alpha=.55,label=f"SPARC: {r['points']} points, {r['galaxies']} galaxies")
ax.loglog(x,x,"--",color="#555",lw=1,label="Newton, no extra (g_obs = g_bar)")
ax.loglog(x,np.sqrt(r["a0_deep"]*x),"-",color="#B4552D",lw=2,label=f"Rung zero (deep): g = √(a₀ g_bar), a₀ = {r['a0_deep']*1e10:.2f}×10⁻¹⁰")
ax.loglog(x,x/(1-np.exp(-np.sqrt(x/r["a0_rar"]))),"-",color="#1F3864",lw=2,label=f"RAR function, a₀ = {r['a0_rar']*1e10:.2f}×10⁻¹⁰ m/s²")
ax.axvline(0.12e-10,color="#999",lw=.8,ls=":"); ax.text(0.13e-10,3e-9,"deep regime ←",fontsize=8,color="#666",ha="left")
ax.set_xlim(1e-13,1e-8); ax.set_ylim(1e-12,1e-8)
ax.set_xlabel("Acceleration from visible matter, g_bar (m/s²)"); ax.set_ylabel("Observed acceleration, g_obs (m/s²)")
ax.set_title("Uniquea test 1: one a₀ against 153 galaxies (SPARC)",fontsize=11)
ax.legend(fontsize=7.8,loc="upper left",frameon=False); ax.grid(alpha=.25,which="both")
fig.text(0.01,0.005,"Data: Lelli, McGaugh & Schombert 2016 (SPARC). Υ_disk = 0.5, Υ_bulge = 0.7; Q < 3, i ≥ 30°.",fontsize=6.8,color="#555")
plt.tight_layout(rect=(0,0.02,1,1)); plt.savefig(out)
