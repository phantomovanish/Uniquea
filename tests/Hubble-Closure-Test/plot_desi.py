import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import importlib.util,sys
src=open("hz/realfit.py").read().split("# ---------- Pantheon+")[0]
exec(src)
from scipy.optimize import minimize
L=minimize(lambda p:chi2_bao(E_lcdm,p),[0.3],method="Nelder-Mead").x
best=None
for A in [3,30,300,3e3,3e4,3e5]:
  for t in [2.5,3.0,3.1,3.13]:
    r=minimize(lambda p:chi2_bao(E_sine,p),[A,t],method="Nelder-Mead",options={"maxiter":4000})
    if best is None or r.fun<best.fun: best=r
S=best.x
v=m.v.values; sd=np.sqrt(np.diag(C))
def scaled(E,p):
    t=bao_pred(E,p); B=(t@Ci@v)/(t@Ci@t); return B
BL=scaled(E_lcdm,L); BS=scaled(E_sine,S)
z=np.linspace(0.05,2.45,300); h=m.q=="DH_over_rs"
fig,ax=plt.subplots(figsize=(8,5),dpi=150)
ax.errorbar(m.z[h],v[h],sd[h],fmt="o",color="#111",capsize=3,label="DESI DR2 BAO (measured)")
ax.plot(z,BL/E_lcdm(z,L),color="#1f6feb",lw=2,label=f"ΛCDM (Ωm = {L[0]:.3f}), χ² = {chi2_bao(E_lcdm,L):.0f}")
ax.plot(z,BS/E_sine(z,S),color="#d4a017",lw=2,ls="--",label=f"Sine closure, best fit, χ² = {best.fun:.0f}")
ax.set_xlabel("Redshift z (looking back in time →)"); ax.set_ylabel("D_H / r_d  =  c / (H r_d)")
ax.set_title("The sine breath is too flat: DESI DR2 expansion-rate data",fontsize=11)
ax.legend(frameon=False); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig("/mnt/user-data/outputs/Hubble-Closure-vs-DESI.png")
print(L,S,best.fun)
