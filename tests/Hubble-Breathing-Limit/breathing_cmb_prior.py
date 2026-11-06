"""Long-period wobbles with a CMB prior on Ωm (Planck 2018, TT,TE,EE+lowE+lensing: Ωm = 0.3153 ± 0.0073).
Same model and data as breathing_limit.py; adds ((Ωm − 0.3153)/0.0073)² to the joint χ².
Caveat: Planck's Ωm assumes ΛCDM at all redshifts; this is the simplest stand-in for a full CMB likelihood.
Usage: python3 breathing_cmb_prior.py <data_dir> <out.json>"""
import sys, json, numpy as np
out_path=sys.argv[2]; sys.argv=[sys.argv[0],sys.argv[1],"/dev/null"]
sys.path.insert(0,__import__("os").path.dirname(__file__))
import breathing_limit as b
from scipy.optimize import minimize
from multiprocessing import Pool
OM,SIG=0.3153,0.0073
def chi2p(Om,eps,P,phi): return b.chi2(Om,eps,P,phi)+((Om-OM)/SIG)**2
def prof(args):
    eps,P=args; best=None
    for phi0 in np.linspace(0,2*np.pi,8,endpoint=False):
        r=minimize(lambda q:chi2p(q[0],eps,P,q[1]),[OM,phi0],method="Nelder-Mead",options={"xatol":1e-4,"fatol":1e-3,"maxiter":500})
        if best is None or r.fun<best[0]: best=(r.fun,r.x[0],r.x[1])
    return dict(eps=eps,P=P,chi2=best[0],Om=best[1],phi=best[2]%(2*np.pi))
if __name__=="__main__":
    r0=minimize(lambda q:chi2p(q[0],0,1,0),[OM],method="Nelder-Mead"); c0=float(r0.fun)
    Ps=[1.0,2.0,3.0,5.0,8.0]; E=[0,0.01,0.02,0.03,0.04,0.06,0.08,0.12,0.16,0.2,0.3,0.4,0.5]
    with Pool(2) as Pl: res=Pl.map(prof,[(e,P) for P in Ps for e in E])
    out={"prior":[OM,SIG],"LCDM_chi2":c0,"LCDM_Om":float(r0.x[0]),"grid":res,"limits":{}}
    for P in Ps:
        rr=sorted([r for r in res if r["P"]==P],key=lambda r:r["eps"]); ee=np.array([r["eps"] for r in rr])
        ch=np.array([min(r["chi2"],c0) if r["eps"]==0 else r["chi2"] for r in rr]); cm=ch.min(); ib=int(ch.argmin()); lim=None
        for k in range(ib,len(ee)-1):
            if ch[k]-cm<3.84<=ch[k+1]-cm: lim=ee[k]+(3.84-(ch[k]-cm))*(ee[k+1]-ee[k])/((ch[k+1]-cm)-(ch[k]-cm)); break
        oms=[round(r["Om"],3) for r in rr]
        out["limits"][str(P)]=dict(eps95=lim,best_eps=float(ee[ib]),dchi2_best_vs_LCDM=float(cm-c0),Om_range=[min(oms),max(oms)])
        print(f"P={P:<4} best ε={ee[ib]:.3f} Δχ²={cm-c0:+.2f} 95% limit ε<{lim} Ωm range {min(oms)}–{max(oms)}",flush=True)
    print("ΛCDM+prior χ²=%.2f Ωm=%.4f"%(c0,r0.x[0])); json.dump(out,open(out_path,"w"),indent=1)
