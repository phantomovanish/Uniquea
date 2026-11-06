"""Upper limit on how much the universe can still be 'breathing' today (Uniquea fork 7, idea 2).

Model (E = H/H0), a small wobble on flat ΛCDM, periodic in x = ln(1+z) (log of the scale factor):
  E(z) = E_LCDM(z; Ωm) · [1 + ε·(sin(2πx/P + φ) − sin φ)]      (E(0) = 1)
  ε = fractional wobble amplitude in H; P = period in ln(1+z); φ = phase.
Data and likelihoods exactly as closure_test.py (DESI DR2 BAO full covariance, scale profiled;
Pantheon+ STAT+SYS, zHD > 0.01, calibrators excluded, offset profiled). Joint χ² = BAO + SN.
For each P: profile over (Ωm, φ) on a grid of ε; 95% limit where χ²(ε) − min χ² = 3.84.
Usage: python3 breathing_limit.py <data_dir> <out.json>
"""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid
from multiprocessing import Pool
D=sys.argv[1].rstrip("/")+"/"
zg=np.linspace(0,2.5,5001); xg=np.log1p(zg)
m=pd.read_csv(D+"bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",sep=r"\s+",comment="#",names=["z","v","q"])
Ci=np.linalg.inv(np.loadtxt(D+"bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt"))
pp=D+"pp/Pantheon+_Data/4_DISTANCES_AND_COVAR/"
s=pd.read_csv(pp+"Pantheon+SH0ES.dat",sep=r"\s+"); N=len(s)
cov=np.loadtxt(pp+"Pantheon+SH0ES_STAT+SYS.cov",skiprows=1).reshape(N,N)
keep=((s.zHD>0.01)&(s.IS_CALIBRATOR==0)).values
s=s[keep]; Si=np.linalg.inv(cov[np.ix_(keep,keep)]); mu=s.MU_SH0ES.values; one=np.ones(len(mu))
zHD=s.zHD.values; zHEL=s.zHEL.values; oSo=one@Si@one
def Egrid(Om,eps,P,phi):
    return np.sqrt(Om*(1+zg)**3+1-Om)*(1+eps*(np.sin(2*np.pi*xg/P+phi)-np.sin(phi)))
def chi2(Om,eps,P,phi):
    if not (0.05<Om<0.6): return 1e12
    e=Egrid(Om,eps,P,phi)
    if np.any(e<=0): return 1e12
    I=cumulative_trapezoid(1/e,zg,initial=0)
    t=[]
    for z,q in zip(m.z,m.q):
        dm=np.interp(z,zg,I); dh=1/np.interp(z,zg,e)
        t.append({"DM_over_rs":dm,"DH_over_rs":dh,"DV_over_rs":(z*dm*dm*dh)**(1/3)}[q])
    t=np.array(t); v=m.v.values; B=(t@Ci@v)/(t@Ci@t); r=v-B*t; cb=float(r@Ci@r)
    r=mu-5*np.log10((1+zHEL)*np.interp(zHD,zg,I)); r=r-(one@Si@r)/oSo; cs=float(r@Si@r)
    return cb+cs
def prof(args):
    eps,P=args; best=None
    for phi0 in np.linspace(0,2*np.pi,6,endpoint=False):
        r=minimize(lambda q:chi2(q[0],eps,P,q[1]),[0.31,phi0],method="Nelder-Mead",options={"xatol":1e-4,"fatol":1e-3,"maxiter":400})
        if best is None or r.fun<best[0]: best=(r.fun,r.x[0],r.x[1])
    return dict(eps=eps,P=P,chi2=best[0],Om=best[1],phi=best[2]%(2*np.pi))
if __name__=="__main__":
    r0=minimize(lambda q:chi2(q[0],0,1,0),[0.31],method="Nelder-Mead"); c0=float(r0.fun)
    Ps=[0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.5,2.0,3.0,5.0]
    E=[0,0.0025,0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.06,0.08,0.12]
    with Pool(2) as Pl: res=Pl.map(prof,[(e,P) for P in Ps for e in E])
    out={"LCDM_joint_chi2":c0,"LCDM_Om":float(r0.x[0]),"grid":res,"limits":{}}
    for P in Ps:
        rr=sorted([r for r in res if r["P"]==P],key=lambda r:r["eps"])
        ch=np.array([min(r["chi2"],c0) if r["eps"]==0 else r["chi2"] for r in rr]); ee=np.array([r["eps"] for r in rr])
        cmin=ch.min(); ib=int(ch.argmin()); lim=None
        for k in range(ib,len(ee)-1):
            if ch[k]-cmin<3.84<=ch[k+1]-cmin:
                lim=ee[k]+(3.84-(ch[k]-cmin))*(ee[k+1]-ee[k])/((ch[k+1]-cmin)-(ch[k]-cmin)); break
        out["limits"][str(P)]=dict(eps95=lim,best_eps=float(ee[ib]),dchi2_best_vs_LCDM=float(cmin-c0))
        print(f"P={P:<5} best ε={ee[ib]:.4f} Δχ²(best−ΛCDM)={cmin-c0:+.2f}  95% limit ε<{lim}",flush=True)
    print("ΛCDM joint χ²=%.2f Ωm=%.4f"%(c0,r0.x[0]))
    json.dump(out,open(sys.argv[2],"w"),indent=1)
