"""Hubble-closure test: the sine closure (Master §3.2, B44) vs flat ΛCDM.

Reproduce:
  git clone https://github.com/CobayaSampler/bao_data.git            (commit bb0c1c9)
  git clone --filter=blob:none --sparse https://github.com/PantheonPlusSH0ES/DataRelease.git pp  (commit c447f0f)
  cd pp && git sparse-checkout set "Pantheon+_Data/4_DISTANCES_AND_COVAR" && cd ..
  python3 closure_test.py <data_dir>

Models (E = H/H0):
  ΛCDM (flat):   E(z) = sqrt(Ωm (1+z)^3 + 1 − Ωm)                          1 free parameter
  Sine closure:  cos θ = cos θ0 + ln(1+z)/A,  E(z) = sin θ / sin θ0         2 free parameters (A, θ0)
Likelihoods (Gaussian, χ² = rᵀ C⁻¹ r):
  DESI DR2 BAO: 13 values of DM/rd, DH/rd, DV/rd with the published full covariance.
                The overall scale c/(H0 rd) is fitted analytically (profiled).
  Pantheon+:    MU_SH0ES with the STAT+SYS covariance. Selection zHD > 0.01 and IS_CALIBRATOR == 0
                → 1,580 light curves of 1,466 distinct supernovae (1,590 light curves if the 10
                calibrator light curves at z > 0.01 are kept, the usual Pantheon+-only cut).
                The absolute-magnitude/H0 offset is fitted analytically (profiled).
Optimizer: Nelder–Mead from several starts; a 161 × 400 grid over (A, θ0) confirms the BAO minimum.
Limitations: best-fit χ² only (no posterior, no parameter errors); inclination and radiation ignored;
the sine best fit lies on the boundary A → ∞, θ0 → π, so its χ² is a limiting value.
"""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid
D=sys.argv[1].rstrip("/")+"/" if len(sys.argv)>1 else "./"
def E_lcdm(z,p): Om=p[0]; return np.sqrt(Om*(1+z)**3+1-Om)
def E_sine(z,p):
    A,t0=p; c=np.cos(t0)+np.log1p(z)/A
    if np.any(np.abs(c)>=1) or np.sin(t0)<=0: return None
    return np.sqrt(1-c**2)/np.sin(t0)
zg=np.linspace(0,2.5,5001)
def chi_int(E,p):
    e=E(zg,p)
    if e is None or np.any(e<=0): return None
    return cumulative_trapezoid(1/e,zg,initial=0)
m=pd.read_csv(D+"bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",sep=r"\s+",comment="#",names=["z","v","q"])
Ci=np.linalg.inv(np.loadtxt(D+"bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt"))
def chi2_bao(E,p):
    I=chi_int(E,p)
    if I is None: return 1e12
    t=[]
    for z,q in zip(m.z,m.q):
        dm=np.interp(z,zg,I); dh=1/E(np.array([z]),p)[0]
        t.append({"DM_over_rs":dm,"DH_over_rs":dh,"DV_over_rs":(z*dm*dm*dh)**(1/3)}[q])
    t=np.array(t); v=m.v.values; B=(t@Ci@v)/(t@Ci@t); r=v-B*t; return float(r@Ci@r)
pp=D+"pp/Pantheon+_Data/4_DISTANCES_AND_COVAR/"
s=pd.read_csv(pp+"Pantheon+SH0ES.dat",sep=r"\s+"); N=len(s)
cov=np.loadtxt(pp+"Pantheon+SH0ES_STAT+SYS.cov",skiprows=1).reshape(N,N)
keep=((s.zHD>0.01)&(s.IS_CALIBRATOR==0)).values
s=s[keep]; Si=np.linalg.inv(cov[np.ix_(keep,keep)]); mu=s.MU_SH0ES.values; one=np.ones(len(mu))
def chi2_sn(E,p):
    I=chi_int(E,p)
    if I is None: return 1e12
    r=mu-5*np.log10((1+s.zHEL.values)*np.interp(s.zHD.values,zg,I))
    r=r-(one@Si@r)/(one@Si@one); return float(r@Si@r)
def fit(E,chi2,starts):
    b=None
    for x0 in starts:
        r=minimize(lambda p:chi2(E,p),x0,method="Nelder-Mead",options={"xatol":1e-6,"fatol":1e-6,"maxiter":20000})
        if b is None or r.fun<b.fun: b=r
    return b
out={"N_light_curves":int(len(mu)),"N_distinct_SNe":int(s.CID.nunique()),"N_BAO":int(len(m))}
for name,chi in [("DESI_DR2_BAO",chi2_bao),("PantheonPlus",chi2_sn)]:
    L=fit(E_lcdm,chi,[[0.3]]); S=fit(E_sine,chi,[[A,t] for A in [1.5,3,6,12,30] for t in [1.7,2.2,2.7,3.0]])
    out[name]={"LCDM_chi2":round(L.fun,2),"Om":round(float(L.x[0]),4),"sine_chi2":round(S.fun,2),
               "sine_A":float(S.x[0]),"sine_theta0":float(S.x[1]),"delta_chi2":round(S.fun-L.fun,2)}
    print(name,out[name],flush=True)
json.dump(out,open("closure_results.json","w"),indent=1); print(out)
