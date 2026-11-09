"""Phantom's decisive calculation: Architecture B (Master 3.1, parameters B34) with a defined Lila-damping
term on the z-mode, its lost energy booked into the Alt Field, the cycle and Floquet multipliers recomputed.
Lila term: dz/dtau gains -gamma_L z with gamma_L = zeta * Omega (amplitude damping rate of an oscillator of
frequency Omega with damping ratio zeta = eps_L). Ledger for E_z = r^2/2:
  dE_z/dtau = alpha0 r^2 - [beta_LC r^4 + Gamma(1+beta_U r^2) r^2 + gamma_L r^2]
  dE_Alt/dtau = -alpha0 r^2 + [same losses]          so E_z + E_Alt is conserved.
Usage: python3 lila_on_z.py"""
import numpy as np
from scipy.integrate import solve_ivp
bU,bLC,a0,G0,Om0,Om1,kJ,gq,kq,e1,e2=0.5,1,1,0.1,1,0.5,0.5,1,0.2,0.2,0.1
EPS_L=0.0043444445555666779
def rhs(t,s,zeta,book=False):
    J,q,z1,z2=s[:4]; r2=z1*z1+z2*z2; B8=abs(J-1)**8*(1+bU*r2)**4; G=G0*B8*(1+bU*r2); Om=Om0+Om1*B8
    gL=zeta*Om; A=a0-bLC*r2-G*(1+bU*r2)-gL
    d=[-G*(1+bU*r2)**2*r2-kJ*(J-1), -gq*q+kq*(J-1)+e1*z1+e2*z2, A*z1-Om*z2, A*z2+Om*z1]
    if book:
        loss=bLC*r2*r2+G*(1+bU*r2)*r2+gL*r2; d+= [-a0*r2+loss, gL*r2]   # E_Alt, and the Lila share alone
    return d
def cycle(zeta):
    s=solve_ivp(rhs,[0,400],[1.2,0.1,0.6,0.0],args=(zeta,),rtol=1e-12,atol=1e-14).y[:,-1]
    T=2*np.pi                      # on the orbit J=1, B8=0, Omega=Om0=1
    return s,T
def floquet(zeta,s,T,h=1e-7):
    M=np.zeros((4,4))
    for i in range(4):
        dp=np.zeros(4); dp[i]=h
        a=solve_ivp(rhs,[0,T],s+dp,args=(zeta,),rtol=1e-12,atol=1e-14).y[:,-1]
        b=solve_ivp(rhs,[0,T],s-dp,args=(zeta,),rtol=1e-12,atol=1e-14).y[:,-1]
        M[:,i]=(a-b)/(2*h)
    return np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
for name,zeta in [("No Lila damping (archive B34)",0.0),("Lila damping, zeta = eps_L",EPS_L),("Extinction test, zeta = 1.2",1.2)]:
    s,T=cycle(zeta); r=np.hypot(s[2],s[3])
    mult=floquet(zeta,s,T) if r>1e-6 else None
    full=solve_ivp(rhs,[0,T],list(s)+[0,0],args=(zeta,True),rtol=1e-12,atol=1e-14).y[:,-1]
    Ez0=r*r/2; Ez1=(full[2]**2+full[3]**2)/2; ledger=(Ez1+full[4])-Ez0
    print(f"{name}")
    print(f"   cycle radius r* = {r:.6f}  (predicted sqrt(1 - zeta) = {np.sqrt(max(1-zeta,0)):.6f})   J = {s[0]:.6f}")
    if mult is not None:
        print("   Floquet multipliers:", ", ".join(f"{m:.6g}" for m in mult))
        print(f"   Lila loss per cycle = {full[5]:.5f} = {100*full[5]/Ez0:.2f}% of E_z; injection balances it; ledger residual over one cycle {ledger:+.1e}")
    else:
        print("   the breath is extinguished: r -> 0")
