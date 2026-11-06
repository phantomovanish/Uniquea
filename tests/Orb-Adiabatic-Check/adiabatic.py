import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
A2=1.0; eps=0.004344; A4=-(A2/6)*(1-eps); T1=2*A2+4*A4
dv=lambda y,m: -2*A2*y**-3-4*A4*y**-5+T1*y**-2+2*m*y          # I ∝ chi^2, magnetic term m*y^2
v =lambda y,m: A2*y**-2+A4*y**-4-T1/y+m*y*y
def eq(m):
    r=np.roots([2*m,0,0,T1,-2*A2,0,-4*A4]) if m!=0 else np.roots([T1,-2*A2,0,-4*A4])
    return sorted(x.real for x in r if abs(x.imag)<1e-9 and 0.9<x.real<1.1)
lo,hi=0,1e-4
for _ in range(100):
    mid=(lo+hi)/2; (lo,hi)=(mid,hi) if len(eq(mid))==2 else (lo,mid)
mc=lo; print("m_c =",mc)
lip0,bot0=eq(0); K0=(dv(bot0+1e-7,0)-dv(bot0-1e-7,0))/2e-7; w0=np.sqrt(K0)
print("cup: bottom",bot0,"lip",lip0,"omega0",w0,"period",2*np.pi/w0,"barrier",v(lip0,0)-v(bot0,0))
def run(E0,T):
    rate=mc/T
    y0=bot0; p0=np.sqrt(2*E0)                    # start at bottom with kinetic energy E0
    f=lambda t,s:[s[1],-dv(s[0],rate*t)]
    ev=lambda t,s: s[0]-0.985; ev.terminal=True; ev.direction=-1
    sol=solve_ivp(f,[0,3*T],[y0,p0],method="DOP853",rtol=1e-11,atol=1e-14,events=ev,max_step=2.0)
    te=sol.t_events[0][0]; return rate*te
print("\n E0        T(ramp)   (m_c - m_esc)/m_c")
res={}
for T in []:
    for E0 in [1e-9,3e-10,1e-10,3e-11,1e-11]:
        me=run(E0,T); res[(T,E0)]=mc-me; print(f"{E0:.0e}  {T:.0e}   {(mc-me)/mc:+.4e}")
for T in []:
    Es=[1e-9,3e-10,1e-10,3e-11,1e-11]; d=[res[(T,E)] for E in Es]
    if min(d)>0: print("T=%.0e slope of log(m_c-m_esc) vs log J: %.3f"%(T,np.polyfit(np.log(Es),np.log(d),1)[0]))
print("slower ramp")
for E0 in [1e-10,3e-11,1e-11]:
    me=run(E0,5e6); print(f"{E0:.0e}  5e+06   {(mc-me)/mc:+.4e}")
