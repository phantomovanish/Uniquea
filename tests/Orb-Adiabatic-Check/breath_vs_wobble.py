"""Damp the wobble without damping away the breath (Stuart-Landau, the archive's limit-cycle core, 3.1).
dz/dt = (alpha0 - gamma + i*omega) z - beta |z|^2 z
  alpha0: injection drawn from the Alt Field; gamma: Lila-type loss to the Alt Field; beta: saturation loss.
Breath = the cycle (|z| = r*); wobble = a deviation of |z| from r*.  Ledger: d(E + E_Alt)/dt = 0 with E = |z|^2/2."""
import numpy as np
from scipy.integrate import solve_ivp
omega=2*np.pi; zeta=0.004344; gamma=2*zeta*omega; beta=1.0
def run(alpha0, r0, T=400.0):
    def f(t,s):
        x,y,Ealt=s; r2=x*x+y*y; mu=alpha0-gamma
        dx=mu*x-omega*y-beta*r2*x; dy=mu*y+omega*x-beta*r2*y
        dEalt=(gamma*r2+beta*r2*r2)-alpha0*r2        # losses in, injection out
        return [dx,dy,dEalt]
    s=solve_ivp(f,[0,T],[r0,0,0],rtol=1e-10,atol=1e-13,dense_output=True)
    x,y,Ea=s.y[:,-1]; r=np.hypot(x,y); E0=r0**2/2; E=r**2/2
    return r, E+Ea-E0
rstar=lambda a: np.sqrt(max(a-gamma,0)/beta)
for label,a in [("alpha0 = 2 x loss (sustained)",2*gamma),("alpha0 = 1.2 x loss (sustained, small)",1.2*gamma),("alpha0 = 0.5 x loss (breath dies)",0.5*gamma)]:
    for r0 in [0.5*max(rstar(a),0.05),1.5*max(rstar(a),0.05)]:
        r,led=run(a,r0)
        print(f"{label:40s} start r={r0:.4f} -> after 400 cycles r={r:.4f} (predicted r*={rstar(a):.4f}); ledger residual {led:+.1e}")
a=2*gamma; print("Transverse (wobble) decay per cycle, exp(-2(alpha0-gamma)T) =",np.exp(-2*(a-gamma)*1.0))
print("Archive multiplier 0.0432 implies (alpha0-gamma)*T =",np.log(1/0.0432)/2,"(pi/2 =",np.pi/2,")")
