# Option A (no relaxation: wobble action J carried, family of outcomes) vs
# Option B (weak damping gamma lets J relax; lost energy booked to the Alt Field).
import numpy as np, sys, json
from scipy.integrate import solve_ivp
from multiprocessing import Pool
A2=1.0; eps=0.004344; A4=-(A2/6)*(1-eps); T1=2*A2+4*A4
dv=lambda y,m: -2*A2*y**-3-4*A4*y**-5+T1*y**-2+2*m*y
v =lambda y,m: A2*y**-2+A4*y**-4-T1/y+m*y*y
mc=4.765797337118178e-06; bot0=1.0; T=1e6; rate=mc/T
def run(args):
    E0,gT=args; g=gT/T
    # state: y, p, Wramp (work done by contraction), Ealt (energy sent to Alt Field)
    def f(t,s):
        y,p,W,Ea=s
        return [p, -dv(y,rate*t)-g*p, rate*y*y, g*p*p]
    ev=lambda t,s: s[0]-0.985; ev.terminal=True; ev.direction=-1
    s0=[bot0,np.sqrt(2*E0),0.0,0.0]
    sol=solve_ivp(f,[0,3*T],s0,method="DOP853",rtol=1e-11,atol=1e-15,events=ev,max_step=2.0,dense_output=False)
    te=sol.t_events[0][0]; ye=sol.y_events[0][0]; y,p,W,Ea=ye
    E_start=0.5*s0[1]**2+v(bot0,0); E_now=0.5*p*p+v(y,rate*te)
    ledger=E_now-E_start-W+Ea          # should be 0: change = work by contraction - loss to Alt Field
    return dict(E0=E0,gT=gT,shift=(mc-rate*te)/mc,Ealt=Ea,ledger=ledger)
if __name__=="__main__":
    jobs=[(E0,gT) for E0 in [1e-9,3e-10,1e-10] for gT in [0,1,3,10]]
    with Pool(len(jobs)) as P: out=P.map(run,jobs)
    json.dump(out,open(sys.argv[1],"w"),indent=1)
    for r in out: print(f"E0={r['E0']:.0e} gammaT={r['gT']:>4} shift={r['shift']:+.4e} Ealt={r['Ealt']:.3e} ledger={r['ledger']:+.1e}")
