# Lila damping from an action: the Orb's wobble coupled to an explicit bath (the Alt Field).
# Caldeira-Leggett form (Caldeira and Leggett, 1983; Ford, Kac and Mazur, 1965):
#   H = p^2/2 + V(y) + sum_k [ P_k^2/2 + (w_k^2/2) (q_k - c_k y / w_k^2)^2 ]
# Nothing here puts friction in by hand. Energy is conserved exactly; the bath is the ledger.
# Ohmic spectral density J(w) = (pi/2) sum c_k^2/w_k delta(w - w_k) = eta*w gives friction eta,
# so zeta = eta/(2*w0). Choosing eta = 2*eps_L*w0 makes zeta = eps_L (the Lila mapping).
import numpy as np, json, sys
from scipy.integrate import solve_ivp
eps = 43444445555666779/10**19          # eps_L, exact ledger value (3.12)
from scipy.optimize import brentq
A2 = 1.0
def lip(a4):                            # the cup's other stationary point (the lip), for bottom at y = 1
    t1 = 2*A2 + 4*a4; return (-4*a4 + np.sqrt(16*a4**2 - 16*a4*t1))/(2*t1)
MAPPING = sys.argv[2] if len(sys.argv) > 2 else "exact"
if MAPPING == "exact":                  # canonical (Lila audit, Nov 9): lip exactly at 1 - eps_L
    A4 = brentq(lambda a: lip(a) - (1 - eps), -0.1662, -0.1655, xtol=1e-16)
else:                                   # first-order mapping used in the Nov 5 runs
    A4 = -(A2/6)*(1 - eps)
T1 = 2*A2 + 4*A4                        # Lila-mapped cup (3.15), m = 0
V  = lambda y: A2*y**-2 + A4*y**-4 - T1/y
dV = lambda y: -2*A2*y**-3 - 4*A4*y**-5 + T1*y**-2
w0 = np.sqrt(6*A2 + 20*A4 - 2*T1)       # cup curvature at y = 1: exactly 2*eps (first order), 2*eps*1.00145 (exact)
E0 = 1e-9                               # small wobble, as in the adiabatic runs

def run(N, wmax_over_w0, eta, tmax, single=None):
    if single is None:
        dw = wmax_over_w0*w0/N
        wk = (np.arange(N)+0.5)*dw
        ck = np.sqrt((2/np.pi)*eta*wk**2*dw)    # m_k = 1
    else:
        wk = np.array([w0]); ck = np.array([single])
    def f(t, s):
        y, p = s[0], s[1]; q = s[2:2+len(wk)]; P = s[2+len(wk):]
        force = -dV(y) + np.sum(ck*(q - ck*y/wk**2))
        return np.concatenate(([p, force], P, -wk**2*q + ck*y))
    q0 = ck*1.0/wk**2                    # bath starts relaxed around y = 1 (no initial kick)
    s0 = np.concatenate(([1.0, np.sqrt(2*E0)], q0, np.zeros(len(wk))))
    ts = np.linspace(0, tmax, 4000)
    sol = solve_ivp(f, [0, tmax], s0, t_eval=ts, method="DOP853", rtol=1e-10, atol=1e-14)
    y, p = sol.y[0], sol.y[1]; q = sol.y[2:2+len(wk)]; P = sol.y[2+len(wk):]
    Ew   = 0.5*p**2 + V(y) - V(1.0)                                   # wobble (active) energy
    Ealt = np.sum(0.5*P**2 + 0.5*wk[:,None]**2*(q - ck[:,None]*y/wk[:,None]**2)**2, axis=0)
    return sol.t, Ew, Ealt

out = {"mapping": MAPPING, "A4": A4, "lip": lip(A4), "eps_L": eps, "w0": w0, "w0_sq_over_2eps": w0**2/(2*eps)}
# 1) Continuum bath tuned to zeta = eps_L
eta = 2*eps*w0
T = 2*np.pi/w0
NB = int(sys.argv[3]) if len(sys.argv) > 3 else 1200
t, Ew, Ealt = run(N=NB, wmax_over_w0=3.0, eta=eta, tmax=120*T)
m = (t > 5*T) & (t < 110*T)
# fit log of cycle-averaged energy
nb = int((t[m][-1]-t[m][0])//T)
tb = []; eb = []
for i in range(nb):
    sel = m & (t >= t[m][0]+i*T) & (t < t[m][0]+(i+1)*T)
    tb.append(t[sel].mean()); eb.append((Ew[sel]).mean())
slope = np.polyfit(tb, np.log(eb), 1)[0]
zeta_eff = -slope/(2*w0)
ledger = np.max(np.abs((Ew + Ealt) - (Ew[0] + Ealt[0])))/E0
out["continuum"] = dict(N=NB, zeta_eff=zeta_eff, zeta_over_eps=zeta_eff/eps,
                        cycles_per_amplitude_efold=1/(2*np.pi*zeta_eff),
                        frac_wobble_left_at_110T=float(np.interp(110*T, t, Ew))/E0,
                        ledger_rel=ledger, recurrence_cycles=NB/3.0)
# 2) One partner, same coupling strength scale: energy goes out and comes back (no damping)
c1 = np.sqrt(2*eta*w0**2*w0/np.pi*0.05)  # weak single coupling
t1, Ew1, Ea1 = run(N=1, wmax_over_w0=1, eta=0, tmax=400*T, single=c1)
out["single_partner"] = dict(min_wobble_frac=float(Ew1.min()/E0), max_wobble_frac_after_min=float(Ew1[np.argmin(Ew1):].max()/E0),
                             ledger_rel=float(np.max(np.abs(Ew1+Ea1-Ew1[0]-Ea1[0]))/E0))
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "results.json", "w"), indent=1)
print(json.dumps(out, indent=1))
