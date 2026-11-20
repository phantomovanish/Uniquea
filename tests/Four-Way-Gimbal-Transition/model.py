#!/usr/bin/env python3
"""Synthetic bilateral rocker; A and B differ ONLY by the rocker restoring law.

Run: python3 model.py
No network, external data, fitted parameters, or repository changes are involved.
All quantities are dimensionless model units.  See README.md for the equations.
"""
from __future__ import annotations

import csv
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'


@dataclass(frozen=True)
class Parameters:
    # q = [rocker theta, orb y, left fulcrum, right fulcrum, 3 upper scopes, trim]
    masses: tuple = (2., 1., .6, .6, .30, .30, .30, .30)
    damping: tuple = (.18, .16, .18, .18, .12, .12, .12, .12)
    rocker_k: float = 1.0
    nonlinear_gain: float = 1.5
    transition_lo: float = .10
    transition_hi: float = .24
    load_torque: float = .28
    fulcrum_k: tuple = (.8, .8)
    fulcrum_lever: tuple = (-.7, .7)
    orb_coupling_k: tuple = (.8, 1.15)
    fulcrum_ground_k: float = .22
    orb_ground_k: float = .30
    scope_k: tuple = (1.3, 1.3, 1.3, .5)
    scope_lever: tuple = (1., 1., 1., 1.)
    kp: tuple = (4., 4., 4., 16.)
    kd: tuple = (.65, .65, .65, .8)
    target_bias: tuple = (.25, .25, .25, 0.)
    target_amplitude: tuple = (.35, .35, .35, .035)
    target_phase: tuple = (0., math.pi/2, math.pi, 3*math.pi/2)
    omega: float = .75
    cycles: int = 20
    samples_per_cycle: int = 240
    rtol: float = 1e-9
    atol: float = 1e-11
    max_step: float = .08


P = Parameters()
N = 8
# Augmented state [q(8),v(8),W_pos(4),W_neg(4),Q_passive(1)]
POS = slice(16, 20)
NEG = slice(20, 24)
LOSS = 24


def transition(theta, p=P):
    """Smoothstep on |theta|: zero below lo, one above hi, C2 at ends."""
    h = np.clip((abs(theta)-p.transition_lo)/(p.transition_hi-p.transition_lo), 0., 1.)
    return h**3 * (10. + h*(-15. + 6.*h))


def rocker_force(theta, option, p=P):
    """Restoring torque; option B is hardening and is identical to A below lo."""
    if option not in ('A', 'B'):
        raise ValueError('option must be A or B')
    return -p.rocker_k*theta*(1. + (p.nonlinear_gain*transition(theta, p) if option == 'B' else 0.))


def rocker_energy(theta, option, p=P):
    base = .5*p.rocker_k*theta**2
    if option == 'A':
        return base
    if option != 'B':
        raise ValueError('option must be A or B')
    r = abs(theta)
    lo, hi = p.transition_lo, p.transition_hi
    w = hi-lo
    h = np.clip((r-lo)/w, 0., 1.)
    # Exact integral from lo to min(r,hi) of x*S((x-lo)/w) dx.
    integ = w*(lo*(2.5*h**4-3.*h**5+h**6)
               + w*(2.*h**5-2.5*h**6+(6./7.)*h**7))
    integ += .5*max(r*r-hi*hi, 0.)
    return base+p.rocker_k*p.nonlinear_gain*integ


def targets(t, p=P):
    return np.asarray(p.target_bias)+np.asarray(p.target_amplitude)*np.sin(p.omega*t+np.asarray(p.target_phase))


def conservative_force(q, option, p=P):
    """Gradient forces, including constant load torque and reciprocal couplings."""
    th, y, ql, qr = q[:4]
    fulcra = q[2:4]
    b, kf, ko = np.asarray(p.fulcrum_lever), np.asarray(p.fulcrum_k), np.asarray(p.orb_coupling_k)
    a, ks = np.asarray(p.scope_lever), np.asarray(p.scope_k)
    df = fulcra-b*th
    do = y-fulcra
    ds = q[4:8]-a*th
    f = np.zeros(N)
    f[0] = rocker_force(th, option, p)-p.load_torque+np.dot(kf*b, df)+np.dot(ks*a, ds)
    f[1] = -p.orb_ground_k*y-np.dot(ko, do)
    f[2:4] = -kf*df+ko*do-p.fulcrum_ground_k*fulcra
    f[4:8] = -ks*ds
    return f


def energy(q, v, option, p=P):
    th, y = q[:2]
    f = q[2:4]
    ds = q[4:8]-np.asarray(p.scope_lever)*th
    df = f-np.asarray(p.fulcrum_lever)*th
    do = y-f
    return float(.5*np.dot(p.masses, v*v)+rocker_energy(th, option, p)+p.load_torque*th
                 +.5*p.orb_ground_k*y*y+.5*p.fulcrum_ground_k*np.dot(f, f)
                 +.5*np.dot(p.fulcrum_k, df*df)+.5*np.dot(p.orb_coupling_k, do*do)
                 +.5*np.dot(p.scope_k, ds*ds))


def controller(q, v, t, p=P, enabled=True):
    if not enabled:
        return np.zeros(4)
    return np.asarray(p.kp)*(targets(t, p)-q[4:8])-np.asarray(p.kd)*v[4:8]


def rhs(t, z, option, p=P, enabled=True):
    q, v = z[:N], z[N:2*N]
    forces = controller(q, v, t, p, enabled)
    f = conservative_force(q, option, p)-np.asarray(p.damping)*v
    f[4:8] += forces
    power = forces*v[4:8]
    dz = np.zeros_like(z)
    dz[:N] = v
    dz[N:2*N] = f/np.asarray(p.masses)
    dz[POS] = np.maximum(power, 0.)
    dz[NEG] = np.maximum(-power, 0.)
    dz[LOSS] = np.dot(p.damping, v*v)
    return dz


def initial_state(perturb=False):
    z = np.zeros(25)
    z[:8] = [.06, .035, -.015, .025, .20, .20, .20, 0.]
    if perturb:
        z[0] += .04
        z[1] += .04
    return z


def run(option, p=P, perturb=False, enabled=True, *, rtol=None, atol=None, max_step=None):
    T = 2*math.pi/p.omega
    t = np.linspace(0., p.cycles*T, p.cycles*p.samples_per_cycle+1)
    sol = solve_ivp(lambda t,z: rhs(t,z,option,p,enabled), (0.,t[-1]), initial_state(perturb),
                    method='DOP853', t_eval=t, rtol=p.rtol if rtol is None else rtol,
                    atol=p.atol if atol is None else atol,
                    max_step=p.max_step if max_step is None else max_step)
    if not sol.success or len(sol.t) != len(t) or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f'{option} solver failed or nonfinite: {sol.message}')
    return sol.t, sol.y.T


def observed_forces(q, option, p=P):
    theta = q[0]
    scope = np.asarray(p.scope_k)*np.asarray(p.scope_lever)*(q[4:8]-np.asarray(p.scope_lever)*theta)
    fulcrum = np.asarray(p.fulcrum_k)*np.asarray(p.fulcrum_lever)*(q[2:4]-np.asarray(p.fulcrum_lever)*theta)
    return float(np.sum(scope[:3])), float(scope[3]), float(np.sum(fulcrum)), float(rocker_force(theta, option, p))


def wrap_deg(x):
    return (x+180.)%360.-180.


def harmonic_phase(t, value, p=P):
    # Best-fit constant + sin + cos: phase is relative to sin(omega*t).
    X = np.column_stack((np.ones(len(t)), np.sin(p.omega*t), np.cos(p.omega*t)))
    c = np.linalg.lstsq(X, value, rcond=None)[0]
    return float(math.hypot(c[1],c[2])), float(wrap_deg(math.degrees(math.atan2(c[2],c[1]))))


def ledger(z, option, p=P):
    H = np.array([energy(row[:8], row[8:16], option, p) for row in z])
    residual = H-H[0]-np.sum(z[:, POS], axis=1)+np.sum(z[:, NEG], axis=1)+z[:, LOSS]
    return H, residual


def stage_rows(t, z, option, p=P):
    n = p.samples_per_cycle//4
    H, _ = ledger(z, option, p)
    rows = []
    for j in range(p.cycles*4):
        i, k = j*n, (j+1)*n
        zz, tt = z[i:k+1], t[i:k+1]
        a = abs(zz[:,0])
        # Same endpoints and quadrature resolution in both options; occupancy is sampled.
        upper, trim, fulcrum, rock = np.asarray([observed_forces(q[:8], option, p) for q in zz]).T
        mean = lambda x: float(np.trapezoid(x, tt)/(tt[-1]-tt[0]))
        dp = z[k,POS]-z[i,POS]; dn = z[k,NEG]-z[i,NEG]; dq = z[k,LOSS]-z[i,LOSS]
        dH = H[k]-H[i]
        rows.append(dict(option=option,cycle=j//4+1,stage=j%4+1,
                         phase_start_deg=(j%4)*90,phase_end_deg=(j%4+1)*90,
                         time_start=float(t[i]),time_end=float(t[k]),
                         theta_start=float(z[i,0]),theta_end=float(z[k,0]),
                         theta_min=float(np.min(zz[:,0])),theta_max=float(np.max(zz[:,0])),
                         orb_min=float(np.min(zz[:,1])),orb_max=float(np.max(zz[:,1])),
                         mean_upper_torque=mean(upper),mean_trim_torque=mean(trim),
                         mean_fulcra_torque=mean(fulcrum),mean_rocker_restoring_torque=mean(rock),
                         load_torque=-p.load_torque,delta_H=float(dH),
                         upper_positive_work=float(np.sum(dp[:3])),trim_positive_work=float(dp[3]),
                         upper_negative_work_magnitude=float(np.sum(dn[:3])),trim_negative_work_magnitude=float(dn[3]),
                         passive_loss=float(dq),
                         ledger_residual=float(dH-np.sum(dp)+np.sum(dn)+dq),
                         fraction_below_transition=float(np.mean(a<=p.transition_lo)),
                         fraction_within_transition=float(np.mean((a>p.transition_lo)&(a<p.transition_hi))),
                         fraction_above_transition=float(np.mean(a>=p.transition_hi))))
    return rows


def crossings(t, value):
    out = []
    for i in np.flatnonzero(value[:-1]*value[1:] < 0.):
        f = -value[i]/(value[i+1]-value[i])
        out.append(float(t[i]+f*(t[i+1]-t[i])))
    return out


def metrics(t, z, option, perturbed, tight, p=P):
    sp = p.samples_per_cycle
    H, resid = ledger(z, option, p)
    last = z[-sp-1:]; lt = t[-sp-1:]
    fit = [harmonic_phase(lt, last[:,i]) for i in range(4,8)]
    phases = [v[1] for v in fit]
    # Initial-value response measured in dimensionless physical-state Euclidean norm.
    diff = np.linalg.norm(perturbed[:,:16]-z[:,:16], axis=1)
    turn_t = crossings(lt, last[:,8])
    zero_t = crossings(lt, last[:,0])
    # max state discrepancy includes both physical trajectories and integrated ledger states
    return dict(option=option, H_initial=float(H[0]),H_final=float(H[-1]),
                delta_H=float(H[-1]-H[0]),positive_work_by_control=z[-1,POS].tolist(),
                negative_work_magnitude_by_control=z[-1,NEG].tolist(),
                net_actuator_work=float(np.sum(z[-1,POS]-z[-1,NEG])),
                passive_loss=float(z[-1,LOSS]),max_abs_ledger_residual=float(np.max(abs(resid))),
                final_ledger_residual=float(resid[-1]),
                last_cycle_theta_min=float(np.min(last[:,0])),last_cycle_theta_max=float(np.max(last[:,0])),
                last_cycle_theta_min_deg=float(np.degrees(np.min(last[:,0]))),last_cycle_theta_max_deg=float(np.degrees(np.max(last[:,0]))),
                last_cycle_orb_min=float(np.min(last[:,1])),last_cycle_orb_max=float(np.max(last[:,1])),
                last_cycle_fulcra_ranges=[[float(np.min(last[:,i])),float(np.max(last[:,i]))] for i in (2,3)],
                last_cycle_scope_half_ranges=[float(np.ptp(last[:,i])/2.) for i in range(4,8)],
                last_cycle_scope_fundamental_amplitudes=[x[0] for x in fit],
                target_phase_deg=[0.,90.,180.,-90.],actual_fundamental_phase_deg=phases,
                actual_side_phases_relative_scope1_deg=[float(wrap_deg(phases[1]-phases[0])),float(wrap_deg(phases[3]-phases[0]))],
                actual_scope2_minus_scope4_phase_deg=float(wrap_deg(phases[1]-phases[3])),
                last_cycle_turning_times=turn_t,last_cycle_zero_crossing_times=zero_t,
                initial_perturbation_norm=float(diff[0]),difference_at_end_cycle_1=float(diff[sp]),
                difference_at_end_cycle_5=float(diff[5*sp]),difference_at_end_cycle_20=float(diff[-1]),
                difference_max_all_samples=float(np.max(diff)),
                last_cycle_perturbation_max=float(np.max(diff[-sp-1:])),
                baseline_last_cycle_closure_norm=float(np.linalg.norm(z[-1,:16]-z[-sp-1,:16])),
                tight_solver_max_physical_state_difference=float(np.max(abs(z[:,:16]-tight[:,:16]))),
                tight_solver_max_ledger_state_difference=float(np.max(abs(z[:,16:]-tight[:,16:]))))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def plot(t, runs, p=P):
    sp = p.samples_per_cycle
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    times = (t[-sp-1:]-t[-sp-1])/(2*math.pi/p.omega)*360
    for option,color in [('A','#1d5f9e'),('B','#d46b24')]:
        z = runs[option][-sp-1:]
        ax[0].plot(times, np.degrees(z[:,0]), color=color,label=f'{option}: rocker')
        ax[1].plot(times, z[:,1],color=color,label=f'{option}: orb')
        H,R = ledger(runs[option],option,p)
        ax[2].plot(times,R[-sp-1:],color=color,label=f'{option}: residual')
    ax[0].axhline(0,color='0.55',linewidth=.7)
    for a in ax:
        for ph in (0,90,180,270,360): a.axvline(ph,color='0.85',linewidth=.7,linestyle='--')
        a.grid(alpha=.16); a.legend(loc='best',fontsize=9)
    ax[0].set_ylabel('Rocker angle (deg)'); ax[1].set_ylabel('Orb displacement'); ax[2].set_ylabel('Energy residual')
    ax[2].set_xlabel('Drive phase in final cycle (degrees)')
    ax[0].set_title('Synthetic driven model: final cycle; phase is a timing marker, not a physical flip')
    label = 'Manus (version not recorded) Uniquea Field Observations UFO'
    fig.get_layout_engine().set(rect=(0, .04, 1, .96))
    fig.text(.5, .015, label, ha='center', va='bottom', fontsize=8)
    fig.savefig(OUT/'final_cycle.png', dpi=150, metadata={'Description': label})
    plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    runs = {}; summary = {}; all_stages = []
    for option in ('A','B'):
        print(f'Running {option}: baseline, perturbation, tighter solver...',flush=True)
        t,z = run(option); _,p = run(option,perturb=True)
        _,tight = run(option,rtol=2e-11,atol=2e-13,max_step=.04)
        runs[option] = z
        summary[option] = metrics(t,z,option,p,tight)
        all_stages.extend(stage_rows(t,z,option))
        H,res = ledger(z,option)
        fields = ['time','theta','orb','fulcrum_L','fulcrum_R','upper_1','upper_2','upper_3','lower_trim']
        fields += ['v_'+x for x in fields[1:]]
        fields += ['W_positive_'+str(i) for i in range(1,5)]+['W_negative_magnitude_'+str(i) for i in range(1,5)]+['Q_passive','H','ledger_residual']
        data = np.column_stack((t,z,H,res))
        with (OUT/f'trajectory_{option}.csv').open('w',newline='') as f:
            w = csv.writer(f); w.writerow(fields); w.writerows(data)
        diff = np.linalg.norm(p[:,:16]-z[:,:16],axis=1)
        with (OUT/f'perturbation_{option}.csv').open('w',newline='') as f:
            w = csv.writer(f); w.writerow(['time','physical_state_difference_norm']); w.writerows(zip(t,diff))
    print('Running no-controller passive check...',flush=True)
    passive = {}
    for option in ('A','B'):
        tp,zp=run(option,enabled=False); Hp,rp=ledger(zp,option)
        passive[option] = dict(H_initial=float(Hp[0]),H_final=float(Hp[-1]),
                               max_energy_increase_per_sample=float(np.max(np.diff(Hp))),
                               passive_loss=float(zp[-1,LOSS]),max_abs_ledger_residual=float(np.max(abs(rp))))
    write_csv(OUT/'stages.csv',all_stages)
    output = dict(status='synthetic, driven, unvalidated',parameters=asdict(P),
                  initial_state=initial_state().tolist(),perturbation_theta_and_orb=[.04,.04],
                  software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
                  solver='DOP853', comparison=summary,passive_control=passive)
    (OUT/'summary.json').write_text(json.dumps(output,indent=2)+'\n')
    plot(t,runs)
    print('Finished. Results:',OUT,flush=True)
    for o in ('A','B'):
        m=summary[o]
        print(f"{o}: theta last [{m['last_cycle_theta_min']:.6g}, {m['last_cycle_theta_max']:.6g}], "
              f"net work {m['net_actuator_work']:.6g}, passive loss {m['passive_loss']:.6g}, "
              f"ledger max {m['max_abs_ledger_residual']:.3g}, final perturb norm {m['difference_at_end_cycle_20']:.3g}",flush=True)


if __name__ == '__main__':
    main()
