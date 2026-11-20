#!/usr/bin/env python3
"""Prospective, driven synthetic hybrid of the two archived rocker laws.

Offline, standalone: the neighboring reference_model.py is an unchanged copy.
See PROTOCOL.md and README.md. No observational claim is made.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import reference_model as ref

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'
P = ref.P
WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)
N = ref.N
POS, NEG, LOSS = ref.POS, ref.NEG, ref.LOSS
TIGHT = dict(rtol=2e-11, atol=2e-13, max_step=.04)
VARIANTS = ('baseline_standard', 'baseline_tight', 'plus_standard', 'plus_tight', 'minus_standard')
COLUMNS = ['theta','orb','fulcrum_L','fulcrum_R','upper_1','upper_2','upper_3','lower_trim']
COLUMNS += ['v_'+x for x in COLUMNS]
COLUMNS += ['W_positive_'+str(i) for i in range(1,5)] + ['W_negative_magnitude_'+str(i) for i in range(1,5)] + ['Q_passive']


def key(w):
    return f'lambda_{int(round(w*100)):03d}'


def check_weight(w):
    if not np.isfinite(w) or w < 0. or w > 1.:
        raise ValueError('lambda must be finite and in [0,1]')


def rocker_force(theta, w, p=P):
    check_weight(w)
    if w == 0.: return ref.rocker_force(theta, 'A', p)
    if w == 1.: return ref.rocker_force(theta, 'B', p)
    return (1.-w)*ref.rocker_force(theta, 'A', p) + w*ref.rocker_force(theta, 'B', p)


def rocker_energy(theta, w, p=P):
    check_weight(w)
    if w == 0.: return ref.rocker_energy(theta, 'A', p)
    if w == 1.: return ref.rocker_energy(theta, 'B', p)
    return (1.-w)*ref.rocker_energy(theta, 'A', p) + w*ref.rocker_energy(theta, 'B', p)


def conservative_force(q, w, p=P):
    check_weight(w)
    if w == 0.: return ref.conservative_force(q, 'A', p)
    if w == 1.: return ref.conservative_force(q, 'B', p)
    f = ref.conservative_force(q, 'A', p)
    f[0] += rocker_force(q[0], w, p) - ref.rocker_force(q[0], 'A', p)
    return f


def energy(q, v, w, p=P):
    check_weight(w)
    if w == 0.: return ref.energy(q, v, 'A', p)
    if w == 1.: return ref.energy(q, v, 'B', p)
    return ref.energy(q, v, 'A', p) + rocker_energy(q[0], w, p) - ref.rocker_energy(q[0], 'A', p)


def rhs(t, z, w, p=P, enabled=True):
    check_weight(w)
    if w == 0.: return ref.rhs(t, z, 'A', p, enabled)
    if w == 1.: return ref.rhs(t, z, 'B', p, enabled)
    q, v = z[:N], z[N:2*N]
    controls = ref.controller(q, v, t, p, enabled)
    f = conservative_force(q, w, p) - np.asarray(p.damping)*v
    f[4:8] += controls
    power = controls*v[4:8]
    dz = np.zeros_like(z)
    dz[:N] = v
    dz[N:2*N] = f/np.asarray(p.masses)
    dz[POS] = np.maximum(power, 0.)
    dz[NEG] = np.maximum(-power, 0.)
    dz[LOSS] = np.dot(p.damping, v*v)
    return dz


def initial_state(perturbation=0):
    if perturbation not in (-1, 0, 1): raise ValueError('perturbation must be -1, 0, or 1')
    z = ref.initial_state()
    z[0] += .04*perturbation
    z[1] += .04*perturbation
    return z


def run(w, p=P, perturbation=0, tight=False, enabled=True):
    check_weight(w)
    period = 2*math.pi/p.omega
    t = np.linspace(0., p.cycles*period, p.cycles*p.samples_per_cycle+1)
    settings = TIGHT if tight else dict(rtol=p.rtol, atol=p.atol, max_step=p.max_step)
    sol = solve_ivp(lambda tt,zz: rhs(tt,zz,w,p,enabled), (0.,t[-1]), initial_state(perturbation),
                    method='DOP853', t_eval=t, **settings)
    if not sol.success or len(sol.t) != len(t) or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f'lambda={w} integration failed: {sol.message}')
    return sol.t, sol.y.T


def ledger(z, w, p=P):
    H = np.array([energy(row[:N], row[N:2*N], w, p) for row in z])
    residual = H-H[0]-np.sum(z[:,POS],axis=1)+np.sum(z[:,NEG],axis=1)+z[:,LOSS]
    return H, residual


def stage_rows(t, z, w, p=P):
    n = p.samples_per_cycle//4
    H, _ = ledger(z,w,p)
    result = []
    for j in range(p.cycles*4):
        i, k = j*n, (j+1)*n
        zz, tt = z[i:k+1], t[i:k+1]
        q = zz[:,:8]
        sf = np.asarray(p.scope_k)*np.asarray(p.scope_lever)*(q[:,4:8]-np.asarray(p.scope_lever)*q[:,0,None])
        ff = np.asarray(p.fulcrum_k)*np.asarray(p.fulcrum_lever)*(q[:,2:4]-np.asarray(p.fulcrum_lever)*q[:,0,None])
        rf = np.array([rocker_force(theta,w,p) for theta in q[:,0]])
        mean = lambda v: float(np.trapezoid(v,tt)/(tt[-1]-tt[0]))
        dp = z[k,POS]-z[i,POS]; dn = z[k,NEG]-z[i,NEG]; dq = z[k,LOSS]-z[i,LOSS]
        dH = H[k]-H[i]; a = np.abs(q[:,0])
        result.append(dict(lambda_weight=w, cycle=j//4+1, stage=j%4+1,
            phase_start_deg=(j%4)*90, phase_end_deg=(j%4+1)*90,
            time_start=float(t[i]),time_end=float(t[k]),theta_start=float(q[0,0]),theta_end=float(q[-1,0]),
            theta_min=float(np.min(q[:,0])),theta_max=float(np.max(q[:,0])),orb_min=float(np.min(q[:,1])),orb_max=float(np.max(q[:,1])),
            mean_upper_torque=mean(np.sum(sf[:,:3],axis=1)),mean_trim_torque=mean(sf[:,3]),
            mean_fulcra_torque=mean(np.sum(ff,axis=1)),mean_rocker_torque=mean(rf),load_torque=-p.load_torque,
            delta_H=float(dH),upper_positive_work=float(np.sum(dp[:3])),trim_positive_work=float(dp[3]),
            upper_negative_work_magnitude=float(np.sum(dn[:3])),trim_negative_work_magnitude=float(dn[3]),
            passive_loss=float(dq),ledger_residual=float(dH-np.sum(dp)+np.sum(dn)+dq),
            fraction_below_or_at_lo=float(np.mean(a<=p.transition_lo)),
            fraction_inside_transition=float(np.mean((a>p.transition_lo)&(a<p.transition_hi))),
            fraction_at_or_above_hi=float(np.mean(a>=p.transition_hi))))
    return result


def pair_metrics(base, pert, sp):
    d = np.linalg.norm(pert[:,:16]-base[:,:16],axis=1)
    return dict(initial=float(d[0]),end_cycle_1=float(d[sp]),end_cycle_5=float(d[5*sp]),
                end_cycle_20=float(d[-1]),maximum=float(np.max(d)),last_cycle_max=float(np.max(d[-sp-1:])))


def metrics(t, runs, w, p=P):
    z = runs['baseline_standard']; tight=runs['baseline_tight']; sp=p.samples_per_cycle
    H,R = ledger(z,w,p)
    lt=t[-sp-1:]; last=z[-sp-1:]; samples=last[:-1,0]; a=np.abs(samples)
    plus=pair_metrics(z,runs['plus_standard'],sp); minus=pair_metrics(z,runs['minus_standard'],sp)
    compare = lambda a,b: dict(max_physical_state_component=float(np.max(np.abs(a[:,:16]-b[:,:16]))),
                                max_ledger_state_component=float(np.max(np.abs(a[:,16:]-b[:,16:]))))
    crossings=ref.crossings(lt,last[:,0]); turns=ref.crossings(lt,last[:,8])
    fit=[ref.harmonic_phase(lt,last[:,i]) for i in range(4,8)]
    return dict(lambda_weight=w,H_initial=float(H[0]),H_final=float(H[-1]),delta_H=float(H[-1]-H[0]),
        positive_work_by_control=z[-1,POS].tolist(),negative_work_magnitude_by_control=z[-1,NEG].tolist(),
        total_positive_work=float(np.sum(z[-1,POS])),total_negative_work_magnitude=float(np.sum(z[-1,NEG])),
        net_actuator_work=float(np.sum(z[-1,POS]-z[-1,NEG])),passive_loss=float(z[-1,LOSS]),
        max_abs_ledger_residual=float(np.max(np.abs(R))),final_ledger_residual=float(R[-1]),
        last_cycle_theta_min=float(np.min(last[:,0])),last_cycle_theta_max=float(np.max(last[:,0])),
        last_cycle_theta_rms=float(np.sqrt(np.trapezoid(last[:,0]**2,lt)/(lt[-1]-lt[0]))),
        last_cycle_theta_mean=float(np.trapezoid(last[:,0],lt)/(lt[-1]-lt[0])),
        last_cycle_theta_half_range=float(np.ptp(last[:,0])/2.),last_cycle_theta_max_abs=float(np.max(np.abs(last[:,0]))),
        last_cycle_orb_min=float(np.min(last[:,1])),last_cycle_orb_max=float(np.max(last[:,1])),
        last_cycle_turning_times=turns,last_cycle_zero_crossing_times=crossings,
        last_cycle_below_or_at_lo_fraction=float(np.mean(a<=p.transition_lo)),
        last_cycle_transition_fraction=float(np.mean((a>p.transition_lo)&(a<p.transition_hi))),
        last_cycle_at_or_above_hi_fraction=float(np.mean(a>=p.transition_hi)),
        last_cycle_scope_fundamental_amplitudes=[f[0] for f in fit],
        last_cycle_scope_fundamental_phases_deg=[f[1] for f in fit],
        baseline_last_cycle_closure_norm=float(np.linalg.norm(z[-1,:16]-z[-sp-1,:16])),
        perturb_plus=plus,perturb_minus=minus,
        normal_vs_tight_baseline=compare(z,tight),
        normal_vs_tight_plus=compare(runs['plus_standard'],runs['plus_tight']))


def compute_one(w):
    runs={}; t=None
    for name,pert,tight in [('baseline_standard',0,False),('baseline_tight',0,True),
                            ('plus_standard',1,False),('plus_tight',1,True),('minus_standard',-1,False)]:
        print(f'{key(w)} {name}: integrating',flush=True)
        t,z=run(w,perturbation=pert,tight=tight); runs[name]=z
    if w in (0.,1.):
        option = 'A' if w==0. else 'B'
        print(f'{key(w)} original {option}: independent endpoint integration',flush=True)
        tr,zref=ref.run(option)
        if not np.array_equal(t,tr) or not np.array_equal(runs['baseline_standard'],zref):
            raise AssertionError(f'full reference endpoint {option} differs')
        runs['reference_endpoint']=zref
    return w,t,runs


def write_csv(path, rows):
    if not rows: raise ValueError('empty CSV')
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def plot_results(t, collected, summ):
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),constrained_layout=True)
    colors=plt.cm.viridis(np.linspace(.05,.95,len(WEIGHTS)))
    sp=P.samples_per_cycle; phase=(t[-sp-1:]-t[-sp-1])/(2*math.pi/P.omega)*360
    for w,c in zip(WEIGHTS,colors):
        zz=collected[w]['baseline_standard'][-sp-1:]
        ax[0,0].plot(phase,zz[:,0],color=c,label=f'λ={w:g}')
        ax[0,1].plot(phase,zz[:,1],color=c,label=f'λ={w:g}')
    for lim in (-P.transition_lo,P.transition_lo): ax[0,0].axhline(lim,color='0.55',lw=.8,ls='--')
    rms=[summ[key(w)]['last_cycle_theta_rms'] for w in WEIGHTS]
    half=[summ[key(w)]['last_cycle_theta_half_range'] for w in WEIGHTS]
    ax[1,0].plot(WEIGHTS,rms,'o-',label='RMS about zero')
    ax[1,0].plot(WEIGHTS,half,'s-',label='Half-range')
    ax[1,1].plot(WEIGHTS,[summ[key(w)]['total_positive_work'] for w in WEIGHTS],'o-',label='Positive work')
    ax[1,1].plot(WEIGHTS,[summ[key(w)]['total_negative_work_magnitude'] for w in WEIGHTS],'s-',label='Negative work magnitude')
    ax[0,0].set(title='Final-cycle rocker (dashed: ±0.10 rad)',ylabel='θ (rad)',xlabel='Drive phase (degrees)')
    ax[0,1].set(title='Final-cycle orb displacement',ylabel='y (model units)',xlabel='Drive phase (degrees)')
    ax[1,0].set(title='Angle metrics are distinct',ylabel='Radians',xlabel='Blend weight λ')
    ax[1,1].set(title='20-cycle mechanical actuator work',ylabel='Model energy units',xlabel='Blend weight λ')
    for a in ax.flat: a.grid(alpha=.2); a.legend(fontsize=8)
    fig.suptitle('Synthetic, externally driven model — not physical validation',fontsize=12)
    fig.savefig(OUT/'overview.png',dpi=160); plt.close(fig)


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def verify_outputs(expected, rows_count):
    path=OUT/'trajectories.npz'
    monotonicity={}
    with np.load(path,allow_pickle=False) as archive:
        if set(archive.files)!=set(expected): raise AssertionError('NPZ keys differ')
        for name, original in expected.items():
            data=archive[name]
            if data.shape!=original.shape or not np.array_equal(data,original) or not np.all(np.isfinite(data)):
                raise AssertionError(f'archive round-trip mismatch: {name}')
            if name!='time':
                if data.shape!=(4801,25): raise AssertionError(f'wrong shape {name}: {data.shape}')
                work_diff=np.diff(data[:,16:24],axis=0)
                passive_diff=np.diff(data[:,24],axis=0)
                i,j=np.unravel_index(np.argmin(work_diff),work_diff.shape)
                monotonicity[name]=dict(
                    minimum_work_increment=float(work_diff[i,j]),
                    minimum_work_increment_time_index=int(i),
                    minimum_work_increment_state_column=int(j+16),
                    work_decreases_exceeding_1e_9=int(np.sum(work_diff < -1e-9)),
                    minimum_passive_increment=float(np.min(passive_diff)),
                    passive_decreases_exceeding_1e_9=int(np.sum(passive_diff < -1e-9)))
    for name, count in rows_count.items():
        with (OUT/name).open(newline='',encoding='utf-8') as f:
            rows=list(csv.reader(f))
        if len(rows)!=count+1 or any(len(r)!=len(rows[0]) for r in rows): raise AssertionError(f'bad CSV {name}')
        if not np.isfinite(np.array(rows[1:],dtype=float)).all(): raise AssertionError(f'nonfinite CSV {name}')
    return dict(roundtrip_arrays=len(expected),n_time_samples=len(expected['time']),csv_rows=rows_count,
                all_arrays_exact_after_roundtrip=True,all_csv_finite_and_rectangular=True,
                work_monotone_with_tolerance_1e_9=all(x['work_decreases_exceeding_1e_9']==0 for x in monotonicity.values()),
                passive_monotone_with_tolerance_1e_9=all(x['passive_decreases_exceeding_1e_9']==0 for x in monotonicity.values()),
                monotonicity_by_run=monotonicity)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    # Files here are entirely generated by this program; clear stale files for an honest inventory.
    for path in OUT.iterdir():
        if path.is_file(): path.unlink()
        else: raise RuntimeError('unexpected output subdirectory; remove it manually')
    collected={}; times={}
    with ProcessPoolExecutor(max_workers=5) as pool:
        futures=[pool.submit(compute_one,w) for w in WEIGHTS]
        for future in as_completed(futures):
            w,t,runs=future.result(); collected[w]=runs; times[w]=t
            print(f'{key(w)} finished',flush=True)
    t=times[0.]
    if any(not np.array_equal(t,times[w]) for w in WEIGHTS): raise AssertionError('time grids differ')
    arrays={'time':t}; summaries={}; stages=[]; rowcounts={}
    for w in WEIGHTS:
        runs=collected[w]; kk=key(w)
        for name,z in runs.items(): arrays[kk+'__'+name]=z
        summaries[kk]=metrics(t,runs,w); stages.extend(stage_rows(t,runs['baseline_standard'],w))
        z=runs['baseline_standard']; H,R=ledger(z,w)
        data=np.column_stack((t,z,H,R)); fn=OUT/f'{kk}_baseline.csv'
        with fn.open('w',newline='',encoding='utf-8') as f:
            writer=csv.writer(f); writer.writerow(['time']+COLUMNS+['H','ledger_residual']); writer.writerows(data)
        rowcounts[fn.name]=len(t)
        for name in VARIANTS:
            _,res=ledger(runs[name],w)
            if np.max(np.abs(res))>=1e-7: raise AssertionError(f'ledger residual too large: {kk} {name}')
        for pair in ('normal_vs_tight_baseline','normal_vs_tight_plus'):
            cmp=summaries[kk][pair]
            if cmp['max_physical_state_component']>=5e-7 or cmp['max_ledger_state_component']>=1e-4:
                raise AssertionError(f'tighter solver comparison failed: {kk} {pair}')
    np.savez_compressed(OUT/'trajectories.npz',**arrays)
    write_csv(OUT/'stages.csv',stages); rowcounts['stages.csv']=len(stages)
    flat=[]
    for w in WEIGHTS:
        d=summaries[key(w)]
        flat.append({name:d[name] for name in ('lambda_weight','last_cycle_theta_min','last_cycle_theta_max','last_cycle_theta_rms',
            'last_cycle_theta_mean','last_cycle_theta_half_range','last_cycle_below_or_at_lo_fraction','last_cycle_transition_fraction',
            'last_cycle_at_or_above_hi_fraction','total_positive_work','total_negative_work_magnitude','net_actuator_work','passive_loss',
            'delta_H','max_abs_ledger_residual')})
    write_csv(OUT/'metrics.csv',flat); rowcounts['metrics.csv']=len(flat)
    plot_results(t,collected,summaries)
    integrity=verify_outputs(arrays,rowcounts)
    output=dict(status='synthetic, externally driven, unvalidated',protocol='PROTOCOL.md',reference_source_sha256=digest(HERE/'reference_model.py'),
        weights=list(WEIGHTS),parameters=asdict(P),initial_state=initial_state().tolist(),
        perturbation_plus=[.04,.04],perturbation_minus=[-.04,-.04],
        solver='DOP853',tighter_solver_settings=TIGHT,
        software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        state_columns=COLUMNS,archive_schema='time plus lambda_NNN__VARIANT arrays, each 4801 by 25',
        comparison=summaries,endpoint_reference_trajectory_bitwise_equal={'A':True,'B':True},integrity=integrity)
    (OUT/'summary.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    paths=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt')
    (OUT/'SHA256SUMS.txt').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in paths),encoding='utf-8')
    print('All runs complete; array/CSV integrity passed.',flush=True)
    if not integrity['work_monotone_with_tolerance_1e_9']:
        worst=min((v['minimum_work_increment'],k) for k,v in integrity['monotonicity_by_run'].items())
        print(f'WARNING: gross cumulative work monotonicity at tolerance 1e-9 FAILED; worst {worst}. '
              'Values are retained unmodified; see summary.json.',flush=True)
    for row in flat:
        print('lambda={lambda_weight:.2f} theta=[{last_cycle_theta_min:.6f},{last_cycle_theta_max:.6f}] rms={last_cycle_theta_rms:.6f} W+={total_positive_work:.6f} W-={total_negative_work_magnitude:.6f} Q={passive_loss:.6f} residual={max_abs_ledger_residual:.3g}'.format(**row),flush=True)


if __name__=='__main__': main()
