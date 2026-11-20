#!/usr/bin/env python3
"""Ten-cycle controller-ablation experiment using model.py in this directory.

No network or repository writes. Execute: python3 ablation.py
Synthetic units; masks remove the entire controller force, not damping alone.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
import numpy as np
import scipy
from scipy.integrate import solve_ivp

import model as m

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'
P = replace(m.P, cycles=10)  # Only the observation horizon differs from model.run.
MASKS = {
    'all_on': (True, True, True, True),
    'lower_only': (False, False, False, True),
    'upper_only': (True, True, True, False),
    'all_off': (False, False, False, False),
}
OPTIONS = ('A', 'B')


def checked_mask(mask):
    """Return a four-entry boolean mask, refusing ambiguous values."""
    a = np.asarray(mask)
    if a.shape != (4,) or a.dtype.kind != 'b':
        raise ValueError('mask must contain exactly four booleans')
    return a


def masked_rhs(t, z, option, p=m.P, mask=MASKS['all_on']):
    """Original equations with only controller forces/power masked.

    Fi = mask_i * (kp_i*(target_i(t)-q_i) - kd_i*v_i).
    Passive damping and all conservative coupling forces are unchanged. Positive
    and negative work rates are split after masking, using the corresponding
    scope's velocity, just as in model.rhs.
    """
    mask = checked_mask(mask)
    q, v = z[:m.N], z[m.N:2*m.N]
    forces = np.where(mask, m.controller(q, v, t, p, enabled=True), 0.)
    f = m.conservative_force(q, option, p) - np.asarray(p.damping)*v
    f[4:8] += forces
    power = forces*v[4:8]
    dz = np.zeros_like(z)
    dz[:m.N] = v
    dz[m.N:2*m.N] = f/np.asarray(p.masses)
    dz[m.POS] = np.maximum(power, 0.)
    dz[m.NEG] = np.maximum(-power, 0.)
    dz[m.LOSS] = np.dot(p.damping, v*v)
    return dz


def run(option, mask, perturb=False, p=P):
    """Match model.run's DOP853 solver and sampling, with a 10-cycle horizon."""
    mask = checked_mask(mask)
    period = 2*math.pi/p.omega
    t = np.linspace(0., p.cycles*period, p.cycles*p.samples_per_cycle+1)
    sol = solve_ivp(lambda time, z: masked_rhs(time, z, option, p, mask),
                    (0., t[-1]), m.initial_state(perturb), method='DOP853',
                    t_eval=t, rtol=p.rtol, atol=p.atol, max_step=p.max_step)
    if not sol.success or len(sol.t) != len(t) or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f'{option}/{mask.tolist()} solver failure: {sol.message}')
    return sol.t, sol.y.T


def rms(values, t):
    return float(np.sqrt(np.trapezoid(values*values, t)/(t[-1]-t[0])))


def trajectory_metrics(t, z, option, p=P):
    """Sampled extrema, continuous-time trapezoidal RMS, and integrated ledger."""
    theta = z[:, 0]
    sl = slice((p.cycles-1)*p.samples_per_cycle, None)
    tf, xf = t[sl], theta[sl]
    mean = float(np.trapezoid(xf, tf)/(tf[-1]-tf[0]))
    H, residual = m.ledger(z, option, p)
    wp, wn = z[-1, m.POS], z[-1, m.NEG]
    return {
        'rocker': {
            'all_10_cycles_sampled_min': float(np.min(theta)),
            'all_10_cycles_sampled_max': float(np.max(theta)),
            'last_cycle_sampled_min': float(np.min(xf)),
            'last_cycle_sampled_max': float(np.max(xf)),
            'last_cycle_mean': mean,
            'last_cycle_rms_about_zero': rms(xf, tf),
            'last_cycle_rms_about_own_mean': rms(xf-mean, tf),
        },
        'energy_ledger': {
            'H_initial': float(H[0]), 'H_final': float(H[-1]),
            'delta_H': float(H[-1]-H[0]),
            'positive_actuator_work_by_control': wp.tolist(),
            'negative_actuator_work_magnitude_by_control': wn.tolist(),
            'positive_actuator_work_total': float(np.sum(wp)),
            'negative_actuator_work_magnitude_total': float(np.sum(wn)),
            'net_actuator_work': float(np.sum(wp-wn)),
            'passive_bath': float(z[-1, m.LOSS]),
            'final_residual': float(residual[-1]),
            'max_abs_residual_all_samples': float(np.max(np.abs(residual))),
        },
    }


def disturbance_metrics(baseline, perturbed, p=P):
    """Euclidean physical-state separation, excluding work/bath coordinates."""
    norm = np.linalg.norm(perturbed[:, :2*m.N]-baseline[:, :2*m.N], axis=1)
    sp = p.samples_per_cycle
    return {
        'initial': float(norm[0]),
        'end_cycle_1': float(norm[sp]),
        'end_cycle_5': float(norm[5*sp]),
        'end_cycle_10': float(norm[10*sp]),
        'max_over_all_samples': float(np.max(norm)),
        'max_sample_index': int(np.argmax(norm)),
    }


def sorted_masks(results, option, kind):
    def value(name):
        case = results[option][name]
        if kind == 'rms_about_zero':
            return case['baseline']['rocker']['last_cycle_rms_about_zero']
        if kind == 'rms_about_own_mean':
            return case['baseline']['rocker']['last_cycle_rms_about_own_mean']
        return case['disturbance'][kind]
    return [{'mask': name, 'value': value(name)} for name in sorted(MASKS, key=value)]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results, arrays = {}, {}
    common_t = None
    for option in OPTIONS:
        results[option] = {}
        for name, mask in MASKS.items():
            print(f'Running {option}/{name}: baseline and perturbation', flush=True)
            t, base = run(option, mask, perturb=False)
            tp, pert = run(option, mask, perturb=True)
            if not np.array_equal(t, tp) or (common_t is not None and not np.array_equal(t, common_t)):
                raise RuntimeError('sampling grids differ')
            common_t = t
            results[option][name] = {
                'enabled_mask_upper1_upper2_upper3_lower': list(mask),
                'baseline': trajectory_metrics(t, base, option),
                'perturbed': trajectory_metrics(t, pert, option),
                'disturbance': disturbance_metrics(base, pert),
            }
            arrays[f'{option}_{name}_baseline'] = base
            arrays[f'{option}_{name}_perturbed'] = pert
    arrays['time'] = common_t
    np.savez_compressed(OUT/'trajectories.npz', **arrays)
    kinds = ('rms_about_zero', 'rms_about_own_mean',
             'end_cycle_1', 'end_cycle_5', 'end_cycle_10', 'max_over_all_samples')
    rankings = {option: {kind: sorted_masks(results, option, kind) for kind in kinds}
                for option in OPTIONS}
    output = {
        'status': 'Synthetic, actively driven where controllers are enabled; not physical validation.',
        'model_copy_sha256': hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
        'units_and_caveats': [
            'All quantities are synthetic model units; theta is the model angle in radians.',
            'Disabling a controller removes its entire PD force, including periodic target and controller damping; this is not a damping-only ablation.',
            'Passive spring couplings, all physical coordinates, masses, passive damping, load and timing of remaining targets are unchanged.',
            'The disturbance is one common theta and orb initial perturbation; norm mixes normalized coordinates and velocities, not a universal physical metric.',
            'Extrema and residual maxima are sampled at 240 samples/cycle; RMS is trapezoidal integration over sampled final cycle including its endpoints.',
            'Actuator work is mechanical work on scope coordinates; negative work fate is unspecified. The passive bath is bookkeeping only.',
        ],
        'parameters': asdict(P),
        'observation': {'cycles': P.cycles, 'period': 2*math.pi/P.omega,
                        'initial_state': m.initial_state(False).tolist(),
                        'perturbed_initial_state': m.initial_state(True).tolist(),
                        'perturbation_q': [0.04, 0.04, 0, 0, 0, 0, 0, 0],
                        'perturbation_v': [0]*8},
        'solver': {'method': 'DOP853', 'rtol': P.rtol, 'atol': P.atol,
                   'max_step': P.max_step, 'samples_per_cycle': P.samples_per_cycle},
        'software': {'python': platform.python_version(), 'numpy': np.__version__,
                     'scipy': scipy.__version__, 'matplotlib': matplotlib.__version__},
        'trajectory_file': {'path': 'trajectories.npz',
                            'location': 'next to ablation.json',
                            'keys': ['time'] + [f'{o}_{name}_{which}' for o in OPTIONS
                                               for name in MASKS for which in ('baseline', 'perturbed')],
                            'columns': ['q_theta', 'q_orb', 'q_fulcrum_L', 'q_fulcrum_R',
                                        'q_upper_1', 'q_upper_2', 'q_upper_3', 'q_lower_trim',
                                        'v_theta', 'v_orb', 'v_fulcrum_L', 'v_fulcrum_R',
                                        'v_upper_1', 'v_upper_2', 'v_upper_3', 'v_lower_trim',
                                        'W_plus_1', 'W_plus_2', 'W_plus_3', 'W_plus_4',
                                        'W_minus_magnitude_1', 'W_minus_magnitude_2',
                                        'W_minus_magnitude_3', 'W_minus_magnitude_4', 'Q_passive']},
        'results': results,
        'ascending_rankings_lower_to_higher': rankings,
    }
    (OUT/'ablation.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(f'Wrote {OUT/"ablation.json"} and {OUT/"trajectories.npz"}', flush=True)


if __name__ == '__main__':
    main()
