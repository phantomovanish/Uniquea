#!/usr/bin/env python3
"""Reproduce a four-point, orb-only damping probe in the synthetic model.

Run from this directory, without network access:
    OPENBLAS_NUM_THREADS=1 python3 -B damping_transition_probe.py

The four gamma values are decay-rate coefficients in inverse model-time units:
    0.0, 0.04, 0.16, 0.40.
The orb mass is 1 in the unmodified model, and its passive force is
-c_y*y_dot = -m_y*gamma*y_dot.  The sweep replaces only c_y.  This script's
local A/B labels refer to linear/nonlinear rocker laws, NOT the archived
no-damping/damping A/B fork.  Zero orb damping is not an undamped system.

The probe gamma*T_drive uses T_drive = 2*pi/omega.  It is not the archived
contraction product gamma*T, does not reproduce that model or its scaling, and
is not evidence for a physical Alt Field.  Bath and controller bookkeeping
are the ones already in model.py.  Results are finite-time synthetic behavior,
not an asymptotic stability proof or an observational validation.
"""
from __future__ import annotations

import json
import math
import platform
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
import numpy as np
import scipy

import model

GAMMAS = (0.0, 0.04, 0.16, 0.40)
OPTIONS = ("A", "B")
OUTPUT = Path(__file__).resolve().parent / "output" / "damping_transition.json"
TIGHT_SOLVER = {"rtol": 2e-11, "atol": 2e-13, "max_step": 0.04}


def parameters_for_gamma(gamma: float, base: model.Parameters = model.P) -> model.Parameters:
    """Replace only the orb's passive damping coefficient c_y = m_y*gamma."""
    gamma = float(gamma)
    if not math.isfinite(gamma) or gamma < 0:
        raise ValueError("gamma must be finite and nonnegative")
    if not math.isfinite(base.masses[1]) or base.masses[1] <= 0:
        raise ValueError("orb mass must be finite and positive")
    damping = list(base.damping)
    damping[1] = gamma * base.masses[1]
    return replace(base, damping=tuple(damping))


def sample_rms(t: np.ndarray, values: np.ndarray) -> float:
    """Trapezoidal time-weighted RMS, rather than a double-counted endpoint mean."""
    return float(np.sqrt(np.trapezoid(values * values, t) / (t[-1] - t[0])))


def evaluate_case(option: str, p: model.Parameters, *, convergence: bool = True) -> dict:
    """Evaluate one rocker law and gamma with common initial and drive settings.

    The perturbed run changes theta(0) and orb y(0) by +0.04 each.  All other
    initial values, parameters, time samples, and normal solver settings match.
    Optional tighter baseline rerun is a numerical sensitivity check, not an
    error bound.  Only the original total passive bath is integrated by the ODE;
    the orb contribution is separately approximated by sampled trapezoidal rule.
    """
    if option not in OPTIONS:
        raise ValueError("option must be A or B")
    t, z = model.run(option, p)
    tp, zp = model.run(option, p, perturb=True)
    if not np.array_equal(t, tp):
        raise RuntimeError("baseline and perturbed sample grids differ")
    if convergence:
        tt, tight = model.run(option, p, **TIGHT_SOLVER)
        if not np.array_equal(t, tt):
            raise RuntimeError("normal and tighter sample grids differ")
        convergence_metrics = {
            "tight_solver": TIGHT_SOLVER.copy(),
            "max_abs_physical_state_difference": float(np.max(np.abs(z[:, :16] - tight[:, :16]))),
            "max_abs_gross_ledger_state_difference": float(np.max(np.abs(z[:, 16:] - tight[:, 16:]))),
            "final_physical_state_difference_norm": float(np.linalg.norm(z[-1, :16] - tight[-1, :16])),
        }
    else:
        convergence_metrics = None

    H, residual = model.ledger(z, option, p)
    Hp, residual_p = model.ledger(zp, option, p)
    diff = np.linalg.norm(zp[:, :16] - z[:, :16], axis=1)
    sp = p.samples_per_cycle
    final = z[-sp - 1:]
    tfinal = t[-sp - 1:]
    abs_theta = np.abs(final[:-1, 0])  # exclude duplicated phase-zero endpoint
    Wplus = z[-1, model.POS]
    Wminus = z[-1, model.NEG]
    orb_dissipation_power = p.damping[1] * z[:, 9] ** 2
    total_dissipation_power = np.sum(np.asarray(p.damping) * z[:, 8:16] ** 2, axis=1)
    orb_loss_sampled = float(np.trapezoid(orb_dissipation_power, t))
    passive_loss_sampled = float(np.trapezoid(total_dissipation_power, t))
    gamma = p.damping[1] / p.masses[1]
    T_drive = 2.0 * math.pi / p.omega
    return {
        "gamma": float(gamma),
        "orb_damping_coefficient": float(p.damping[1]),
        "gamma_times_T_drive": float(gamma * T_drive),
        "local_rocker_option": option,
        "rocker_law": "linear" if option == "A" else "nonlinear smooth hardening",
        "parameters": asdict(p),
        "energy_ledger": {
            "H_initial": float(H[0]),
            "H_final": float(H[-1]),
            "delta_H": float(H[-1] - H[0]),
            "positive_actuator_work_by_controller": Wplus.tolist(),
            "negative_actuator_work_magnitude_by_controller": Wminus.tolist(),
            "total_positive_actuator_work": float(np.sum(Wplus)),
            "total_negative_actuator_work_magnitude": float(np.sum(Wminus)),
            "net_actuator_work": float(np.sum(Wplus - Wminus)),
            "passive_loss_to_bookkeeping_bath": float(z[-1, model.LOSS]),
            "orb_passive_loss_sampled_trapezoid": orb_loss_sampled,
            "all_passive_loss_sampled_trapezoid": passive_loss_sampled,
            "sampled_total_minus_ode_integrated_passive_loss": float(passive_loss_sampled - z[-1, model.LOSS]),
            "final_ledger_residual": float(residual[-1]),
            "max_abs_sampled_ledger_residual": float(np.max(np.abs(residual))),
            "perturbed_delta_H": float(Hp[-1] - Hp[0]),
            "perturbed_net_actuator_work": float(np.sum(zp[-1, model.POS] - zp[-1, model.NEG])),
            "perturbed_passive_loss_to_bookkeeping_bath": float(zp[-1, model.LOSS]),
            "perturbed_max_abs_sampled_ledger_residual": float(np.max(np.abs(residual_p))),
        },
        "final_cycle_sampled_response": {
            "theta_min": float(np.min(final[:, 0])),
            "theta_max": float(np.max(final[:, 0])),
            "theta_time_weighted_rms": sample_rms(tfinal, final[:, 0]),
            "orb_min": float(np.min(final[:, 1])),
            "orb_max": float(np.max(final[:, 1])),
            "orb_time_weighted_rms": sample_rms(tfinal, final[:, 1]),
            "fraction_below_or_equal_lower_threshold": float(np.mean(abs_theta <= p.transition_lo)),
            "fraction_between_thresholds": float(np.mean((abs_theta > p.transition_lo) & (abs_theta < p.transition_hi))),
            "fraction_above_or_equal_upper_threshold": float(np.mean(abs_theta >= p.transition_hi)),
            "physical_state_closure_norm": float(np.linalg.norm(final[-1, :16] - final[0, :16])),
        },
        "perturbation": {
            "initial_physical_state_difference_norm": float(diff[0]),
            "difference_norm_at_cycle_ends": {
                str(cycle): float(diff[cycle * sp])
                for cycle in (1, 5, 10, 20) if cycle <= p.cycles
            },
            "maximum_sampled_difference_norm": float(np.max(diff)),
            "final_cycle_maximum_sampled_difference_norm": float(np.max(diff[-sp - 1:])),
            "final_sampled_difference_norm": float(diff[-1]),
        },
        "numerical_comparison": convergence_metrics,
    }


def build_output() -> dict:
    p = model.P
    T_drive = 2.0 * math.pi / p.omega
    cases = []
    for gamma in GAMMAS:
        pg = parameters_for_gamma(gamma, p)
        for option in OPTIONS:
            print(f"gamma={gamma:.2f}, gamma*T_drive={gamma*T_drive:.9g}, rocker={option} ({'linear' if option == 'A' else 'nonlinear'})", flush=True)
            cases.append(evaluate_case(option, pg))
    return {
        "status": "synthetic, actively driven, finite-time, unvalidated",
        "description": "Four-point probe, not a continuous sweep or interpolation of archived A/B models.",
        "gamma_values_inverse_model_time": list(GAMMAS),
        "gamma_definition": "gamma=c_orb/m_orb, where passive orb force is -c_orb*y_dot and m_orb=1; no other parameter is varied",
        "T_drive": T_drive,
        "dimensionless_probe_definition": "gamma*T_drive = gamma*(2*pi/omega), not archived gamma*T",
        "archived_comparison_caution": "Archived T is the contraction timescale in a different system; these products and the archived shift scaling are not interchangeable.",
        "naming_caution": "Local A/B are linear/nonlinear rocker laws; archived A/B denote no damping/damping. gamma=0 removes only orb-coordinate passive damping, not other damping or controller losses.",
        "bath_caution": "Damping is integrated into the existing passive bookkeeping bath, not a modeled physical Alt Field. Negative actuator work represents extracted mechanical energy; its destination is unspecified.",
        "limitations": [
            "Only four gamma values, one synthetic parameter set and one specified perturbation are tested.",
            "Dynamically coupled driven coordinates need not behave like an isolated damped orb, and no archived tipping or fold is modeled.",
            "20-cycle results, final-cycle closure and perturbation contraction are not proofs of asymptotic stability, a basin of attraction, or self-sustaining motion.",
            "Extrema and threshold occupancy are sampled; threshold occupancy for the linear law is only a reference to the nonlinear law's thresholds.",
            "The trapezoidal orb and total damping estimates use sampled states; only the total bath is integrated as a state. Their difference is a quadrature diagnostic, not a separate reservoir.",
            "The physical-state Euclidean norm mixes normalized coordinates and velocities and is not a unique physical distance. Small late differences can approach numerical resolution.",
            "No measured data, physical field, fitted damping, experimental validation, Floquet analysis or rigorous solver-error bounds are included."
        ],
        "baseline_parameters": asdict(p),
        "initial_augmented_state": model.initial_state().tolist(),
        "perturbed_initial_augmented_state": model.initial_state(perturb=True).tolist(),
        "solver": {
            "method": "DOP853",
            "cycles": p.cycles,
            "samples_per_cycle": p.samples_per_cycle,
            "rtol": p.rtol,
            "atol": p.atol,
            "max_step": p.max_step,
            "tighter_baseline_comparison": TIGHT_SOLVER,
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "cases": cases,
    }


def main() -> None:
    results = build_output()
    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
