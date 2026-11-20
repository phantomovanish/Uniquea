# Synthetic 90° prescribed motion: idealized possibility calculation

**Label:** Manus (version not recorded) Uniquea Field Observations UFO

**Status: synthetic one-axis inverse dynamics; not a design validation or a control simulation.** This is a separate calculation from the driven multi-state linear/nonlinear hybrid. It neither identifies a mechanism for that hybrid nor shows that its previous simulation could have achieved 90° under its own inputs or constraints. Conversely, a prior synthetic simulation not reaching 90° does not establish that 90° physical motion is impossible in general.

## Model and assumptions

The inputs are the invented body and below-COM pivot in [GEOMETRY_REFERENCE.md](GEOMETRY_REFERENCE.md), copied unchanged from the previous local review branch. All five masses are locked together as **one rigid body** rotating about the single fixed, horizontal y-axis. Each component is treated as a point mass at its listed COM, so the component and frame intrinsic inertias are omitted. Gravity is uniform, the pivot is ideal and externally supported, and the only modeled resistance is rotational viscous damping. The external motor torque is allowed to follow the calculated function exactly, with sufficient positive and negative torque and with the initial state set exactly as prescribed. No control law, drive limits or tracking performance are assumed or demonstrated.

- Total mass: `M = 39 kg`; gravity: `g = 9.81 m/s²`.
- Combined COM: `(3/325, -2/325, 23/52) m`; fixed pivot: `(3/325, -2/325, -1/4) m` in the initial body frame.
- COM height above pivot: `h = 9/13 m`; hence `Mgh = 264.87 N·m`.
- Point-mass inertia about y: `I = Σ mᵢ[(xᵢ − xₚ)² + (zᵢ − zₚ)²] = 1041621/32500 = 32.04987692307692 kg·m²`. The source displays the rounded value 32.050; this calculation retains the unrounded value.
- `b = 5 N·m·s/rad`, amplitude `A = π/2 rad`, imposed **phase rate** `ω = 0.15 s⁻¹`. This rate is not the peak angular speed.
- The trajectory is `θ(t) = A sin(ωt)`. Its endpoints are **−90° and +90°**, giving **180° end-to-end travel**, not a total 90° sweep. At the selected `t = 0`, `θ = 0` but `θ̇ = Aω = 0.23561944901923448 rad/s`; it does **not** start from rest.

All masses, dimensions, damping and phase rate are synthetic; these figures are not measurements or equipment specifications.

## Equation and outputs

With positive `θ` in the direction of the below-COM gravity moment, `τ_g = +Mgh sin θ` and `U = Mgh cos θ`. The rotational balance is `Iθ̈ = τ_motor + τ_g − bθ̇`, so

```text
θ̇ = Aω cos(ωt)              θ̈ = −Aω² sin(ωt)
τ_motor = Iθ̈ − Mgh sin θ + bθ̇.
```

Using this prescribed trajectory and torque, the numerical results for one full cycle are:

| Quantity | Result |
|---|---:|
| Period | 41.88790204786391 s |
| Angle minimum / maximum | −90° / +90° |
| Peak absolute angular speed | 0.23561944901923448 rad/s = 13.5°/s |
| Peak absolute angular acceleration | 0.035342917352885174 rad/s² = 2.025°/s² |
| Motor torque minimum / maximum | −266.12509300047674 / +266.12509300047674 N·m |
| Peak absolute motor torque | 266.12509300047674 N·m |
| Positive actuator mechanical work `W+` | +534.4313619702378 J |
| Negative actuator mechanical work `W−` | −528.6176850926815 J |
| Net actuator mechanical work `W+ + W−` | +5.813676877556304 J |
| Viscous loss | +5.813676877556215 J |
| Change in `K + U` over cycle | 0.0 J |
| Energy residual `W+ + W− − Δ(K+U) − loss` | +8.881784197001252×10⁻¹⁴ J |

Numbers are binary floating-point outputs, not exact physical precision. The peak torque is found by solving the torque derivative for its extrema, not merely by evaluating the turning points. The negative minimum occurs at phase `1.7163071911673102 rad`, when `θ ≈ +89.04887646474323°`; the positive maximum occurs one half-cycle later. Over the cycle `W+` integrates intervals with positive actuator mechanical power, and `W−` is the **signed** integral where that power is negative. Negative work requires the actuator or its load to absorb energy; it does not guarantee electrical regeneration. The net positive work compensates damping, so this is **not** a no-energy-input motion.

## Independent analytic ledger

Write `φ = ωt`, `C = Mgh`, `q = bA²ω`, and

```text
E(φ) = ½ I(Aω cos φ)² + C cos(A sin φ)
τ(φ) = −I Aω² sin φ − C sin(A sin φ) + bAω cos φ
∫ bθ̇² dt = q [φ/2 + sin(2φ)/4] + constant.
```

The instantaneous balance is `τθ̇ = dE/dt + bθ̇²`. Over a full cycle `ΔE = 0` exactly and

```text
W_net = loss = bA²ωπ = 3π³/16 J.
```

To check the **separate** signed work integrals, let `r` be the unique zero of `τ(φ)` in `(0, π/2)`, numerically `0.0028238912992171665 rad`, and set `s = r + sin(2r)/2`. The other zero is `π + r`; the velocity zeros are `π/2` and `3π/2`. Integrating the energy identity over the resulting positive- and negative-power intervals gives

```text
W+ =  2E(r) + q(π/2 + s)
W− = −2E(r) + q(π/2 − s).
```

These are exact expressions given the root `r`; its value is computed numerically. Compared with independent adaptive quadrature split at those power sign changes, the differences are approximately `−1.14×10⁻¹³ J` for `W+`, `+2.27×10⁻¹³ J` for `W−`, and zero at displayed precision for the viscous loss. The quadrature-reported sum of error estimates for the work pieces is about `1.18×10⁻¹¹ J`. These figures are numerical checks, not uncertainty bounds on real equipment.

For the extrema, `τ(φ + π) = −τ(φ)`. The derivative is strictly negative on `(0, π/2)`. On `(π/2, π)` there is a unique minimum: writing `y = π − φ`, its stationarity condition is `[IAω² + CA cos(A sin y)] cot y = bAω`, whose left-hand side decreases from infinity to zero. Symmetry gives the global maximum in the other half-cycle.

## What this does and does not show

This is a constructive **mathematical** example of reaching and reversing at ±90° for the assumed fixed-axis rigid-body equation, if the specified time-varying torque and initial conditions are supplied. Under the public Master's HAT rubric, the invented geometry and drive are **Postulated** assumptions and the checkable values are **Derived** inside that model, not measured observations. It is not a dynamically controlled simulation, and it does not establish feasibility of a 3D four-scope mechanism, motor sizing, real travel or clearance, bearings or external reactions, structural strength, safety, control or stability against perturbations, a probability of success, or avoidance of gimbal lock. In particular, the below-COM upright configuration is gravity-unstable; prescribing an exact path does not stabilize departures from it. More detailed geometry, actual inertias and a validated actuator, power, controls and safety analysis would be needed for any physical system.

## Reproduce

From this directory, using Python 3 and SciPy (`requirements.txt`; captured with Python 3.12.3 and SciPy 1.18.1):

```bash
python3 -B calculate_ninety.py      # prints and rewrites results.json in this directory
python3 -B -m unittest -v           # seven independent consistency tests
```

`results.json` records the numerical outputs and analytic comparisons. `test_results.txt` captures the completed test run. The calculation reads no network data and is self-contained after copying the synthetic inputs from the source note.
