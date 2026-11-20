# Damping transition and scope-control checks

**Label:** Manus (version not recorded) Uniquea Field Observations UFO

**Status:** Derived inside the specified synthetic models; physical interpretation Open. The test inputs below are chosen examples, not measured properties. **UFO number: pending.**

## An orthogonal damping probe

The earlier archive's Option A/B fork concerns **no damping versus damping**, whereas the local rocker-law A/B compares **linear versus nonlinear restoring force**. To avoid silently merging the two, `damping_transition_probe.py` independently varies the passive damping of the orb coordinate for each rocker law. All other parameters, initial conditions, 20-cycle external drive and passive damping remain fixed.

The four damping rates are γ = 0, 0.04, 0.16 and 0.40 per model-time unit. With the drive period T_drive = 2π/0.75, the products γT_drive are approximately 0, 0.3351, 1.3404 and 3.3510. The older archived fork's T is a *contraction timescale* in another model; these numbers cannot be substituted into its tipping equations. Even at γ = 0 here, other passive damping and controller damping remain. This is not a fully undamped system or a replication of the archived fork.

For one initial disturbance, the baseline-to-perturbed physical-state norm starts at 0.05657. At the end of cycle 5, the observed norms for the four rates are:

| Orb damping rate γ | Linear rocker A | Nonlinear rocker B |
|---:|---:|---:|
| 0 | 0.0040081 | 0.0038759 |
| 0.04 | 0.0023978 | 0.0022827 |
| 0.16 | 0.0008830 | 0.0006632 |
| 0.40 | 0.0007587 | 0.0004236 |

For these four chosen cases, the sampled norm at cycle 5 decreases as the orb damping rate rises. This observation is not proof of a general monotonic law. The norm combines model-coordinate positions and velocities; it is not a unique physical measure of wobble. Some trajectories first amplify the disturbance. One initial disturbance, one parameter set and a finite run do not establish global or asymptotic stability. Late-cycle values around 1e-9 approach numerical sensitivity.

Passive dissipated energy is recorded in the existing **bookkeeping bath**; no independent physical Alt Field, coupling or reservoir dynamics is modeled. The controller work is tracked separately. Across the baseline damping cases the largest sampled energy-ledger residual is about 1.84e-9 model-energy units. Solver comparisons of gross work are less accurate than the small net residual; the full precision results and convergence diagnostics are in `output/damping_transition.json`.

## Which controllers help?

`ablation.py` compares four configurations: all four controllers on; the three upper controls off with the lower on; the lower off with the three uppers on; and all four off. It runs both rocker laws for ten cycles, with the same initial conditions and one paired perturbation. In each ablation the entire selected proportional–derivative controller is removed, **including its periodic drive**; spring links and passive damping remain. Therefore this does not isolate damping alone.

In this chosen example, the largest sampled difference between paired trajectories is about 0.0926 with all four controls on, versus 0.1434 with all controls off, for both laws. But the all-off system loses its periodic drive and settles near a shifted angle, so its lack of oscillation about its own mean is not evidence of accurate alignment or a preserved breath. The detailed ranking depends on the measurement: at cycle 10, the all-off linear case has a smaller perturbation norm than the all-on case.

Comparing upper-only to all-on, the lower controller reduces the final-cycle RMS angle **about each run's own mean** by approximately 21.9% for the linear law and 11.3% for the nonlinear law. In the same comparison, its cycle-10 perturbation norm is larger by approximately 74.6% and 14.0%, respectively. This illustrates a trade-off in this toy model; it does not establish that the three upper scopes remove most wobble, that the lower scope is a physical gimbal, or that either design is generally superior. A change in the lower controller changes both its feedback and its drive.

The full sampled trajectories are in `output/trajectories.npz`; `output/ablation.json` includes per-controller work, passive losses, baselines and perturbations, settings, rankings and sampled energy residuals. The worst recorded sampled ablation ledger residual is about 1.18e-8 model-energy units. Eight ablation unit tests and nine damping-probe unit tests were run in addition to the nine main model tests.

## Reproduce

From the repository root, with the dependencies in `requirements.txt` installed:

```sh
cd tests/Four-Way-Gimbal-Transition
OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python3 -m unittest -v test_model test_damping_probe test_ablation
OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python3 damping_transition_probe.py
OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python3 ablation.py
```

The last two commands replace their respective files under `output/`. Work on a copy if retaining the recorded run. All files are model inputs or outputs, not measured cosmological or hardware data. A test on a real mechanism needs defined mass and linkage geometry, sensor and actuator limits, a stated disturbance metric, and independent measurements.

**Not claimed:** observational validation, a real 90° gimbal flip, a sustained unpowered oscillation, an experimentally established physical Alt Field, a unique optimally stable design, or a deduction of the cosmos from these synthetic examples.
