# Hybrid rocker experiment — synthetic, externally driven

This is a **prospective numerical toy experiment**, not a physical or observational validation. It compares five fixed blends of the **newer local** linear rocker law A and nonlinear-transition rocker law B in the existing four-scope surrogate. Both endpoints and all three intermediate laws use the same nonzero damping, external controls, model parameters, initial state and 20-cycle solver settings. The underlying model is an independently constructed surrogate, not an archived field equation or a three-dimensional gimbal.

**Do not confuse the forks:** the earlier `Orb-Adiabatic-Check` A/B denotes *no damping* versus *damping with energy booked to an Alt Field*. That is a different model and question. The present A/B denotes the **rocker restoring torque**. Both current cases have damping; their bath is a cumulative bookkeeping variable, not a modeled physical Alt Field or heat reservoir. See [FORK_LEDGER.md](FORK_LEDGER.md).

## Exact blend and transition

For constant blend weight `λ ∈ {0, 0.25, 0.5, 0.75, 1}`, with `k=1`, gain `g=1.5`, angle `θ`, `r=|θ|`, limits `a=0.10`, `b=0.24`, `w=b−a`, and `h=clamp((r−a)/w,0,1)`, set `S(h)=6h⁵−15h⁴+10h³`:

```text
FA = −k θ                      UA = ½ k θ²
FB = −k θ [1 + g S(h)]         UB = UA + k g I(r)
Fλ = (1−λ) FA + λ FB          Uλ = (1−λ) UA + λ UB
   = −k θ [1 + λg S(h)]           = ½ k θ² + λkg I(r)
```

The reference's exact integral is

```text
I(r) = w [a (2.5h⁴−3h⁵+h⁶) + w (2h⁵−2.5h⁶+(6/7)h⁷)]
       + ½ max(r²−b², 0).
```

Thus `Fλ=−dUλ/dθ`; it is not an ad-hoc force-only mixture. **B already equals A at `|θ|≤0.10` rad**, becomes nonlinear in the smooth `0.10–0.24` rad transition, and is linear in θ again above `0.24` with stiffness factor 2.5 rather than 1. The upper-branch stiffness factor of a hybrid is `1+1.5λ`. At λ=0 there is no nonlinear transition; occupancy relative to the same thresholds is shown only for comparison. Force and its first two derivatives match at the boundaries. All coordinates, coefficients and time are dimensionless model quantities except that θ is interpreted as radians.

Each simulation keeps λ **fixed in time**. If an external mechanism changes it, the mechanical energy identity acquires the additional parameter-work term `λdot (UB−UA)`. No switching schedule, λ-actuator, or such work is simulated here; a time-varying blend cannot be inferred from this fixed-blend comparison.

## Dynamics and accounting

The dynamic coordinates and velocities are `q=(θ,y,fL,fR,s1,s2,s3,s4)` and `qdot`, for rocker, orb, two fulcra, three upper scopes and lower trim. Spring couplings act in both directions, not as prescribed orb/fulcrum motions. The inherited common potential is

```text
U = Uλ + Lθ + ky y²/2 + kq(fL²+fR²)/2
    + Σj kf,j(fj−bjθ)²/2 + Σj ko,j(y−fj)²/2
    + Σi ks,i(si−aiθ)²/2.
M qddot = −∇U − C qdot + Fcontroller(t,q,qdot).
```

The reference coefficients are: masses `(2,1,0.6,0.6,0.3,0.3,0.3,0.3)`; viscous damping `(0.18,0.16,0.18,0.18,0.12,0.12,0.12,0.12)`; `L=0.28`, `ky=0.30`, `kq=0.22`; fulcrum levers `(-0.7,+0.7)`, `kf=(0.8,0.8)`, `ko=(0.8,1.15)`; scope levers `(1,1,1,1)` and `ks=(1.3,1.3,1.3,0.5)`. External controls act only on the four scopes:

```text
ri = biasi + amplitudei sin(0.75 t + phasei)
Fi = kp,i (ri−si) − kd,i sdot,i
```

`bias=(0.25,0.25,0.25,0)`, `amplitude=(0.35,0.35,0.35,0.035)`, `phase=(0,π/2,π,3π/2)`, `kp=(4,4,4,16)`, and `kd=(0.65,0.65,0.65,0.8)`. Targets are time-phase markers, not physical gimbal angles. The conservative load's potential is already included in `H`, so it is not also booked as separate work.

`H` is kinetic plus potential energy. For each of the four controls, signed instantaneous mechanical power is `pi=Fi sdot_i` and the integrated ledger records `W+i=∫max(pi,0)dt`, `W−i=∫max(−pi,0)dt` and passive loss `Q=∫Σ Cj qdot_j²dt`. The checked balance is

```text
H(t)−H(0) = Σ W+i − Σ W−i − Q.
```

Controller damping is included in signed controller work, not counted again as passive `Q`. Negative work is extracted from the model; what happens to it is not specified. The oscillations are supplied by the external controls, not by a demonstrated self-powering mechanism. Work is not electrical consumption or efficiency.

## Method and output

The reference initial positions are `(0.06,0.035,−0.015,0.025,0.20,0.20,0.20,0)`; velocities and all ledger variables start at zero. The primary pair changes both `θ0` and `y0` by `+0.04`, as in the reference. The extra mirror pair changes both by `−0.04`. Each pair uses the **same** drive as its own baseline; differences are finite-time deterministic responses, not probability estimates or a full stability survey.

We use the reference 20 periods of `2π/0.75`, 240 samples per period (4,801 samples), DOP853 `rtol=1e−9`, `atol=1e−11`, `max_step=0.08`. Baseline and `+` perturbation at every λ are repeated with `rtol=2e−11`, `atol=2e−13`, `max_step=0.04`; direct endpoint runs independently invoke the unchanged original integrator. The protocol and predetermined checks are in [PROTOCOL.md](PROTOCOL.md). [RESULTS.md](RESULTS.md) reports the observed values and limitations.

- `reference_model.py` is a byte-for-byte copy of the archived source. `hybrid_experiment.py` is the new standalone experiment; `test_hybrid_experiment.py` provides mathematical, simulation and saved-output tests. `requirements.txt` records the exact installed versions used.
- `output/trajectories.npz` contains complete sampled arrays for five λ × five runs: `baseline_standard`, `baseline_tight`, `plus_standard`, `plus_tight`, `minus_standard`; two additional arrays are direct reference endpoint reruns. The `time` array has shape `(4801,)`; each state array has shape `(4801,25)`. Array keys are `lambda_NNN__VARIANT`, e.g. `lambda_025__plus_tight`. Load with `np.load(..., allow_pickle=False)`. Column names are in `output/summary.json`.
- Five `output/lambda_NNN_baseline.csv` files additionally give time, all 25 state/ledger columns, energy and residual. `output/metrics.csv`, `output/stages.csv` (400 quarter-cycle records), and `output/summary.json` contain summaries, per-controller works, perturbation norms, numerical comparisons and parameters. `output/overview.png` is a plot. `output/SHA256SUMS.txt` is a local integrity inventory, not a project-standard manifest.

## Reproduce offline

Python 3.12.3 and the exact versions in `requirements.txt` were used. With those packages already available locally, run:

```sh
cd /home/ubuntu/gimbal_hybrid_work
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg python3 hybrid_experiment.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest -v
(cd output && sha256sum -c SHA256SUMS.txt)
```

The first command **replaces generated files under this work directory's `output/`**. No internet access is needed or attempted by the code. If packages must be installed, provide matching local wheels and use `pip --no-index --find-links LOCAL_WHEEL_DIRECTORY -r requirements.txt`; no download is built into this workflow. Different software/platform arithmetic can change the last digits, and a bit-for-bit checksum is guaranteed only for the saved files as delivered, not across environments.

## Limitations and source preservation

No measured quantities, observations, calibration, independent experimental replication, actuator limits, impacts, 3D mechanism, physical bath or explicit scale mapping are present. Sampled extrema and threshold occupancy have sampling error; solver reruns provide numerical comparisons rather than certified bounds. Short finite-time perturbations and one nearly repeating final cycle do not establish a limit cycle or global stability. The House Rule 11 text and actual HAT rubric were not supplied or identified in this working scope; no HAT category, archive filing or publication is asserted. This isolated package is staged only, and the original repository was not edited, committed or published.
