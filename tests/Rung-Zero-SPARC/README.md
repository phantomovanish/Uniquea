# Uniquea test 1: one a₀ against galaxy rotation curves (rung zero)

**Status: Tested (provisional), November 9, 2026; revised the same day after review by Phantom.**

## The relation tested
Rung zero of the Nash ladder is the deep-MOND law, F = m·a²/a₀ (credit: Milgrom, 1983). In the low-acceleration
regime it predicts

```
g_obs = √(a₀ · g_bar)
```

- **g_bar** is the Newtonian acceleration computed from the baryons alone (gas, stellar disk, bulge):
  g_bar = (V_gas|V_gas| + Υ_disk V_disk² + Υ_bulge V_bulge²) / R.
- **g_obs** is the acceleration inferred from the measured rotation curve: g_obs = V_obs² / R.

## Data
SPARC (Lelli, McGaugh & Schombert 2016): 175 disk galaxies with Spitzer 3.6 µm photometry and HI/Hα
rotation curves. The official SPARC site is not reachable from this workspace, so the files were taken from
public GitHub copies (Rotmod_LTG: Cisneros-Galaxy/RCFM, commit 4f9e239; galaxy table: JuliaAstroSim/AstroIC.jl,
commit 79c9997). The official SPARC page (astroweb.case.edu/SPARC) is online and provides the sample and
mass-model data; it was not reachable from this analysis workspace. Spot check: NGC 2403 distance 3.16 Mpc and
inclination 63° match SPARC. **Next step: rerun against the official files.**

**Selection, reconciled with the published count.** The published radial-acceleration-relation analysis
(McGaugh, Lelli & Schombert 2016, *PRL* 117, 201101) reports 2,693 points in 153 galaxies. Applying its cuts
step by step to these files (point-level mask: point_mask.json):

| Step | Points | Galaxies |
|---|---|---|
| All SPARC rotation-curve points | 3,390 | 175 |
| Quality Q < 3 | 3,270 | 163 |
| Inclination ≥ 30° | 3,167 | 153 |
| Positive g_bar (used in this test's first run) | 3,165 | 153 |
| Velocity uncertainty below 10% (the published point cut) | 2,695 | 147 with points left |

The first run omitted the 10% velocity-error cut, which accounts for the difference (472 points). With it, the
count is 2,695 against the published 2,693; the remaining 2 are not identified here.
Mass-to-light ratios fixed at Υ_disk = 0.5, Υ_bulge = 0.7. Errors: rotation-velocity errors plus a 0.1 dex floor
for distance, inclination and mass-to-light systematics.

## Results
Both selections are shown. Uncertainties marked "stat" are statistical under this fit only.

| Test | First run (3,165 points) | Published cut (2,695 points) |
|---|---|---|
| Deep regime, one a₀ in g_obs = √(a₀ g_bar) | a₀ = 1.31 × 10⁻¹⁰ m/s² (± 0.02 stat; 1,240 points) | a₀ = 1.39 × 10⁻¹⁰ m/s² (957 points) |
| All points, standard RAR function | a₀ = 1.12 × 10⁻¹⁰ m/s² | a₀ = 1.17 × 10⁻¹⁰ m/s² (published: 1.20 ± 0.02 stat ± 0.24 sys) |

**Systematic check (mass-to-light ratio, published cut):** Υ_disk = 0.4, 0.5, 0.6 (Υ_bulge = 1.4 Υ_disk) gives
deep-regime a₀ = 1.60, 1.39, 1.22 × 10⁻¹⁰ and RAR a₀ = 1.54, 1.17, 0.93 × 10⁻¹⁰. So the total uncertainty is
**roughly ± 15–30%**, far larger than the statistical ± 0.02, in line with the published ± 0.24. Distance and
inclination are galaxy-level nuisance parameters too; this test does not marginalize them.

### Detail from the first run
| Test | Result |
|---|---|
| Deep regime only (g_bar < 1.2 × 10⁻¹¹ m/s², 1,240 points), rung zero g = √(a₀ g_bar) | **a₀ = (1.31 ± 0.02) × 10⁻¹⁰ m/s²** (statistical); mean residual −0.025 dex, scatter 0.20 dex |
| All points, the standard RAR function g_bar / (1 − e^(−√(g_bar/a₀))) | a₀ = (1.12 ± 0.01) × 10⁻¹⁰ m/s² (published: 1.20 ± 0.02 stat ± 0.24 sys) |
| One a₀ per galaxy (92 galaxies with ≥ 5 deep points) | log a₀ = −9.92, spread 0.35 dex at fixed mass-to-light ratios |
| Archive's lead cH₀/2π (H₀ = 67.4) | 1.04 × 10⁻¹⁰ m/s²: within 7–20% of the fitted values, not equal |

## Reading
- **Ensemble result and universality are different claims.** A shared relation across the whole sample does
  not prove that one a₀ holds for every galaxy.
- **One a₀ describes the ensemble.** In the deep regime the points follow the square-root law with a
  single a₀ near 1.2–1.3 × 10⁻¹⁰ m/s², with no offset in the mean. The test as stated in the audit packet does
  not fail at this level.
- **Universality galaxy by galaxy is not settled here.** With fixed mass-to-light ratios the per-galaxy a₀
  values spread by a factor of about 2. The literature disputes whether that spread is real (Rodrigues et al.
  2018) or comes from distance, inclination and mass-to-light uncertainties (Li et al. 2018, who fit those per
  galaxy). This test does not marginalize them, so it cannot decide.
- **Rung zero inherits MOND's record.** It shares MOND's success on rotation curves and its known problems:
  galaxy clusters, the Bullet Cluster and the cosmic microwave background (open in the audit packet).
- **cH₀/2π stays a lead**, not a derivation: close, but the fitted a₀ is 7–20% higher.

## Run
```
python3 rung_zero_sparc.py <Rotmod_LTG dir> <SPARC_Lelli2016c.mrt> results.json
python3 plot_rar.py <Rotmod_LTG dir> <SPARC_Lelli2016c.mrt> results.json rung_zero_sparc.png
```

## For kids (Stevo's version)
Galaxies spin like merry-go-rounds. The stars at the edge should go slower, but they don't. Either invisible
stuff pulls them (dark matter), or gravity works a little differently when it is very weak. We tested the
second idea on 153 real galaxies, and one number made the rule work for the whole bunch, like one recipe that
bakes 153 different cakes right. Not a win yet, but it passed a real test from the sky.

*En español:* Probamos una idea sobre la gravedad con 153 galaxias reales, y un solo número funcionó para todas
juntas. No es la última palabra, pero pasó la prueba del cielo.
