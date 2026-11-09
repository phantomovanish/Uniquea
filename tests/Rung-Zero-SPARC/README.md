# Uniquea test 1: one a₀ against galaxy rotation curves (rung zero)

**Status: Tested (provisional), November 9, 2026.** Rung zero of the Nash ladder is the deep-MOND law,
F = m·a²/a₀, so in the low-acceleration regime it predicts g_obs = √(a₀·g_bar). Credit: Milgrom (1983).

## Data
SPARC (Lelli, McGaugh & Schombert 2016): 175 disk galaxies with Spitzer 3.6 µm photometry and HI/Hα
rotation curves. The official SPARC site is not reachable from this workspace, so the files were taken from
public GitHub copies (Rotmod_LTG: Cisneros-Galaxy/RCFM, commit 4f9e239; galaxy table: JuliaAstroSim/AstroIC.jl,
commit 79c9997). Spot check: NGC 2403 distance 3.16 Mpc and inclination 63° match SPARC. A rerun against the
official files is recommended.

Cuts (as McGaugh, Lelli & Schombert 2016): quality Q < 3, inclination ≥ 30° → **153 galaxies, 3,165 points**.
Mass-to-light ratios fixed at Υ_disk = 0.5, Υ_bulge = 0.7. Errors: rotation-velocity errors plus a 0.1 dex floor
for distance, inclination and mass-to-light systematics.

## Results
| Test | Result |
|---|---|
| Deep regime only (g_bar < 1.2 × 10⁻¹¹ m/s², 1,240 points), rung zero g = √(a₀ g_bar) | **a₀ = (1.31 ± 0.02) × 10⁻¹⁰ m/s²** (statistical); mean residual −0.025 dex, scatter 0.20 dex |
| All points, the standard RAR function g_bar / (1 − e^(−√(g_bar/a₀))) | a₀ = (1.12 ± 0.01) × 10⁻¹⁰ m/s² (published: 1.20 ± 0.02 stat ± 0.24 sys) |
| One a₀ per galaxy (92 galaxies with ≥ 5 deep points) | log a₀ = −9.92, spread 0.35 dex at fixed mass-to-light ratios |
| Archive's lead cH₀/2π (H₀ = 67.4) | 1.04 × 10⁻¹⁰ m/s²: within 7–20% of the fitted values, not equal |

## Reading
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
