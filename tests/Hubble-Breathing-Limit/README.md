# How much can the universe still be breathing today? (Uniquea fork 7, idea 2)

**Status: Tested (provisional), November 5, 2026.** Same data, cuts and likelihoods as Hubble-Closure-Test.

**Model:** a small wobble on flat ΛCDM, periodic in ln(1+z):
E(z) = E_ΛCDM(z; Ωm) · [1 + ε(sin(2π ln(1+z)/P + φ) − sin φ)].
Joint χ² = DESI DR2 BAO + Pantheon+ (1,580 light curves). For each period P, Ωm and φ are profiled;
the 95% limit is where χ² rises 3.84 above its minimum over ε. ΛCDM alone: χ² = 1,399.80, Ωm = 0.3044.

| Period P in ln(1+z) | 0.1 | 0.15 | 0.2 | 0.3 | 0.4 | 0.6 | 0.8 | 1.0 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 95% limit on ε | 2.2% | 2.1% | 2.3% | 1.7% | 2.2% | 2.5% | 2.6% | 2.5% | 3.6% | 6.9% | 18%* | none* |
| Δχ² best wobble − ΛCDM | −0.8 | −3.4 | −5.8 | −3.5 | −2.9 | −4.0 | −6.1 | −7.3 | −6.1 | −5.9 | −6.1 | −6.3 |

\* Long periods trade against Ωm (it wanders to its bounds), so they are not constrained without a CMB prior on Ωm.

**With a CMB prior on Ωm** (Planck 2018, 0.3153 ± 0.0073; breathing_cmb_prior.py):

| P | 1.0 | 2.0 | 3.0 | 5.0 | 8.0 |
|---|---|---|---|---|---|
| 95% limit on ε | 2.4% | 2.6% | 4.5% | 12% | 29% |
| Δχ² best − ΛCDM | −6.0 | −5.0 | −4.6 | −4.2 | −4.2 |

Caveat: Planck's Ωm assumes ΛCDM; a full CMB likelihood is the next refinement.

**Result:** over periods from 0.1 to 1.5 in ln(1+z), any wobble in the expansion rate today is smaller
than about 2–3.6% (95%). The best wobbles improve χ² by at most 7.3 for three extra parameters
(ε, P, φ), about p ≈ 0.06 before accounting for scanning many periods: not a detection. It is in the
same direction as DESI DR2's reported mild preference for evolving dark energy.

**For the archive:** the sine closure's own breathing (period ≳ 2π·10.3 ≈ 65 in ln(1+z)) is far too slow to
appear as a wobble; it shows only as a change in the smooth shape, where the sine closure fit badly (3.2).

**Limitations:** profile likelihood on a coarse grid, six phase starts; no CMB; no look-elsewhere correction.
Run: `python3 breathing_limit.py <data_dir> breathing_results.json` (about 1 minute, 2 cores).
