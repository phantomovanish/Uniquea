# Hubble-closure test: reproducibility package

**Status: Tested (provisional).** Run November 2, 2026, and re-run from scratch November 3, 2026.
The results matched exactly.

## What is tested
The archive's sine closure (Master §3.2, B44) is compared with flat ΛCDM:

- **Sine closure:** cos θ = cos θ₀ + ln(1+z)/A and H/H₀ = sin θ / sin θ₀. Two free parameters, A and θ₀.
- **Flat ΛCDM:** one free parameter, Ωm.

## Data (public, pinned)

| Data | Source | Commit |
|---|---|---|
| DESI DR2 BAO (13 values, full covariance) | github.com/CobayaSampler/bao_data, `desi_bao_dr2/` | bb0c1c9 (June 26, 2025) |
| Pantheon+ distances and STAT+SYS covariance | github.com/PantheonPlusSH0ES/DataRelease, `Pantheon+_Data/4_DISTANCES_AND_COVAR/` | c447f0f (Dec 21, 2022) |

## Supernova sample: what "1,580" means
- The Pantheon+ file has **1,701 light curves**. The data-release paper (Scolnic et al. 2022, arXiv:2112.03863) describes **1,550 distinct supernovae**; the cosmology paper is Brout et al. 2022 (arXiv:2202.04077). The file has 1,543 unique CID labels, a small difference in how names are counted.

| Step | Light curves | Distinct CIDs | Distinct objects (names merged within 3″ and 20 days of peak) |
|---|---|---|---|
| All rows in Pantheon+SH0ES.dat | 1,701 | 1,543 | 1,533 |
| Keep zHD > 0.01 | 1,590 | 1,473 | 1,465 |
| Also drop IS_CALIBRATOR == 1 (10 rows) | 1,580 | 1,466 | 1,458 |

**Open: the distinct count.** The papers report 1,550 distinct supernovae; this file gives 1,543 by CID and 1,533 by object. No CID is shared by two different supernovae. Ten supernovae carry two names each (for example 2016coj / SN2016coj). The 1,550 figure is not reproduced from this file. This does not affect the fit, which uses all light curves with the published covariance.
- This test keeps **zHD > 0.01** and excludes the SH0ES calibrators (IS_CALIBRATOR = 0). That leaves **1,580 light curves of 1,466 distinct supernovae**. Some supernovae were observed by more than one survey, and the covariance accounts for the duplicates.
- The usual Pantheon+-only cut (z > 0.01, calibrators kept) gives 1,590 light curves. The 10 calibrator light curves were dropped here because their distances come from Cepheids, not from the Hubble flow.

## Likelihood
- **Form:** Gaussian, χ² = rᵀC⁻¹r, using the full published covariance in both datasets.
- **BAO:** the overall scale c/(H₀ r_d) is fitted analytically.
- **Supernovae:** the absolute-magnitude/H₀ offset is fitted analytically.
- **Optimizer:** Nelder–Mead from 20 starts. A 161 × 400 grid over (A, θ₀) confirms the BAO minimum (grid minimum χ² = 1,252.4).

## Results (closure_results.json)

| Data | ΛCDM χ² (Ωm) | Sine best χ² | Δχ² |
|---|---|---|---|
| DESI DR2 BAO, 13 points | 10.27 (0.2975) | 1,251.0 | +1,240.8 |
| Pantheon+, 1,580 light curves | 1,386.4 (0.3326) | 1,420.2 | +33.8 |

The ΛCDM values agree with the published DESI DR2 (Ωm ≈ 0.297) and Pantheon+ (Ωm ≈ 0.334) results. That checks the pipeline.

## Reported result
The minimal sine closure fits the stated data substantially worse than flat ΛCDM. Independent replication is pending. If reproduced, this would disfavor this specific waveform, not all breathing-universe models.

## Limitations, and why the label is provisional
- Only best-fit χ² values are reported: no posterior, no parameter uncertainties, and no joint BAO + supernova fit.
- The sine best fit lies on the edge of its parameter space (A → ∞, θ₀ → π), so its χ² is a limiting value.
- Radiation and orbital details are ignored. These are negligible at z < 2.5.
- The test has not yet been checked by anyone independent. It moves from provisional to Tested when an outside reviewer reruns it.

## Run it
```
git clone https://github.com/CobayaSampler/bao_data.git
git clone --filter=blob:none --sparse https://github.com/PantheonPlusSH0ES/DataRelease.git pp
cd pp && git sparse-checkout set "Pantheon+_Data/4_DISTANCES_AND_COVAR" && cd ..
python3 closure_test.py .
```
Requires numpy, scipy and pandas. Runtime is a few minutes.
