# Test 3 Recount: The Review, the Rerun and the Results

**Label:** Claude 5.2 Uniquea Field Observations UFO

*Uniquea · November 17, 2026 · Stephen M. Griffin (Phantom.O.Vanish) · ORCID 0009-0005-8386-3999 · DOI 10.5281/zenodo.24212673 · UFO-008*

> **The finding, in one line:** Phantom's tetrahedra still show crowding among real galaxies, but at about **z ≈ 3.5** once the several cuts are allowed for, not the **z = +4.3** first reported. Crowding stays a **Candidate**; the voids stay **Open**.

## 1. The review

A review by Manus raised four points about a PDF summary made outside this archive. Each was checked against the record.

| # | Manus's point | Verdict | What was done |
|---|---|---|---|
| 1 | The PDF says "175 galaxies, 2MRS," but Test 3 used the 175 SPARC galaxies. 2MRS is a survey of about 44,000 galaxies | **Correct.** That PDF named the wrong catalog; this record never used 2MRS | Noted in the public Master |
| 2 | z = +4.3 needs explaining: if 1 of 300 shuffles matched, the count gives p ≈ 0.007, about z ≈ 2.5. Several cuts were tried, which weakens any single result | **Correct on both counts.** The z was a bell-curve reading, and 12 numbers were read | The test was rerun with 100,000 shuffles and counted (sections 2 to 4) |
| 3 | π: the first 50 digits sum to 247 and the running total passes 100 at the 21st digit; "p = 2.2 × 10⁻¹⁶" is the smallest number R prints; the tenfold shrink per digit is how decimals work | **Arithmetic checked and correct.** "Climbs to 100" and "2.2 × 10⁻¹⁶" are not in this record | A plain line added to the data explainer: a teaching demo, not evidence about nature |
| 4 | The Higgs field doesn't make galaxies glow; most ordinary mass comes from the strong force; the switch from "rain" to "interaction" holds up | **Already agreed.** Filed as "Phantom's Higgs rain: not supported; recast" | None needed |

## 2. What went wrong with z = +4.3

The first run (November 12) compared the real sky with **300** random skies, made by shuffling which distance goes with which galaxy. It reported

z = (real value − average of the shuffles) ÷ spread of the shuffles.

That formula assumes the shuffles fall in a bell curve. Far out in the tail the assumption cannot be checked with 300 shuffles: **none of the 300** came out as extreme as the real sky, so a count alone says only "rarer than about 1 in 300," which is **z ≈ 2.7 at most**. The reported 4.3 went beyond what the run could show.

## 3. The rerun

- **Same data, same method:** 175 SPARC galaxies (Lelli, McGaugh and Schombert 2016), sky positions from NASA/IPAC NED, 3D Delaunay tetrahedra, distances shuffled among the galaxies.
- **Counted, not assumed:** for each ruler, the number of shuffles at least as extreme as the real sky, in the predicted direction. Chance = (that number + 1) ÷ (shuffles + 1), which is never zero and cannot overstate.
- **100,000 shuffles** (seed 2026, about ten minutes). The first 300-shuffle run (seed 42) reproduces the original results exactly.
- **12 comparisons allowed for:** four rulers in three cuts. Each chance is multiplied by 12 (the Bonferroni correction); because the rulers overlap, this is cautious.

## 4. The results

| Measure | Cut | Shuffles as extreme, of 100,000 | Chance, counted | z, counted | z, after 12 comparisons | First run: z as reported / most a 300 count can show |
|---|---|---|---|---|---|---|
| **Spread of sizes** | **All 175** | **1** | 2.0 × 10<sup>-5</sup> | **4.11** | **3.49** | +4.27 / 2.71 |
| Share of tiny tetrahedra | All 175 | 2,131 | 0.021 | 2.03 | not significant | +2.15 / 1.93 |
| Median nearest neighbor | All 175 | 2,976 | 0.03 | 1.88 | not significant | -1.69 / 1.84 |
| Largest 10% share (voids) | All 175 | 31,594 | 0.316 | 0.48 | not significant | +0.45 / 0.30 |
| Spread of sizes | Within 30 Mpc | 4,103 | 0.041 | 1.74 | not significant | +1.80 / 1.84 |
| Share of tiny tetrahedra | Within 30 Mpc | 36,545 | 0.365 | 0.34 | not significant | +0.31 / 0.41 |
| Median nearest neighbor | Within 30 Mpc | 31 | 3.2 × 10<sup>-4</sup> | 3.41 | 2.67 | -2.83 / 2.71 |
| Largest 10% share (voids) | Within 30 Mpc | 6,547 | 0.065 | 1.51 | not significant | +1.61 / 1.61 |
| Spread of sizes | Without Ursa Major | 377 | 3.8 × 10<sup>-3</sup> | 2.67 | 1.69 | +2.51 / 2.71 |
| Share of tiny tetrahedra | Without Ursa Major | 54 | 5.5 × 10<sup>-4</sup> | 3.26 | 2.48 | +3.85 / 2.71 |
| Median nearest neighbor | Without Ursa Major | 26,339 | 0.263 | 0.63 | not significant | -0.58 / 0.51 |
| Largest 10% share (voids) | Without Ursa Major | 11,864 | 0.119 | 1.18 | not significant | +1.16 / 1.18 |

*z, counted:* the bell-curve z that has the same chance. *Not significant:* the chance after 12 comparisons is above 5%.

## 5. What the results say

- **The headline survives, smaller.** One shuffle in 100,000 matched the real spread of tetrahedron sizes: z ≈ 4.1 by direct count, and **z ≈ 3.5 after allowing for the 12 numbers read**.
- **Part of it is one group.** SPARC places the galaxies of the Ursa Major group at one shared distance, which can fake a clump. Without them the spread falls to z ≈ 2.7 (1.7 after 12), and the share of tiny tetrahedra carries the result instead (z ≈ 3.3; 2.5 after 12).
- **In the near neighborhood** (within 30 Mpc), galaxies' nearest companions are about 30% closer than chance (z ≈ 3.4; 2.7 after 12).
- **The voids are not seen** in any cut. 175 galaxies picked for their rotation curves are too few to outline regions tens of millions of light-years across.
- **Labels:** crowding **Candidate**; voids **Open**. *Not claimed:* that this sample maps the cosmic web, or that these cuts were the only reasonable ones.

## 6. What would settle it

- **Crowding:** the same test on an independent sample with better distances (for example BIG-SPARC, about 3,882 galaxies, when public).
- **Voids:** a full redshift survey of tens of thousands of galaxies, tetrahedralized the same way.

## 7. What changed in the record

- The Test 3 write-up gained "How the z was computed, and a counted check," with the script `tetrahedra_ranks.py` and its results.
- UFO-008 now reads z ≈ 4.1 counted, ≈ 3.5 after 12 comparisons, and notes that z = +4.3 was first reported from 300 shuffles.
- The data explainer carries the correction and a plain line on the π demo.
- The public Master answers all four review points, with credit to Manus.
- Spreadsheet with live formulas: *Test 3 Recount: 100,000 Shuffles* (Claude Sheets), and `Test3-Recount-100000-Shuffles.xlsx`.

*Rule:* outside review is welcome and is filed with its corrections, whoever gives it. The highest seal is a stranger's count.
