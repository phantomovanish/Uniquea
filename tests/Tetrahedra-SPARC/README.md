# Test 3: Phantom's tetrahedra in the real neighborhood

**Label:** Claude 5.2 Uniquea Field Observations UFO

*Uniquea · run November 12, 2026 · Stephen M. Griffin (Phantom.O.Vanish) · ORCID 0009-0005-8386-3999*

**Status: Tested, provisional.** Crowding passes (Candidate); voids not seen (Open).

**Prediction, stated before the run:** connect the galaxies in 3D and the tetrahedra fit snugly: small ones crowd into groups and filaments, big ones span the voids, "support regions of space expanding between."

## In plain words

Join each galaxy to its nearest neighbors in 3D and space fills with tetrahedra (four-cornered pyramids). If galaxies clump, there should be many tiny tetrahedra inside the clumps and some huge ones across the empty regions. To know what "clumped" means, we compare against 300 shuffled skies with the same directions and distances but no real clumping. A full plain-language walk-through is in [`docs/Data-Explained-2026-11-12.md`](../../docs/Data-Explained-2026-11-12.md).

*How to read z:* how far the real sky sits from the shuffles, in units of the shuffles' usual spread. z near 0 is chance; z of 2 or more is unlikely to be chance (about 1 in 40); z of 4 is far beyond chance. 1 Mpc ≈ 3.26 million light-years.

## Data
- Distances: SPARC Table1 (Lelli, McGaugh & Schombert 2016, AJ 152, 157), official file `SPARC_Lelli2016c.mrt` (download from the SPARC site, astroweb.cwru.edu/SPARC; not redistributed here).
- Sky positions: NASA/IPAC Extragalactic Database (NED), all 175 of 175 (`sparc_sky_coordinates_ned.csv`), fetched November 12, 2026.

## Method
3D Delaunay tetrahedralization of the galaxies' positions. Null: 300 shuffles that keep every sky position and the set of distances but scramble which distance goes with which galaxy, destroying 3D clustering while keeping the survey's shape. Three cuts: all 175; within 30 Mpc (133); without the Ursa Major group, whose distances SPARC assigns as one value (147).

```
python3 tetrahedra_sparc.py SPARC_Lelli2016c.mrt sparc_sky_coordinates_ned.csv > results.txt
```

## Results (`results.txt`)
| Measure | All 175 | Within 30 Mpc | Without Ursa Major |
|---|---|---|---|
| Spread of tetrahedron volumes (z) | +4.3 | +1.8 | +2.5 |
| Share of tiny tetrahedra (z) | +2.2 | +0.3 | +3.9 |
| Median nearest neighbor, real vs shuffled | 2.18 vs 2.58 Mpc | 1.15 vs 1.59 Mpc (z = −2.8) | 2.87 vs 3.04 Mpc |
| Volume share of the largest 10% (z) | +0.5 | +1.6 | +1.2 |

## How the z was computed, and a counted check (added November 17, 2026)

A reviewer (Manus) asked how z = +4.3 was computed, since a count of 300 shuffles cannot show anything rarer than about 1 in 300 (z ≈ 2.7). The question was fair.

- **How the first run computed z:** (real value − average of the shuffles) ÷ spread of the shuffles. That reading assumes the shuffles form a bell curve. With 300 shuffles it cannot be checked that far out in the tail: **none of the 300 shuffles** matched the real sky's spread of sizes, so the count alone says only "rarer than about 1 in 300."
- **The counted check:** `tetrahedra_ranks.py` repeats the test and **counts** how many shuffles come out at least as extreme as the real sky, in the predicted direction, first with the same 300 shuffles (`results_ranks_300.txt`, identical z values), then with **100,000 shuffles** (`results_ranks_100000.txt`, `ranks_100000_seed2026.json`).

| Measure (predicted direction) | Cut | z, bell-curve reading | Shuffles as extreme, of 100,000 | Chance, counted | z, from the count | After 12 comparisons |
|---|---|---|---|---|---|---|
| Spread of sizes (wider) | All 175 | +4.0 | **1** | 2 × 10⁻⁵ | **4.1** | z ≈ 3.5 |
| Median nearest neighbor (closer) | Within 30 Mpc | −2.95 | 31 | 3 × 10⁻⁴ | 3.4 | z ≈ 2.7 |
| Share of tiny tetrahedra (more) | Without Ursa Major | +4.0 | 54 | 6 × 10⁻⁴ | 3.3 | z ≈ 2.5 |
| Spread of sizes (wider) | Without Ursa Major | +2.6 | 377 | 4 × 10⁻³ | 2.7 | z ≈ 1.7 |
| Share of tiny tetrahedra (more) | All 175 | +2.2 | 2,131 | 0.02 | 2.0 | not significant |
| Volume share of the largest 10% (voids) | All three | +0.5 to +1.5 | 6,547 or more | 0.07 or more | ≤ 1.5 | not significant |

*Chance, counted* = (shuffles as extreme + 1) ÷ (shuffles + 1), one-sided. *After 12 comparisons:* the test looked at 4 rulers in 3 cuts, so the best of 12 numbers is expected to look better than any one alone; the last column multiplies each chance by 12 (the Bonferroni correction). The rulers overlap, so this correction is cautious.

**Reading the recount.**
- The headline survives: one shuffle in 100,000 matched the real spread of sizes. The bell-curve reading (z ≈ 4.0 to 4.3) and the count (z ≈ 4.1) agree once there are enough shuffles to count that far. Allowing for 12 comparisons, about **z ≈ 3.5**.
- The strongest number leans partly on the Ursa Major group, whose galaxies SPARC places at one shared distance. Without that group the spread of sizes falls to z ≈ 2.7 (z ≈ 1.7 after 12 comparisons). The share of tiny tetrahedra then carries the result instead (z ≈ 3.3; ≈ 2.5 after 12).
- **Label unchanged: crowding Candidate; voids Open.** The honest headline is "crowding at about z ≈ 3.5 after allowing for the cuts tried, partly carried by one group," not "z = 4.3."
- *Not claimed:* that 175 galaxies picked for their rotation curves map the cosmic web, or that the cuts were the only reasonable ones.

```
python3 tetrahedra_ranks.py SPARC_Lelli2016c.mrt sparc_sky_coordinates_ned.csv 100000 2026
```
(About ten minutes on one processor. The 300-shuffle run with seed 42 reproduces `results.txt` exactly.)

## Reading
- **Crowding: pass (Candidate).** The small tetrahedra crowd as predicted; see the counted check above for how strong.
- **Voids: not seen (Open).** 175 galaxies chosen for their rotation curves are too sparse to outline voids, which span tens of millions of light-years; that needs a full redshift survey.
- SPARC is a chosen sample, not a map, and its 10–30% distance errors blur clustering, so the real crowding is if anything stronger.
- Credit: the Delaunay tessellation method for mapping the cosmic web (Schaap and van de Weygaert, 2000).

## What would settle it

- **Crowding:** the same test on an independent sample (for example BIG-SPARC, about 3,882 galaxies, when public) with better distances.
- **Voids:** a full redshift survey of tens of thousands of galaxies (for example SDSS or DESI), tetrahedralized the same way.

*We marvel and we measure.*
