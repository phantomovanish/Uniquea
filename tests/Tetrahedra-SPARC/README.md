# Test 3: Phantom's tetrahedra in the real neighborhood

*Uniquea · run November 12, 2026 · Stephen M. Griffin (Phantom.O.Vanish) · ORCID 0009-0005-8386-3999*

**Status: Tested, provisional.** Crowding passes (Candidate); voids not seen (Open).

**Prediction, stated before the run:** connect the galaxies in 3D and the tetrahedra fit snugly: small ones crowd into groups and filaments, big ones span the voids, "support regions of space expanding between."

## Data
- Distances: SPARC Table1 (Lelli, McGaugh & Schombert 2016, AJ 152, 157), official file `SPARC_Lelli2016c.mrt`.
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

## Reading
- **Crowding: pass (Candidate).** The small tetrahedra crowd as predicted.
- **Voids: not seen (Open).** 175 galaxies chosen for their rotation curves are too sparse to outline voids, which span tens of millions of light-years; that needs a full redshift survey.
- SPARC is a chosen sample, not a map, and its 10–30% distance errors blur clustering, so the real crowding is if anything stronger.
- Credit: the Delaunay tessellation method for mapping the cosmic web (Schaap and van de Weygaert, 2000).

*We marvel and we measure.*
