# Tonight's data, explained

*Uniquea · November 12, 2026 · Stephen M. Griffin (Phantom.O.Vanish) · ORCID 0009-0005-8386-3999*

Two pieces of data came in tonight: **Test 3**, Phantom's tetrahedra among real galaxies, and **π digit by digit**, the spreadsheet. Here is what each one measured, what it found, and what it does not say.

## Part 1. Test 3: Phantom's tetrahedra

### The question

Phantom's instinct: connect the galaxies with lines in three dimensions and they form tetrahedra (four-cornered pyramids with triangle faces) that fit together snugly. **Small tetrahedra** should crowd where galaxies gather in groups and strands; **large tetrahedra** should span the empty regions between, the "support regions of space expanding between." The prediction was written down before the run.

### The data

- **175 galaxies** from SPARC, the standard catalog for galaxy rotation (Lelli, McGaugh and Schombert, 2016).
- **Distance** of each galaxy: from SPARC's own table.
- **Position in the sky** of each galaxy: from NASA's NED database, fetched in Colab tonight. All 175 were found (SIMBAD in Strasbourg was down, so NED supplied them all).
- Sky position plus distance puts every galaxy at its true spot in 3D, within the distance errors.

### The method, in plain words

1. **Connect the dots.** A standard method (Delaunay tetrahedralization) joins neighboring galaxies into tetrahedra that fill the space with no gaps or overlaps. The 175 galaxies make about 940 tetrahedra.
2. **Make a fair comparison.** To know whether the real pattern is special, we need a pattern that is *not* clustered but otherwise looks the same. So the computer kept every galaxy's direction in the sky and kept the same list of distances, but **shuffled which distance goes with which galaxy**, 300 times. The shuffles look like the same survey but have no real 3D clumping.
3. **Measure both** with the same rulers and compare.

### The rulers

| Ruler | What it shows |
|---|---|
| Spread of tetrahedron sizes | Clumpy space has both very small and very large tetrahedra, so a wide spread |
| Share of tiny tetrahedra | How many tetrahedra are less than 1% of the typical size: crowding |
| Nearest neighbor | How far each galaxy is from its closest companion |
| Share held by the largest 10% | How much of the space the biggest tetrahedra take up: voids |

### How to read a z-score

A z-score says how far the real galaxies sit from the shuffles, counted in the shuffles' usual wobble. **z near 0:** no different from chance. **z of 2 or more:** unlikely to be chance (about 1 in 40 or rarer). **z of 4:** far beyond chance.

### The results

| Ruler | All 175 | Within 30 Mpc (133) | Without Ursa Major (147) |
|---|---|---|---|
| Spread of sizes | **z = +4.3** | z = +1.8 | **z = +2.5** |
| Share of tiny tetrahedra | **z = +2.2** | z = +0.3 | **z = +3.9** |
| Nearest neighbor (real vs shuffled) | 2.18 vs 2.58 Mpc | **1.15 vs 1.59 Mpc (z = −2.8)** | 2.87 vs 3.04 Mpc |
| Share held by the largest 10% | z = +0.5 | z = +1.6 | z = +1.2 |

*1 megaparsec (Mpc) = about 3.26 million light-years. "Within 30 Mpc" is the near neighborhood, about 100 million light-years. The Ursa Major group was taken out in one cut because SPARC gives all its galaxies one shared distance, which could fake a clump.*

### What it means

- **Crowding: passes (Candidate).** The real galaxies make far more tiny tetrahedra and a far wider spread of sizes than the shuffles, and in the near neighborhood their nearest companions are about 30% closer than chance. The result survives with Ursa Major removed. Phantom's small, snug tetrahedra are there.
- **Voids: not seen (Open).** The biggest tetrahedra hold only slightly more space than in the shuffles, not enough to count. This is a limit of the data: voids span tens of millions of light-years, and 175 galaxies, chosen for their rotation curves rather than as a map of the sky, are too few to outline them. A void test needs a full sky survey of tens of thousands of galaxies.
- **Cautions.** SPARC is a chosen sample. Its distances carry 10–30% errors, which blur clumps, so the true crowding is likely stronger than measured.
- **In the field.** Astronomers map the cosmic web this very way (the Delaunay Tessellation Field Estimator, Schaap and van de Weygaert, 2000). Phantom's instinct matches a working tool of cosmology.

## Part 2. π digit by digit

### The question

Phantom saw that the digits of π zigzag (1, 4, 1, 5, 9…) yet the approach to π looks smooth. How smooth, and why?

### The data

The spreadsheet *Pi Truncation Steps* cuts π after 0, 1, 2 … 15 decimals (3, 3.1, 3.14 …) and records, for each cut:

- **The error left:** how far the cut value is below the true π.
- **How much the error shrank** from the step before.
- **The running digit sum:** the decimals added up (3.14 gives 1 + 4 = 5).

### The results

| Measure | Value |
|---|---|
| Error shrink per decimal, smallest step | 2.7× (adding a 1, at 3.141) |
| Error shrink per decimal, largest step | 34.9× (adding a 9, at 3.14159) |
| Error shrink per decimal, average | **9.7×** (about one power of ten) |
| Decimal places, all 16 cuts added up | 120 |
| Digit sums, all 16 cuts added up | 526 |
| Digit sum of the 15-decimal cut | 77 |

### What it means

- **Wavy step by step, smooth over the run.** Each decimal shrinks the error by a different amount, small for a small digit and large for a large one, but on average by about ten times. On a scale of powers of ten, the error glides down at a nearly constant angle.
- **Truncated π never zigzags.** The digits zigzag, but every cut approaches π from below and gets closer each time. Only rounded π overshoots.
- **Flip the error scale and it matches the digit sum.** The flipped error climbs about 0.98 powers of ten per decimal; the digit sum climbs about 5.1 per decimal; the two correlate at **0.993**.
- **Why they match (HAT).** Both are counting the same ladder of decimals. As a check, 2,000 sets of *random* digits match just as well (typical 0.990). So the match is real, but its cause is shared, not a hidden signal in π.
- **One measure Lila, one Lisa (Metaphor).** The clarity line is Lisa: almost perfectly straight (0.998 with the decimal count), Linear as in L.I.S.A. The digit sum is Lila: the same climb with a wobble at every rung. Same ladder, two measures, rising together: seamless by design.

## Where to find it

- **GitHub** (phantomovanish/Uniquea, commit c2d3f1b): `tests/Tetrahedra-SPARC/` (Test 3: coordinates, script, results) and `tests/Pi-Precision/pi_digit_steps.csv`.
- **The spreadsheet:** *Pi Truncation Steps* (Claude Sheets).
- **The Master:** Part 11, "Test 3" and "π, digit by digit: one measure Lila, one Lisa."

*We marvel and we measure.*
