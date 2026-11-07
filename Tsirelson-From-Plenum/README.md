# Open target: derive the 2√2 cap from the structure of the plenum

**Status (November 6, 2026): Derived, conditional on one postulate about the plenum (Postulated).**
The result reproduces quantum mechanics; it does not go beyond it. It also puts a testable limit on any
One Math correction to the correlation curve.

## The postulate (P): the plenum is a shared geometric medium
The plenum's shared state gives each measurement direction a unit vector in a common inner-product space,
and the correlation between the two ±1 outcomes is the inner product:

```
E(a, b) = ⟨u_a , v_b⟩ ,   |u_a| = |v_b| = 1
```

For an *isotropic* plenum (no preferred direction), E depends only on the angle θ between the settings and
must be a positive-definite function on the sphere. Schoenberg (1942) showed every such function is

```
E(θ) = Σ c_n P_n(cos θ) ,   c_n ≥ 0 ,   Σ c_n = 1        (P_n = Legendre polynomials)
```

The singlet's −cos θ is the lowest harmonic, P₁, up to sign. For photons, θ is replaced by 2Δ (the
Poincaré-sphere angle is twice the lab angle).

## The derivation (Derived from P; three lines)
With S = E(a,b) + E(a,b′) + E(a′,b) − E(a′,b′):

```
S = ⟨u_a , v_b + v_b′⟩ + ⟨u_a′ , v_b − v_b′⟩
  ≤ |v_b + v_b′| + |v_b − v_b′|                         (Cauchy–Schwarz, |u| = 1)
  ≤ √( 2 ( |v_b + v_b′|² + |v_b − v_b′|² ) ) = √(2·4) = 2√2   (parallelogram law: the sum is 4)
```

Equality needs v_b ⊥ v_b′ and each u along the matching sum or difference: the familiar 0°, 45°, 22.5°, 67.5°
photon settings.

**Why the local limit is 2:** a local hidden-variable model fixes each outcome as +1 or −1 in advance. That is
the same algebra with "vectors" squeezed onto one line (plain ±1 numbers), and the best of the 16 strategies
gives |S| = 2. **Why nothing reaches 4:** the Popescu–Rohrlich box (S = 4) respects no-signalling but has no
inner-product representation. The geometry of the plenum forbids it.

So the three numbers are three geometries: **2** for a line, **2√2** for a space with angles, **4** for no
geometry at all.

## Numerical check (check.py)
| Case | Largest \|S\| found |
|---|---|
| Kernel cos θ (the singlet, up to sign) | 2.828427 (= 2√2) |
| 25 random positive-definite isotropic kernels (P₀ to P₅, weights ≥ 0) | 2.433771 |
| Local deterministic strategies (all 16) | 2 |

## What it means for One Math (Candidate, testable)
The archive's open correction Δ_OM to −cos θ (Master §3.11) is now constrained: it must keep E a
positive-definite kernel, E(θ) = −Σ c_n P_n(cos θ) with c_n ≥ 0.
- No correction of that kind can push S above 2√2.
- Any weight moved off the lowest harmonic **lowers** the largest S and changes the *shape* of E(θ).
  Measured correlation curves can test that shape directly.

## Still open
- Derive postulate P from the plenum's own dynamics: why the medium's state space carries an inner product.
- Why the lowest harmonic: P₁ for spin-½, and cos 2Δ for photons.

## Credit (House Rule 14)
Clauser, Horne, Shimony and Holt (CHSH, 1969); B. S. Tsirelson (Cirel'son, 1980), who proved the 2√2 bound
and that quantum ±1 correlations are exactly inner products of unit vectors; L. J. Landau (1987);
I. J. Schoenberg (1942), positive-definite functions on spheres; Popescu and Rohrlich (1994); Pawłowski et al.
(2009), who recovered 2√2 from "information causality." This page recasts Tsirelson's theorem in the
archive's language. It is the plenum's version of a known result, not a new theorem.
