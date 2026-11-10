# Open target: derive the 2√2 cap from the structure of the plenum

**Status: the quantum 2√2 cap is Established. The plenum derivation is Open.**
*Correction (November 6, 2026, outside review accepted):* an earlier version of this page labeled the result
"Derived, conditional on one postulate about the plenum." That postulate, correlations as inner products of unit
vectors, is the quantum correlation structure itself, so it imports the answer. The earlier label is Superseded
and kept below for the record (House Rule 11).

## The established cap, re-derived (Established)
Write the CHSH operator as

```
Ŝ = A₀ ⊗ (B₀ + B₁) + A₁ ⊗ (B₀ − B₁)
```

with binary local outcomes, A_i² = B_j² = I, and opposite-wing measurements commuting, [A_i, B_j] = 0. Then

```
Ŝ² = 4I − [A₀, A₁] ⊗ [B₀, B₁]
```

(the squared terms give 2I + {B₀,B₁} + 2I − {B₀,B₁} = 4I; the cross terms give −[A₀,A₁] ⊗ [B₀,B₁]).
Each commutator has operator norm at most 2, so ‖Ŝ²‖ ≤ 8 and |⟨Ŝ⟩| ≤ ‖Ŝ‖ ≤ 2√2 (Tsirelson, 1980).
That is the standard quantum bound, not yet a plenum derivation.

## The plenum bridge (Open)
**Bell-cap route through the plenum — Open.** If the plenum is a local classical field with measurement-independent settings and Bell-local outcomes, it obeys the CHSH bound (\|S\| ≤ 2), regardless of curvature; that model class cannot reproduce the observed quantum violations. A quantum plenum could recover the 2√2 Tsirelson cap if it derives the appropriate dichotomic observables, spacelike-separated observable commutation, entangled state, and measurement correlations. A two-state detector and photon polarization's doubled angle are useful components or analogies, not complete derivations. GW propagation near c does not establish plenum microcausality. A new One Math result requires plenum dynamics fixed independently of Bell data and a distinct testable prediction.

What the plenum model must derive, not assume:

- **Binary local readouts:** how each analyzer produces +1 or −1. A bistable "cup with two minima" is a helpful
  measurement analogy but does not by itself derive the Born probabilities or the Bell correlations.
- **The effective state:** the plenum's equations must yield the joint (entangled) state or correlation rule.
- **Measurement structure:** why local observables have the required algebra, and why observables localized in
  spacelike-separated regions commute (as in a microcausal relativistic quantum field theory). A measured
  gravitational-wave speed near c does not establish that a proposed plenum has this property.
- **The correlation limit:** why derived correlations cannot exceed 2√2, while reaching it for the right state
  and settings. For photons the established correlation is E(a, b) = −cos 2(a − b); the doubled angle is the
  right polarization geometry (linear polarization is unchanged by a 180° turn), but inserting that correlation
  would reuse quantum mechanics, not derive it.

A curved metric by itself does not supply this bridge. The bridge must not reintroduce the retired local
shared-seed explanation or a propagating field signal between the wings.

**Candidate route:** specify the plenum's dynamics independently of Bell data, derive the measurement algebra and
correlations, recover the bound without importing it, and state a distinct testable prediction.

---

## Earlier version (Superseded, kept for the record)

### The postulate (P): the plenum is a shared geometric medium
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

### The derivation from P (inequality steps correct; P imports the quantum structure)
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

### Numerical check (check.py)
| Case | Largest \|S\| found |
|---|---|
| Kernel cos θ (the singlet, up to sign) | 2.828427 (= 2√2) |
| 25 random positive-definite isotropic kernels (P₀ to P₅, weights ≥ 0) | 2.433771 |
| Local deterministic strategies (all 16) | 2 |

### What it means for One Math (holds inside quantum theory)
The archive's open correction Δ_OM to −cos θ (Master §3.11) is now constrained: it must keep E a
positive-definite kernel, E(θ) = −Σ c_n P_n(cos θ) with c_n ≥ 0.
- No correction of that kind can push S above 2√2.
- Any weight moved off the lowest harmonic **lowers** the largest S and changes the *shape* of E(θ).
  Measured correlation curves can test that shape directly.

### Still open
- Derive postulate P from the plenum's own dynamics: why the medium's state space carries an inner product.
- Why the lowest harmonic: P₁ for spin-½, and cos 2Δ for photons.

## Credit (House Rule 14)
Clauser, Horne, Shimony and Holt (CHSH, 1969); microcausality in relativistic quantum field theory (Haag and Kastler, *J. Math. Phys.* 5, 1964); B. S. Tsirelson (Cirel'son, 1980), who proved the 2√2 bound
and that quantum ±1 correlations are exactly inner products of unit vectors; L. J. Landau (1987);
I. J. Schoenberg (1942), positive-definite functions on spheres; Popescu and Rohrlich (1994); Pawłowski et al.
(2009), who recovered 2√2 from "information causality." The corrected statuses and the bridge requirements follow an outside review (November 6, 2026).

## The cap as a dial (sharpened November 10, 2026; Established, re-checked)
With ±1 readouts and opposite wings commuting, Ŝ² = 4I − [A₀,A₁] ⊗ [B₀,B₁] (Landau, Phys. Lett. A 120 (1987) 54).
By Jordan's lemma (Jordan, 1875) any pair of ±1 observables splits into 2 × 2 blocks, so the largest CHSH value
over all states is exactly

```
S_max = √(4 + ‖[A₀,A₁]‖ · ‖[B₀,B₁]‖) = 2 √(1 + sin α · sin β),     sin α = ‖[A₀,A₁]‖ / 2
```

- α = 0 (wing A's two settings commute, as for any classical readout): cap 2, the CHSH bound.
- α = β = 90° on the Bloch sphere (anticommuting settings; 45° between polarizers): cap 2√2, Tsirelson's bound.
- Check (angle_check.py): 800 random sets of ±1 observables in dimensions 2 to 6 agree with the formula to
  3 × 10⁻¹⁵; the qubit dial at 0°, 30°, 45°, 60°, 90° gives 2, 2.2361, 2.4495, 2.6458, 2.8284.

**What it changes for the plenum (Open).** The target reduces to one local question per wing: do the plenum's
dynamics give two non-commuting binary readouts, with which α, while readouts at spacelike separation commute?
2√2 needs α = β = 90° and the matching entangled state. A Δ_OM correction would show as α < 90°, a cap below 2√2
(Candidate; testable against measured CHSH values). Deriving α still needs plenum dynamics fixed independently of
Bell data.
