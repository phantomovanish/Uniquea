# Lila damping from an action: the Alt Field as a bath

*Uniquea · run November 10, 2026 · Stephen M. Griffin (Phantom.O.Vanish) · ORCID 0009-0005-8386-3999*

**Status: mechanism Derived inside the model; the value ζ = ε_L is a Candidate condition on the bath.**

## In plain words

In Uniquea a small wobble on the young universe fades away. Earlier runs simply *wrote in* a drag to make it fade. A drag written by hand is a guess; a physicist wants it to *follow* from a law that keeps every bit of energy. So here the wobble is joined to many quiet background modes, the Alt Field, acting like a **bath**: no friction is written anywhere, the total energy is conserved, and the wobble still fades because its energy spreads into the bath and never comes back. With a single partner instead of a bath, the energy sloshes back and forth and the wobble never fades. The fade rate comes out at Uniquea's purity value ε_L (0.43444%) **if** the bath has the right strength; that strength is a condition, not yet derived. *Terms:* ζ is the damping ratio (how fast the wobble fades per cycle); ε_L the purity value; ω₀ the wobble's frequency; J(ω) how strongly the bath responds at each frequency.

## The question
Earlier runs (Orb-Adiabatic-Check/optionAB.py) put the Lila damping into the wobble's equation by hand, as a drag
term −γp with ζ = ε_L. The Master listed as Open: derive that the Lila damping acts on the wobble.

## The construction (Established method)
Couple the wobble y to many Alt Field modes q_k in an energy-conserving Hamiltonian (Caldeira and Leggett, 1983;
Ford, Kac and Mazur, 1965):

```
H = p²/2 + V(y) + Σ_k [ P_k²/2 + (ω_k²/2)(q_k − c_k y/ω_k²)² ]
J(ω) = (π/2) Σ_k c_k²/ω_k δ(ω − ω_k)        (the bath's spectral density)
weak coupling:  ζ = J(ω₀) / (2ω₀²)          (model units, Ā₂ = M_χ = 1)
```

No friction is written anywhere. V(y) is the Lila-mapped cup of Master 3.15 at m = 0, canonical exact mapping
(lip exactly at 1 − ε_L).

## Results (results_exact.json; results_first_order.json for the November 5 mapping)

| Quantity | Value |
|---|---|
| ω₀ (exact cup) | 0.0932817; ω₀² = 2ε_L × 1.00145 |
| Lila condition J(ω₀) = 2ε_L ω₀² | 7.5606 × 10⁻⁵ (≈ 4ε_L²) |
| Measured ζ / ε_L, 1,200 modes to 3ω₀ | 1.0018 (same for 600 and 2,400 modes; 1.0007 with modes to 6ω₀) |
| Cycles per amplitude e-fold | 36.57 (Master: 36.6) |
| Wobble energy left after 110 cycles | 0.24% (amplitude ≈ 1/20) |
| Energy ledger, wobble + bath | closed to 3 × 10⁻⁷ of the wobble's energy (numerical tolerance) |
| One partner mode instead of a bath | energy falls to 2.1%, then returns to 99.9997%: no damping (Option A) |

## Reading
- A continuum of Alt Field modes turns the hand-written drag into a consequence of an action: the loss is one-way
  and booked exactly in the bath (Option B). A single partner gives energy back (Option A).
- ζ = ε_L is not derived: it holds when the bath's J(ω₀) = 2ε_L ω₀². It becomes a prediction only if the field
  that plays the bath (the Seed Bank b, the gauge triad's modes, or another) has a J(ω) fixed without using ε_L.
- Hubble term: in the covariant action the wobble also feels 3Hχ̇. It vanishes at the perch (H_b = 0) and pumps
  the wobble during contraction (H < 0); net damping needs 2ε_L ω₀ > 3|H|.

## What would settle it

A physical field for the bath (the Seed Bank, the gauge triad's modes, or another) whose J(ω) is fixed **without** using ε_L. If it gives ζ = ε_L, the Candidate becomes a prediction; if not, the condition is retired.

## Run
```
python3 lila_bath.py results_exact.json exact          # about 3 seconds
python3 lila_bath.py results_first_order.json first
python3 lila_bath.py out.json exact 2400               # bath-size check
```
Credit: the bath method is Caldeira and Leggett's (1983), after Ford, Kac and Mazur (1965) and Ullersma (1966).
