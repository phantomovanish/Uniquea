# Three upper scopes and fourth gimbal: result summary

**Label:** Manus (version not recorded) Uniquea Field Observations UFO

**Status:** Derived inside synthetic examples; Postulated / Open as physical interpretation. **UFO number: pending.** The numerical geometry and the dimensionless dynamics are separate illustrative examples, not a fitted design or observational data.

## 1. Candidate fourth point

For a synthetic 39 kg assembly of three upper scope-and-load pairs, a lower trim scope and a frame, the combined center of mass is `(0.009231, -0.006154, 0.442308)` m in the specified body frame. A pivot vertically aligned with it has zero gravity torque at the upright pose, but this condition does not select a unique height or guarantee stability. In the example a pivot **below** the center of mass is unstable to small tilts; a pivot **above** it gives a restoring gravity torque. Active control or a constraint can stabilize a below-COM system.

For the chosen 3 kg trim mass and ±10° trim swing, the greatest near-upright static balance angle is approximately 0.552° with the below-COM pivot and 1.246° with the above-COM pivot. These are limits of that particular hypothetical trim geometry, not universal rotation limits. Geometry, gravity torque and the actual assumptions are in [`geometry/technical_note.md`](geometry/technical_note.md).

## 2. Axis and 90° check

A two-axis gimbal and a three-axis gimbal provide different degrees of freedom. For one conventional serial Z-Y-X arrangement, the angular-rate Jacobian in `axes_radix_checks.py` has rank three at 0°, 45° and 89° middle-axis pitch, and rank two at 90°. Its condition number increases to about 114.6 at 89°. This illustrates a configuration-specific gimbal lock: a quaternion can avoid a numerical coordinate singularity, but cannot undo alignment of physical joint axes. A different real mechanism must be assessed separately.

## 3. Rocking, transition and alternatives

The separate 16-state mechanical surrogate includes one rocker angle, one orb coordinate, two dynamic fulcra, three upper scope coordinates and one lower trim coordinate, with a velocity for each. The scope targets are actively driven with 0°, 90°, 180° and 270° **time phases**. They are not verified spatial field directions, and these targets do not establish a 90° physical turn.

In the final sampled cycle, the linear restoring law has a rocker range of approximately 0.14° to 10.54°; the smoothly nonlinear law ranges from −0.82° to 9.66°. Both reverse smoothly in the example. Neither reaches a physical 90° tilt, and the linear case does not cross zero during its final cycle. The upper response half-ranges are about 0.27–0.29 dimensionless units, compared with about 0.032 for the lower control; this small trim response follows from the chosen targets and parameters, not a demonstrated machine.

The fork ledger preserves both linear and nonlinear candidates. It also distinguishes the earlier archive's damping/no-damping Option A/B, which is a different fork. A four-value orb-damping probe and four controller ablations are recorded in [`EXTENDED_TESTS.md`](EXTENDED_TESTS.md). In the ablations, removing a controller also removes its periodic drive: the resulting trade-offs do not prove that three scopes alone remove most wobble.

The driven-model energy ledger tracks mechanical work supplied and extracted by the controllers and passive losses to a bookkeeping bath. Maximum sampled ledger residuals in the main 20-cycle examples are about 5.60e-12 and 1.83e-9 model-energy units. The finite trajectories are externally driven; no free or self-sustaining energy source is demonstrated. Twenty-six unit tests passed, the geometry output reproduced, and the main 20-cycle metrics were reproduced in a clean directory. These tests check the implementation, not a cosmic or hardware claim.

## Notation and open questions

`axes_radix_checks.py` verifies exact round trips between decimal and a positional radix-64 representation for several integers. A change of numeral base alone does not change a number or the physical equations. If “base 64” denotes new operations or physical rules in this project, those must be specified and tested separately; no speed or measurement-precision comparison was made.

The design still needs the actual locations, masses and inertias; the relationship between frame and trim motion; the desired 90° range; the crossing and bridge topology; and what “fractal” and “Factal” mean quantitatively. The exact meaning of “380 barrier” has not been confirmed. No mapping to observed quantum or cosmological data is supplied.

**Publication safeguard:** The project’s House Rule 16 requires a private withheld-terms scan of all tracked text, PDFs and workbooks before every push. This report was prepared without access to that private scan or term list. Do not push until the actual required scan has completed successfully; do not publish the private list. Retain the personal HAT draft separately until its owner selects material for release.

**Not claimed:** a working 3D gimbal design, physical 90° tip/flip, verified field or graviton, a self-powering device, or empirical confirmation of a cosmological hypothesis.
