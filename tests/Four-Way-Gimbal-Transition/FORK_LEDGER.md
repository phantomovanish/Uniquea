# Four-scope transition: fork ledger addendum

**Label:** Manus (version not recorded) Uniquea Field Observations UFO

**Status:** Derived inside the stated numerical models; Postulated / Open as a physical interpretation. This addendum does not replace the archive's workbook or assign a UFO number. **UFO number: pending.**

## Fork 1: rocker restoring response

- **Option A — linear:** a rocker restoring torque proportional to its angle throughout the stated simulation.
- **Option B — nonlinear:** the same law at small angle, smoothly increasing in stiffness through the chosen transition interval. It is an independently specified candidate, not established physics.
- **Status:** Both run, with identical synthetic initial conditions, external targets, dampers, loads and all other parameters. See [`README.md`](README.md), `model.py` and `output/summary.json`.
- **What decides it:** An independently specified physical system, its dimensions and measured torque–angle law; data acquired across and outside the proposed transition, including uncertainty and comparison with alternatives.
- **Result:** Both rock in the finite, externally driven examples. In the reported final cycle their sampled angles are approximately 0.14° to 10.54° for A and −0.82° to 9.66° for B. Neither makes a 90° physical flip. The assumed nonlinear transition is sampled, but the model does not establish which response nature follows. A frequency probe shows that the relative net-work ordering changes with the chosen drive frequency.

## Relationship to the earlier A/B fork

The archive's earlier [`Orb-Adiabatic-Check`](../Orb-Adiabatic-Check/README.md) uses **Option A for no damping and Option B for damping with lost energy booked to an Alt Field**. That is a different question in a different model; the new labels above do not supersede it. In this new model passive losses are integrated into an ideal bookkeeping variable, not a physical Alt Field with its own evolution. The actuators also exchange work with an external source. See the damping probe for the distinction between varying orb damping and eliminating every dissipative channel.

## Fork 2: pivot placement

- **Below the center of mass:** zero gravity torque at the exactly upright pose in the synthetic example, but that pose is passively unstable under small tilts. Active control or constraints could stabilize it; this example does not size a real controller.
- **Above the center of mass:** zero gravity torque at the aligned pose and a restoring gravity torque for small tilt, assuming the locked assembly and fixed pivot described in the geometry note.
- **Status:** Both calculated; not tested on hardware. The synthetic geometry is independent of the dimensionless rocking simulation.
- **What decides it:** Actual three-scope and frame geometry, individual masses and inertias, allowed travel, joint layout, actuator limits, friction and a measured response to disturbances. See [`geometry/technical_note.md`](geometry/technical_note.md).

## Other open forks and safeguards

- The four drive phases are timing markers, not four verified spatial field directions, a criss-crossing field, two actual bridge channels or a demonstrated infinity-shaped path. The precise connections, gimbal axes and meaning of the owner's QL/CL labels remain Open.
- A 90° stage of phase is not necessarily a 90° mechanical rotation. The Euler-axis test illustrates gimbal lock for one conventional three-axis arrangement; other mechanisms require their own analysis.
- The 380,000-year CMB/recombination record in [`../../docs/The-Birth-Canal-at-380000-Years.md`](../../docs/The-Birth-Canal-at-380000-Years.md) is separate from this arbitrary model-time transition. Whether the owner's phrase “380 barrier” refers to this event is not established here; no time or physical scale has been fitted.
- A positional base change can encode the same number in a different notation; it does not alone change the equations or their experimental accuracy. A distinct base-64 *physical model* would require its new rules and predictions to be stated and compared with data. `axes_radix_checks.py` verifies exact integer round trips only; it does not benchmark speed or hardware.
- “Fractal” and “Factal” are not defined or tested by this experiment. Demonstrating self-similarity or a scale-dependent geometric law would require specified observables, a scale range, an estimator and suitable data.

**Not claimed:** a validated 3D gimbal, observed quantum–cosmic field, literal full 90° tilt, graviton, independent energy source, mathematical proof of a cosmic phase transition, or an observational SPARC fit. Preserve older archive versions and complete the project's private pre-push scan before any publication.

**Owner's one-time exception, 2026-11-20:** The owner subsequently authorized one specified public branch push without the private scan. The exception and continuing future safeguard are recorded in [PUBLICATION_NOTE.md](PUBLICATION_NOTE.md); the previous wording is retained above.
