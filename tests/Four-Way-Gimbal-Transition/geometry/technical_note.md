# Balancing a four-scope, two-axis rocking assembly

**Label:** Manus (version not recorded) Uniquea Field Observations UFO

**Status: synthetic example, not a design validation.** This note tests the physics of three upper scope assemblies with loads and a fourth, smaller moving scope below them. The masses, dimensions, trim linkage, bearing geometry, and motion below are invented for illustration. They are not measurements of the owner's assembly and do not establish that a real mechanism will fit, be safe, or operate as described.

## Main conclusion

A pivot on the vertical plumb line through the combined center of mass (COM) has **zero total gravity torque at that particular pose**, whether it is below or above the COM. That does not establish stability:

- With the pivot **below** the COM, the vertically upright pose is an inverted pendulum: a small tilt makes gravity increase the tilt. It needs active stabilization or a mechanical constraint to stay upright.
- With the pivot **above** the COM, a small tilt produces a restoring gravity torque: it is passively stable in tilt. Without damping, it oscillates rather than settling. On a free third rotation about vertical, gravity provides no yaw restoring torque.
- With the pivot **at** the COM, the ideal rigid assembly is gravity-neutral at every attitude, not self-righting.

Thus, a point on the plumb line below the COM is only a *geometric candidate* for a pivot, not proof of a stable or practical gimbal. Actual pivots, travel clearances, load ratings, bearing friction, moving components, actuator torque, and controls must still be verified.

## Synthetic masses and positions

Coordinates are in meters in a body-fixed frame, with +z upward in the starting pose. Each upper scope's mass **includes its attached load**. The frame mass is represented at its own COM. The example treats each component as a point mass for inertia; the scopes' and frame's actual intrinsic inertia would add to the values calculated here.

| Component | Mass (kg) | Body-frame COM (x, y, z), m |
|---|---:|---|
| Upper scope 1 + load | 8 | (0.600, 0, 0.800) |
| Upper scope 2 + load | 8 | (-0.300, +0.3√3, 0.800) |
| Upper scope 3 + load | 8 | (-0.300, -0.3√3, 0.800) |
| Lower trim scope | 3 | (0.120, -0.080, -0.650) |
| Frame and hardware | 12 | (0, 0, 0) |

For component masses `m_i` at component COM positions `r_i`,

```text
M = Σ m_i
r_COM = (Σ m_i r_i) / M
```

Here `M = 39 kg` and `r_COM = (0.009231, -0.006154, 0.442308) m`. For a real scope or attached load with its own local COM `c_i`, first transform it into the common frame: `r_i = t_i + R_i c_i`. Include every load, mount, cable contribution if materially relevant, and moving scope; re-evaluate positions and COM when their configuration changes. This COM formula gives exact net gravity torque for distributed bodies in a uniform gravity field, but component COMs alone do not determine rotational inertia.

## Torque and stability check

Let the main pivot `P` be fixed in space; let `R` map the body frame to the world frame; let world gravity be `g_vec = (0, 0, -9.81) m/s²`. For a locked, rigid configuration,

```text
r = R (r_COM - P)
τ_gravity = r × (M g_vec)
U = M g z_COM_world
```

At the starting pose, selecting `P.x = r_COM.x` and `P.y = r_COM.y` makes `r` parallel to gravity, so the **entire gravity torque is zero**. This line does not, by itself, specify a unique pivot height or a realizable mechanism. More generally, a constrained joint need only have zero torque component along its free axis; the bearing may react the other components.

Compare two illustrative points on that plumb line:

| Pivot | Pivot P (m) | Signed COM height above pivot h (m) | Gravity torque at +10° pitch (N·m) | Gravity torque at +90° pitch (N·m) | Starting-pose tilt behavior |
|---|---|---:|---:|---:|---|
| Below COM | (0.009231, -0.006154, -0.250000) | +0.692308 | +45.994 | +264.870 | Unstable, inverted |
| Above COM | (0.009231, -0.006154, +0.750000) | -0.307692 | -20.442 | -117.720 | Restoring, pendulum-like |

The signed pitch `θ` is rotation about the body/world y axis from the starting pose. For these initially vertical configurations,

```text
h = z_COM - z_P
a = M g h
τ_y(θ) = a sin(θ)
U(θ) - U(0) = a [cos(θ) - 1]
τ_y = -dU/dθ
```

Near zero, `τ_y ≈ a θ`. Positive `a` reinforces a positive tilt; negative `a` opposes it. With the simplified point-mass inertia about y, `I_y = Σ m_i[(x_i-x_P)² + (z_i-z_P)²]`, the below-COM example has `I_y = 32.050 kg·m²` and a small-angle unstable growth rate of `2.875 s⁻¹`. The above-COM example has `I_y = 17.050 kg·m²` and an undamped small-angle frequency of `2.628 rad/s` (period `2.391 s`). These rates assume motion is restricted to this one axis and exclude intrinsic component inertia, damping, bearing effects, and actuator dynamics; they must not be used as real hardware specifications.

**Holding torque versus motion.** With the frame stationary at any of the nonzero pitch angles in the table, a main-axis motor or brake would need to supply the opposite of the displayed gravity torque, in addition to overcoming other loads. A moving trajectory also requires inertia, acceleration, velocity-dependent effects, friction, safety margin, and the dynamics of moving trim masses. The inverted case can be actively stabilized with suitable sensing, control authority, power, and safety provisions; a brake or hard support can constrain it, but does not make it passively self-righting.

## What a small-moving trim scope can and cannot do

A small trim **angle** is not the same thing as a small main-frame angle. A motor, gearing, linkage, or dynamic reaction can move the frame through a large angle even when the trim scope moves only a little; there is no general kinematic rule equating the two. Conversely, small trim travel does **not** automatically cancel a large gravity torque. It must be checked with the actual mass and lever geometry.

As a specific test, put the trim scope's body-fixed swing joint at `(0.120, -0.080, -0.150) m`. Its 3 kg COM hangs `L = 0.500 m` below the joint in the nominal pose. Allow a trim rotation `δ` about body y limited to ±10°. Its COM displacement from nominal is

```text
Δx = -L sin(δ)
Δz = L [1 - cos(δ)]
τ_y(θ, δ) = g [(M h + m_trim Δz) sin(θ) + m_trim Δx cos(θ)]
```

At the trim limit, `|Δx| = 0.086824 m` and `Δz = 0.007596 m`. For this geometry, the maximum magnitude of the **near-upright static-balance frame angle** is only **0.552°** for the below-COM pivot, or **1.246°** for the above-COM pivot. Those numbers are balance limits for the stated linkage and masses, not mechanical frame-travel limits. The below-COM equilibria remain inverted/unstable if trim is fixed; the above-COM equilibria are restoring in tilt.

At `θ = +90°`, the horizontal trim displacement has no world-horizontal projection (`cos θ = 0`). The small vertical displacement changes the y-axis gravity torque by at most **+0.224 N·m** here, leaving **+265.094 N·m** for the below-COM pivot or **−117.496 N·m** for the above-COM pivot at either ±10° trim limit. A different linkage, larger travel, counterweight, spring, or external drive changes the answer. The calculation is not a universal impossibility proof; it is a no-go for *passive static cancellation by this particular small trim* at 90°. A moving trim mass can temporarily exchange angular momentum with the rest of the system, but a static hold requires the net external moment to be balanced; internal motion alone is not a substitute for a continuously sustainable external reaction or a genuinely balanced geometry.

## Four smooth phase stages are not a 90° physical tilt

As one unambiguous example, let a phase `φ = ωt` and command

```text
roll  α = 10° cos(φ)
pitch β = 10° sin(φ)
R = R_y(β) R_x(α)
```

These are smooth periodic commands with a 90° **time-phase difference**. The stages below are separated by 90° of phase, while the physical axes each move only ±10°. `ω` determines speed; it is not specified here. The table gives world-axis **gravity torque**, not the full moving-actuator torque.

| Phase φ | Roll α | Pitch β | Gravity torque (τx, τy), below-COM pivot, N·m | Gravity torque (τx, τy), above-COM pivot, N·m |
|---:|---:|---:|---|---|
| 0° | +10° | 0° | (+45.994, 0) | (−20.442, 0) |
| 90° | 0° | +10° | (0, +45.994) | (0, −20.442) |
| 180° | −10° | 0° | (−45.994, 0) | (+20.442, 0) |
| 270° | 0° | −10° | (0, −45.994) | (0, +20.442) |

For an actual 90° *frame* rotation, use the 90° torque calculations above, not the phase-stage table. Smooth mathematical commands do not ensure smooth hardware motion: actuators must cover the required trajectory and torque, with practical clearance to stops and allowances for backlash, flex, and inertia. The equations also assume a fixed main pivot; offset or translating joint centers require a full linkage and external-reaction analysis.

## Bilateral orb, double fulcrum, and gimbal limitations

- A pair of opposed bearings/fulcrums on a common line ordinarily defines **one revolute axis**, not two independently steerable axes. A second tilt axis requires a nested joint, suitable spherical/flexure arrangement, or other compliant or articulated design. The word “orb” alone does not establish the actual joint degrees of freedom.
- For the simple fixed-point model above, the intended rotational axes must intersect at the same physical pivot center and supports must permit the required motion. If they do not intersect, the center moves and the simple fixed-pivot model is inadequate; nonconcurrent linkages can still work, but need their own kinematic and force analysis. Conflicting rigid supports can bind or overload bearings.
- A two-axis gimbal can provide two controlled pointing/tilt degrees of freedom; it does not independently control all three body orientations. A three-axis gimbal adds rotation about a third axis. Gravity alone cannot provide restoring yaw torque about world vertical.
- With three serial axes in a conventional Z-Y-X Euler/gimbal arrangement, a ±90° middle-axis angle can align the first and third axes and lose an independent instantaneous control direction. Quaternions avoid the corresponding **software coordinate** singularity, but do not remove actual physical joint-axis alignment. A 90° time-phase increment is not this condition. Other layouts can have different singularities; for two-axis pointing, azimuth may become undefined when the pointing direction aligns with its azimuth axis even though the physical joint axes need not align.
- Check bearings, orbital or spherical-joint articulation, cable wrap, scope collisions, motor limits, and hard stops through the **entire** commanded trajectory. A stop cannot be passed; reaching and gently reversing at a limit is theoretically possible, but relies on adequate control and clearance rather than an impact. Do not infer full 90° travel from a static plumb-line drawing.

## Reproduction

From this directory run `python3 calculate_gimbal.py`; the script uses only the Python standard library and prints the component-COM calculation, torque/stability cases, trim result, four phase stages, and numerical cross-checks. `calculations.log` is the captured machine output. All values depend on the synthetic inputs and assumptions listed above.
