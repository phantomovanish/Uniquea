#!/usr/bin/env python3
"""Synthetic gimbal example. Run: python3 calculate_gimbal.py

Coordinates, component masses and trim linkage are illustrative, not measurements.
The model treats each component as a point mass and ignores motor/frame intrinsic
rotational inertia, friction, and flex. SI units throughout; angles printed in deg.
"""

import math

G = 9.81
DEG = math.pi / 180.0

# Fixed-frame body coordinates. Each upper mass includes its attached load.
# The three upper scopes form an equilateral horizontal triangle about the z axis.
# The frame's mass is represented by its COM, not by an inertia tensor.
BODIES = [
    ("upper 1 + load", 8.0, (0.6, 0.0, 0.8)),
    ("upper 2 + load", 8.0, (-0.3, 0.3 * math.sqrt(3.0), 0.8)),
    ("upper 3 + load", 8.0, (-0.3, -0.3 * math.sqrt(3.0), 0.8)),
    ("trim scope", 3.0, (0.12, -0.08, -0.65)),
    ("frame and hardware", 12.0, (0.0, 0.0, 0.0)),
]
TRIM_MASS = 3.0
TRIM_ARM = 0.5  # horizontal trim-swing pivot is 0.5 m above nominal trim COM
TRIM_LIMIT = 10.0 * DEG
FRAME_AMPLITUDE = 10.0 * DEG


def vadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(k, a):
    return tuple(k * x for x in a)


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def rotate_x(a, v):
    x, y, z = v
    c, s = math.cos(a), math.sin(a)
    return (x, c*y - s*z, s*y + c*z)


def rotate_y(a, v):
    x, y, z = v
    c, s = math.cos(a), math.sin(a)
    return (c*x + s*z, y, -s*x + c*z)


def total_torque(pivot, roll=0.0, pitch=0.0):
    # Active body-to-world orientation R = Ry(pitch) Rx(roll).
    tau = (0.0, 0.0, 0.0)
    for _, mass, pos in BODIES:
        r = rotate_y(pitch, rotate_x(roll, vsub(pos, pivot)))
        tau = vadd(tau, cross(r, (0.0, 0.0, -mass*G)))
    return tau


def potential(pivot, pitch):
    return sum(mass*G*rotate_y(pitch, vsub(pos, pivot))[2]
               for _, mass, pos in BODIES)


def inertia_y_point_mass(pivot):
    return sum(mass*((pos[0]-pivot[0])**2 + (pos[2]-pivot[2])**2)
               for _, mass, pos in BODIES)


def fmt_vec(v, digits=6):
    return '(' + ', '.join(f'{x:+.{digits}f}' for x in v) + ')'


M = sum(mass for _, mass, _ in BODIES)
COM = tuple(sum(mass*pos[j] for _, mass, pos in BODIES)/M for j in range(3))
CASES = [('below_COM', (COM[0], COM[1], -0.25)),
         ('above_COM', (COM[0], COM[1], 0.75))]

assert abs(M - 39.0) < 1e-12
assert abs(COM[0] - 0.36/39.0) < 1e-12
assert abs(COM[1] + 0.24/39.0) < 1e-12
assert abs(COM[2] - 17.25/39.0) < 1e-12

print('SYNTHETIC GIMBAL CALCULATION (SI units)')
print(f'M = {M:.6f} kg; g = {G:.6f} m/s^2')
print(f'COM = {fmt_vec(COM)} m')
print('Upper 2 and 3 use y = +/-0.3*sqrt(3) m.')
print()
for label, p in CASES:
    h = COM[2] - p[2]  # signed vertical COM position relative to pivot
    K = M*h            # kg m
    torque_scale = M*G*h
    Iy = inertia_y_point_mass(p)
    rate = math.sqrt(abs(torque_scale)/Iy)
    print(f'CASE {label}:')
    print(f'  Pivot = {fmt_vec(p)} m; h = {h:+.6f} m; M*h = {K:+.6f} kg m')
    print(f'  Zero-angle gravity torque = {fmt_vec(total_torque(p))} N m')
    print(f'  tau_y(theta) = ({torque_scale:+.6f}) * sin(theta) N m')
    print(f'  U(theta)-U(0) = ({torque_scale:+.6f}) * (cos(theta)-1) J')
    print(f'  Point-mass I_y = {Iy:.6f} kg m^2')
    if h > 0:
        print(f'  Linearized unstable exponential rate = {rate:.6f} s^-1')
    else:
        print(f'  Linearized undamped natural frequency = {rate:.6f} rad/s')
        print(f'  Point-mass small-angle period = {2*math.pi/rate:.6f} s')
    print(f'  At +10 deg pitch: gravity torque = {fmt_vec(total_torque(p, pitch=10*DEG))} N m')
    print(f'  At +90 deg pitch: gravity torque = {fmt_vec(total_torque(p, pitch=90*DEG))} N m')

    # Trim swing about its own horizontal body-y joint: delta z positive either way.
    # Delta r_trim = (-L sin(delta), 0, L (1 - cos(delta))).
    dx = -TRIM_ARM*math.sin(TRIM_LIMIT)
    dz = TRIM_ARM*(1-math.cos(TRIM_LIMIT))
    # For this example K is far from zero; the equilibrium branch near theta=0
    # is monotone in |delta| over the chosen range, so its maximum is at 10 deg.
    theta_near_upright = math.atan(TRIM_MASS*TRIM_ARM*math.sin(TRIM_LIMIT) /
                                   (K + TRIM_MASS*dz))
    extra_90 = TRIM_MASS*G*dz
    print(f'  Trim limit |delta| = 10 deg gives |delta_x| = {abs(dx):.6f} m, delta_z = {dz:.6f} m')
    print(f'  Largest |static-balanced frame pitch| near zero = {abs(theta_near_upright)/DEG:.6f} deg')
    print(f'  Trim torque change at +90 deg pitch (at |delta|=10 deg) = {extra_90:+.6f} N m')
    print(f'  Remaining gravity tau_y at +90 deg pitch = {torque_scale + extra_90:+.6f} N m')
    print()

print('FOUR TEMPORAL PHASE STAGES: roll=10*cos(phi) deg, pitch=10*sin(phi) deg')
print('Gravity torques are world-axis components; stationary holding torque is their negative.')
for phi_deg in (0, 90, 180, 270):
    phi = phi_deg*DEG
    roll, pitch = FRAME_AMPLITUDE*math.cos(phi), FRAME_AMPLITUDE*math.sin(phi)
    values = [total_torque(p, roll, pitch) for _, p in CASES]
    print(f'  phi={phi_deg:3d} deg: roll={roll/DEG:+.3f}, pitch={pitch/DEG:+.3f} deg; '
          f'below tau={fmt_vec(values[0], 4)}; above tau={fmt_vec(values[1], 4)} N m')

# Independent checks: sum-of-particles vs COM resultant and virtual work.
for _, p in CASES:
    for roll, pitch in ((0.0, 0.0), (0.07, -0.11), (0.0, math.pi/2)):
        r_com = rotate_y(pitch, rotate_x(roll, vsub(COM, p)))
        t_com = cross(r_com, (0.0, 0.0, -M*G))
        t_sum = total_torque(p, roll, pitch)
        assert max(abs(a-b) for a,b in zip(t_com, t_sum)) < 1e-10
    theta = 0.2
    step = 1e-6
    dU = (potential(p, theta+step) - potential(p, theta-step))/(2*step)
    assert abs(total_torque(p, pitch=theta)[1] + dU) < 1e-6
    h = COM[2] - p[2]
    for delta in (-TRIM_LIMIT, 0.0, TRIM_LIMIT):
        dx = -TRIM_ARM*math.sin(delta)
        dz = TRIM_ARM*(1-math.cos(delta))
        # Independently recompute all point-mass torques with the trim moved.
        def moved_torque(pitch):
            total = (0.0, 0.0, 0.0)
            for name, mass, pos in BODIES:
                if name == 'trim scope':
                    pos = vadd(pos, (dx, 0.0, dz))
                r = rotate_y(pitch, vsub(pos, p))
                total = vadd(total, cross(r, (0.0, 0.0, -mass*G)))
            return total[1]
        check_pitch = 0.2
        expression = G*((M*h + TRIM_MASS*dz)*math.sin(check_pitch) +
                         TRIM_MASS*dx*math.cos(check_pitch))
        assert abs(moved_torque(check_pitch) - expression) < 1e-10
        theta_balance = math.atan(-TRIM_MASS*dx/(M*h + TRIM_MASS*dz))
        assert abs(moved_torque(theta_balance)) < 1e-10
print('\nChecks passed: particle/COM torques agree; tau_y = -dU/dtheta;')
print('trim torque and the near-upright balance formula agree with particle sums.')
