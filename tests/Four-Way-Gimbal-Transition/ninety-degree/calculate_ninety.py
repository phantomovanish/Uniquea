#!/usr/bin/env python3
"""One-axis prescribed-motion inverse dynamics; synthetic inputs only.

Run: python3 calculate_ninety.py
Writes results.json next to this source and prints the same JSON to stdout.
Uses scipy.integrate.quad and scipy.optimize.brentq; does not run a controller.
"""

from fractions import Fraction as F
import json
import math
from pathlib import Path
import platform

import scipy
from scipy.integrate import quad
from scipy.optimize import brentq


# Synthetic component COMs from technical_note.md. The two upper y coordinates
# are +0.3*sqrt(3) and -0.3*sqrt(3), so their mass-weighted y terms cancel.
# x,z suffice for point-mass inertia about y. Intrinsic inertia is NOT included.
COMPONENTS_XZ = (
    (8, F(3, 5), F(4, 5)),
    (8, F(-3, 10), F(4, 5)),
    (8, F(-3, 10), F(4, 5)),
    (3, F(3, 25), F(-13, 20)),
    (12, F(0), F(0)),
)
M = sum(m for m, _, _ in COMPONENTS_XZ)
X_COM = sum(m * x for m, x, _ in COMPONENTS_XZ) / M
Y_COM = (3 * F(-2, 25)) / M  # upper y contributions cancel
Z_COM = sum(m * z for m, _, z in COMPONENTS_XZ) / M
P_Z = F(-1, 4)
H = Z_COM - P_Z
I_Y = sum(m * ((x - X_COM) ** 2 + (z - P_Z) ** 2)
          for m, x, z in COMPONENTS_XZ)
G = F(981, 100)
C = M * G * H
B = 5.0                         # N m s/rad
OMEGA = 0.15                    # imposed phase rate, s^-1
A = math.pi / 2.0               # radians; +/- 90 degrees
I = float(I_Y)
C_FLOAT = float(C)
K = I * A * OMEGA**2           # acceleration-torque coefficient
D = B * A * OMEGA              # velocity-torque coefficient


def theta(phi):
    return A * math.sin(phi)


def velocity(phi):
    return A * OMEGA * math.cos(phi)


def acceleration(phi):
    return -A * OMEGA**2 * math.sin(phi)


def gravity_torque(phi):
    return C_FLOAT * math.sin(theta(phi))


def motor_torque(phi):
    return I * acceleration(phi) - gravity_torque(phi) + B * velocity(phi)


def torque_phase_derivative(phi):
    return (-K * math.cos(phi)
            - C_FLOAT * A * math.cos(A * math.sin(phi)) * math.cos(phi)
            - D * math.sin(phi))


def mechanical_energy(phi):
    return 0.5 * I * velocity(phi)**2 + C_FLOAT * math.cos(theta(phi))


def actuator_work_integrand(phi):
    # dW = tau * theta_dot * dt; dt = dphi / omega.
    return motor_torque(phi) * velocity(phi) / OMEGA


def loss_integrand(phi):
    return B * velocity(phi)**2 / OMEGA


def integrate(fun, lo, hi):
    return quad(fun, lo, hi, epsabs=2e-12, epsrel=1e-13, limit=100)


def compute():
    # Exactly one tau=0 in (0, pi/2); the other is a half-cycle later.
    root = brentq(motor_torque, 0.0, math.pi / 2.0,
                  xtol=1e-14, rtol=1e-14)
    # Exactly one torque minimum in (pi/2, pi); symmetry gives the maximum.
    torque_min_phase = brentq(torque_phase_derivative, math.pi / 2.0, math.pi,
                              xtol=1e-14, rtol=1e-14)
    torque_max_phase = (torque_min_phase + math.pi) % (2.0 * math.pi)

    cuts = [0.0, root, math.pi / 2.0, math.pi + root,
            3.0 * math.pi / 2.0, 2.0 * math.pi]
    pieces = [integrate(actuator_work_integrand, lo, hi)
              for lo, hi in zip(cuts[:-1], cuts[1:])]
    # The power signs on these intervals are +, -, +, -, + respectively.
    work_positive = sum(pieces[i][0] for i in (0, 2, 4))
    work_negative = sum(pieces[i][0] for i in (1, 3))
    work_net = work_positive + work_negative
    loss, loss_error = integrate(loss_integrand, 0.0, 2.0 * math.pi)
    energy_delta = mechanical_energy(2.0 * math.pi) - mechanical_energy(0.0)
    residual = work_net - energy_delta - loss

    # Independent exact-antiderivative ledger: E(phi)=I*v^2/2+C*cos(theta),
    # int b*v^2 dt = b*A^2*omega*(phi/2+sin(2*phi)/4).
    # At phi=pi/2 (and 3*pi/2), E=0 in this potential convention.
    root_energy = mechanical_energy(root)
    loss_factor = B * A**2 * OMEGA
    phase_adjust = root + math.sin(2.0 * root) / 2.0
    analytic_positive = (2.0 * root_energy
                         + loss_factor * (math.pi / 2.0 + phase_adjust))
    analytic_negative = (-2.0 * root_energy
                         + loss_factor * (math.pi / 2.0 - phase_adjust))
    analytic_loss = loss_factor * math.pi

    return {
        "case": "synthetic fixed-axis prescribed +/-90-degree motion; no controller simulated",
        "source": "/home/ubuntu/uniquea_gimbal_publish/tests/Four-Way-Gimbal-Transition/geometry/technical_note.md",
        "dependencies": {"python": platform.python_version(), "scipy": scipy.__version__},
        "input": {
            "mass_kg": M,
            "gravity_m_per_s2": float(G),
            "com_m": [float(X_COM), float(Y_COM), float(Z_COM)],
            "pivot_m": [float(X_COM), float(Y_COM), float(P_Z)],
            "h_m": float(H),
            "h_exact": str(H),
            "point_mass_inertia_about_y_kg_m2": I,
            "inertia_exact": str(I_Y),
            "intrinsic_component_inertia_included": False,
            "mgh_Nm": C_FLOAT,
            "mgh_exact": str(C),
            "amplitude_rad": A,
            "phase_rate_per_s": OMEGA,
            "viscous_coefficient_Nm_s_per_rad": B,
            "initial_angle_rad": theta(0.0),
            "initial_angular_speed_rad_per_s": velocity(0.0),
        },
        "results": {
            "period_s": 2.0 * math.pi / OMEGA,
            "angle_min_rad": -A,
            "angle_max_rad": A,
            "angle_min_deg": -90.0,
            "angle_max_deg": 90.0,
            "end_to_end_angular_span_deg": 180.0,
            "peak_abs_angular_speed_rad_per_s": A * OMEGA,
            "peak_abs_angular_speed_deg_per_s": math.degrees(A * OMEGA),
            "peak_abs_angular_acceleration_rad_per_s2": A * OMEGA**2,
            "peak_abs_angular_acceleration_deg_per_s2": math.degrees(A * OMEGA**2),
            "motor_torque_min_Nm": motor_torque(torque_min_phase),
            "motor_torque_max_Nm": motor_torque(torque_max_phase),
            "motor_torque_peak_abs_Nm": abs(motor_torque(torque_min_phase)),
            "torque_min_phase_rad": torque_min_phase,
            "torque_max_phase_rad": torque_max_phase,
            "torque_min_angle_deg": math.degrees(theta(torque_min_phase)),
            "torque_zero_first_phase_rad": root,
            "positive_actuator_work_J": work_positive,
            "negative_actuator_work_J_signed": work_negative,
            "negative_actuator_work_magnitude_J": -work_negative,
            "net_actuator_work_J": work_net,
            "viscous_loss_J": loss,
            "mechanical_energy_change_J": energy_delta,
            "energy_balance_residual_J": residual,
            "energy_balance_residual_definition": "W_positive + W_negative - delta_mechanical_energy - viscous_loss",
            "quadrature_sum_error_estimates_J": sum(p[1] for p in pieces),
            "quadrature_loss_error_estimate_J": loss_error,
        },
        "independent_analytic_ledger": {
            "positive_work_J": analytic_positive,
            "negative_work_J_signed": analytic_negative,
            "viscous_loss_J": analytic_loss,
            "positive_quadrature_minus_analytic_J": work_positive - analytic_positive,
            "negative_quadrature_minus_analytic_J": work_negative - analytic_negative,
            "loss_quadrature_minus_analytic_J": loss - analytic_loss,
            "loss_exact_form": "b*A^2*omega*pi = 3*pi^3/16 J",
        },
    }


def main():
    result = compute()
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    path = Path(__file__).resolve().parent / "results.json"
    path.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
