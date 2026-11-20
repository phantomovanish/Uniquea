"""Reproducibility and conservation checks for calculate_ninety.py."""
import math
from fractions import Fraction as F
import unittest

from scipy.integrate import quad
import calculate_ninety as calc


class PrescribedMotionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = calc.compute()
        cls.results = cls.data['results']
        cls.ledger = cls.data['independent_analytic_ledger']

    def test_exact_geometry_and_point_mass_inertia(self):
        self.assertEqual(calc.M, 39)
        self.assertEqual(calc.X_COM, F(3, 325))
        self.assertEqual(calc.Y_COM, F(-2, 325))
        self.assertEqual(calc.Z_COM, F(23, 52))
        self.assertEqual(calc.H, F(9, 13))
        self.assertEqual(calc.I_Y, F(1041621, 32500))
        self.assertEqual(calc.C, F(26487, 100))

    def test_angles_kinematics_and_nonzero_initial_velocity(self):
        r = self.results
        self.assertEqual((r['angle_min_deg'], r['angle_max_deg']), (-90., 90.))
        self.assertEqual(r['end_to_end_angular_span_deg'], 180.)
        self.assertAlmostEqual(calc.theta(math.pi / 2), math.pi / 2, places=15)
        self.assertAlmostEqual(calc.theta(3 * math.pi / 2), -math.pi / 2, places=15)
        self.assertAlmostEqual(r['period_s'], 2 * math.pi / calc.OMEGA, places=13)
        self.assertAlmostEqual(r['peak_abs_angular_speed_deg_per_s'], 13.5, places=13)
        self.assertAlmostEqual(r['peak_abs_angular_acceleration_deg_per_s2'], 2.025, places=13)
        self.assertGreater(calc.velocity(0), 0)
        self.assertAlmostEqual(calc.velocity(0), calc.A * calc.OMEGA, places=15)

    def test_equation_of_motion_and_energy_rate(self):
        for n in range(101):
            phi = 2 * math.pi * n / 100
            v = calc.velocity(phi)
            acc = calc.acceleration(phi)
            tau = calc.motor_torque(phi)
            grav = calc.gravity_torque(phi)
            self.assertAlmostEqual(calc.I * acc, tau + grav - calc.B * v, places=11)
            dE_dt = calc.I * v * acc - calc.C_FLOAT * math.sin(calc.theta(phi)) * v
            self.assertAlmostEqual(tau * v, dE_dt + calc.B * v * v, places=11)

    def test_signed_gravity_and_turning_point_motor_torque(self):
        phi = math.pi / 2
        self.assertAlmostEqual(calc.gravity_torque(phi), 264.87, places=12)
        self.assertAlmostEqual(calc.motor_torque(phi), -calc.C_FLOAT - calc.K, places=12)
        self.assertAlmostEqual(calc.gravity_torque(3 * math.pi / 2), -264.87, places=12)
        self.assertGreater(calc.gravity_torque(0.01), 0)

    def test_torque_derivative_and_global_extrema(self):
        for phi in (0.2, 0.7, 1.4, 2.2, 3.4, 5.7):
            eps = 1e-6
            numeric = (calc.motor_torque(phi + eps) - calc.motor_torque(phi - eps)) / (2 * eps)
            self.assertAlmostEqual(calc.torque_phase_derivative(phi), numeric, delta=5e-8)
        p = self.results['torque_min_phase_rad']
        self.assertLess(abs(calc.torque_phase_derivative(p)), 1e-10)
        self.assertAlmostEqual(self.results['motor_torque_peak_abs_Nm'], 266.12509300047674, places=9)
        self.assertAlmostEqual(self.results['motor_torque_min_Nm'], -self.results['motor_torque_max_Nm'], places=11)
        peak = self.results['motor_torque_peak_abs_Nm']
        sampled_peak = max(abs(calc.motor_torque(2 * math.pi * n / 20000)) for n in range(20001))
        self.assertLessEqual(sampled_peak, peak + 1e-10)
        self.assertLess(peak - sampled_peak, 1e-5)

    def test_power_signs_and_exact_zero(self):
        root = self.results['torque_zero_first_phase_rad']
        self.assertLess(abs(calc.motor_torque(root)), 1e-10)
        cuts = [0, root, math.pi / 2, math.pi + root, 3 * math.pi / 2, 2 * math.pi]
        for (lo, hi), sign in zip(zip(cuts[:-1], cuts[1:]), (1, -1, 1, -1, 1)):
            self.assertGreater(sign * calc.actuator_work_integrand((lo + hi) / 2), 0)

    def test_numerical_and_analytic_energy_ledger(self):
        r, a = self.results, self.ledger
        self.assertAlmostEqual(r['positive_actuator_work_J'], 534.4313619702378, places=9)
        self.assertAlmostEqual(r['negative_actuator_work_J_signed'], -528.6176850926815, places=9)
        self.assertGreater(r['positive_actuator_work_J'], 0)
        self.assertLess(r['negative_actuator_work_J_signed'], 0)
        self.assertAlmostEqual(r['viscous_loss_J'], 3 * math.pi**3 / 16, places=13)
        self.assertAlmostEqual(r['mechanical_energy_change_J'], 0, places=13)
        self.assertLess(abs(r['energy_balance_residual_J']), 1e-10)
        self.assertAlmostEqual(r['net_actuator_work_J'], r['viscous_loss_J'], places=10)
        for key in ('positive_quadrature_minus_analytic_J',
                    'negative_quadrature_minus_analytic_J',
                    'loss_quadrature_minus_analytic_J'):
            self.assertLess(abs(a[key]), 1e-10)
        whole_cycle_work, _ = quad(calc.actuator_work_integrand, 0, 2 * math.pi,
                                   epsabs=1e-9, epsrel=1e-12)
        self.assertAlmostEqual(whole_cycle_work, r['net_actuator_work_J'], places=10)


if __name__ == '__main__':
    unittest.main(verbosity=2)
