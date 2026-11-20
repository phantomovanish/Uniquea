#!/usr/bin/env python3
"""Offline unit tests for controller masks; run python3 -m unittest -v test_ablation."""
import unittest
from dataclasses import replace

import numpy as np

import ablation as a
import model as m


class AblationTests(unittest.TestCase):
    def setUp(self):
        self.z = m.initial_state()
        self.z[:8] = [.18, .1, -.11, .14, .2, .3, .24, -.01]
        self.z[8:16] = [.08, -.11, .12, -.04, .14, -.07, .05, .02]
        self.z[16:] = np.arange(9)/7

    def test_full_on_rhs_equals_original(self):
        for option in ('A', 'B'):
            for time in (0., .7, 3.1, 14.8):
                with self.subTest(option=option, time=time):
                    actual = a.masked_rhs(time, self.z, option, m.P, a.MASKS['all_on'])
                    expected = m.rhs(time, self.z, option, m.P, enabled=True)
                    np.testing.assert_array_equal(actual, expected)

    def test_all_off_rhs_equals_original_and_no_work(self):
        for option in ('A', 'B'):
            actual = a.masked_rhs(.7, self.z, option, m.P, a.MASKS['all_off'])
            expected = m.rhs(.7, self.z, option, m.P, enabled=False)
            np.testing.assert_array_equal(actual, expected)
            np.testing.assert_array_equal(actual[m.POS], np.zeros(4))
            np.testing.assert_array_equal(actual[m.NEG], np.zeros(4))

    def test_individual_masked_forces_and_power_at_several_times(self):
        for option in ('A', 'B'):
            for name, mask in a.MASKS.items():
                for time in (.1, 1.3, 7.):
                    with self.subTest(option=option, mask=name, time=time):
                        active = np.asarray(mask)
                        off = m.rhs(time, self.z, option, m.P, enabled=False)
                        actual = a.masked_rhs(time, self.z, option, m.P, mask)
                        controller = m.controller(self.z[:8], self.z[8:16], time, m.P)
                        expected_force = np.where(active, controller, 0.)
                        acceleration_difference = actual[12:16] - off[12:16]
                        np.testing.assert_allclose(acceleration_difference,
                                                   expected_force / np.asarray(m.P.masses)[4:8],
                                                   rtol=2e-15, atol=1e-15)
                        power = expected_force*self.z[12:16]
                        np.testing.assert_allclose(actual[m.POS], np.maximum(power, 0.))
                        np.testing.assert_allclose(actual[m.NEG], np.maximum(-power, 0.))
                        self.assertEqual(actual[m.LOSS], off[m.LOSS])
                        np.testing.assert_array_equal(actual[8:12], off[8:12])

    def test_energy_identity_for_every_mask(self):
        step = 1e-6
        for option in ('A', 'B'):
            for name, mask in a.MASKS.items():
                with self.subTest(option=option, mask=name):
                    f = a.masked_rhs(.7, self.z, option, m.P, mask)
                    def H(z):
                        return m.energy(z[:8], z[8:16], option, m.P)
                    derivative = (H(self.z+step*f)-H(self.z-step*f))/(2*step)
                    booked = sum(f[m.POS])-sum(f[m.NEG])-f[m.LOSS]
                    self.assertAlmostEqual(derivative, booked, delta=3e-10)

    def test_full_on_trajectory_equals_original(self):
        p = replace(m.P, cycles=1)
        for option in ('A', 'B'):
            for perturb in (False, True):
                with self.subTest(option=option, perturb=perturb):
                    t, z = a.run(option, a.MASKS['all_on'], perturb, p)
                    expected_t, expected_z = m.run(option, p, perturb, enabled=True)
                    np.testing.assert_array_equal(t, expected_t)
                    np.testing.assert_array_equal(z, expected_z)

    def test_all_off_trajectory_has_no_work(self):
        p = replace(m.P, cycles=1)
        for option in ('A', 'B'):
            for perturb in (False, True):
                with self.subTest(option=option, perturb=perturb):
                    t, z = a.run(option, a.MASKS['all_off'], perturb, p)
                    np.testing.assert_array_equal(z[:, m.POS], np.zeros((len(t), 4)))
                    np.testing.assert_array_equal(z[:, m.NEG], np.zeros((len(t), 4)))
                    H, residual = m.ledger(z, option, p)
                    self.assertGreater(z[-1, m.LOSS], 0.)
                    self.assertLess(np.max(np.diff(H)), 1e-8)
                    self.assertLess(np.max(np.abs(residual)), 1e-7)

    def test_initial_conditions_and_target_signals_unchanged(self):
        self.assertEqual(a.P.cycles, 10)
        self.assertEqual(replace(a.P, cycles=m.P.cycles), m.P)
        base, perturb = m.initial_state(False), m.initial_state(True)
        np.testing.assert_allclose(perturb[:8]-base[:8], [.04, .04, 0, 0, 0, 0, 0, 0], atol=1e-16)
        np.testing.assert_array_equal(perturb[8:], base[8:])
        for time in (0., 1., 20.):
            np.testing.assert_array_equal(m.targets(time, a.P), m.targets(time, m.P))

    def test_reject_malformed_mask(self):
        for bad in ([True, True], [1, 1, 1, 1], [True]*5, 'all_on'):
            with self.subTest(mask=bad), self.assertRaises(ValueError):
                a.checked_mask(bad)


if __name__ == '__main__':
    unittest.main(verbosity=2)
