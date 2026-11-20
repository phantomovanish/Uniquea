#!/usr/bin/env python3
"""Tests for the independent damping probe; no network or pytest required.

Run: OPENBLAS_NUM_THREADS=1 python3 -B -m unittest -v test_damping_probe
No files are written by this test module.  Short one-cycle test runs do not
replace the 20-cycle experiment's parameters or output.
"""
import json
import math
import unittest
from dataclasses import asdict, replace
from unittest.mock import patch

import numpy as np

import damping_transition_probe as probe
import model


class DampingParameterTests(unittest.TestCase):
    def test_explicit_grid_and_only_orb_damping_changes(self):
        self.assertEqual(probe.GAMMAS, (0.0, 0.04, 0.16, 0.40))
        baseline = asdict(model.P)
        for gamma in probe.GAMMAS:
            with self.subTest(gamma=gamma):
                p = probe.parameters_for_gamma(gamma)
                d = asdict(p)
                for key in baseline:
                    if key != "damping":
                        self.assertEqual(d[key], baseline[key], key)
                self.assertEqual(p.damping[:1] + p.damping[2:], model.P.damping[:1] + model.P.damping[2:])
                self.assertEqual(p.damping[1], gamma * p.masses[1])
                self.assertAlmostEqual(p.damping[1] / p.masses[1], gamma)
        self.assertEqual(probe.parameters_for_gamma(0.16), model.P)

    def test_gamma_is_rate_even_if_orb_mass_differs(self):
        masses = list(model.P.masses)
        masses[1] = 2.5
        base = replace(model.P, masses=tuple(masses))
        p = probe.parameters_for_gamma(0.04, base)
        self.assertAlmostEqual(p.damping[1], 0.10)
        self.assertEqual(p.masses, base.masses)
        self.assertEqual(p.damping[0], base.damping[0])

    def test_rejects_invalid_gamma_and_mass(self):
        for gamma in (-0.01, float("nan"), float("inf")):
            with self.subTest(gamma=gamma):
                with self.assertRaises(ValueError):
                    probe.parameters_for_gamma(gamma)
        bad = replace(model.P, masses=(2.0, 0.0) + model.P.masses[2:])
        with self.assertRaises(ValueError):
            probe.parameters_for_gamma(0.04, bad)

    def test_gamma_zero_keeps_other_losses_and_enabled_drive(self):
        p0 = probe.parameters_for_gamma(0)
        self.assertEqual(p0.damping[1], 0)
        self.assertTrue(all(d > 0 for d in p0.damping[:1] + p0.damping[2:]))
        self.assertTrue(all(k > 0 for k in p0.kd))
        self.assertTrue(any(a > 0 for a in p0.target_amplitude))
        self.assertGreater(p0.omega, 0)
        z = model.initial_state()
        z[8:16] = [0.1, 0.2, -0.1, 0.1, 0.2, -0.2, 0.1, 0.1]
        self.assertGreater(model.rhs(0.3, z, "A", p0)[model.LOSS], 0)
        self.assertGreater(np.linalg.norm(model.controller(z[:8], z[8:16], 0.3, p0)), 0)

    def test_rhs_damping_delta_is_orb_force_and_bath_only(self):
        p0 = probe.parameters_for_gamma(0)
        pg = probe.parameters_for_gamma(0.4)
        z = model.initial_state()
        z[8:16] = [0.1, -0.2, -0.1, 0.1, 0.2, -0.2, 0.1, 0.1]
        for option in probe.OPTIONS:
            delta = model.rhs(0.3, z, option, pg) - model.rhs(0.3, z, option, p0)
            expected = np.zeros(25)
            expected[9] = -0.4 * z[9]
            expected[model.LOSS] = pg.damping[1] * z[9] ** 2
            np.testing.assert_allclose(delta, expected, atol=1e-14, rtol=1e-14)


class DampingShortRunTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Reducing duration/sampling only in tests keeps the suite quick.
        cls.base = replace(model.P, cycles=1, samples_per_cycle=120)
        cls.cases = []
        original_run = model.run
        with patch.object(probe.model, "run", wraps=original_run) as spy:
            for gamma in probe.GAMMAS:
                p = probe.parameters_for_gamma(gamma, cls.base)
                for option in probe.OPTIONS:
                    cls.cases.append(probe.evaluate_case(option, p, convergence=False))
        cls.calls = spy.call_args_list

    def test_both_rockers_get_identical_paired_setup_and_solver(self):
        self.assertEqual(len(self.calls), 16)
        for i, gamma in enumerate(probe.GAMMAS):
            p = probe.parameters_for_gamma(gamma, self.base)
            for j, option in enumerate(probe.OPTIONS):
                baseline_call = self.calls[4*i + 2*j]
                perturbed_call = self.calls[4*i + 2*j + 1]
                self.assertEqual(baseline_call.args, (option, p))
                self.assertEqual(baseline_call.kwargs, {})
                self.assertEqual(perturbed_call.args, (option, p))
                self.assertEqual(perturbed_call.kwargs, {"perturb": True})
            ca, cb = self.cases[2*i:2*i+2]
            self.assertEqual(ca["parameters"], cb["parameters"])
            self.assertEqual([ca["local_rocker_option"], cb["local_rocker_option"]], ["A", "B"])
        a, b = model.initial_state(), model.initial_state(perturb=True)
        expected = np.zeros(25)
        expected[0] = expected[1] = 0.04
        np.testing.assert_allclose(b-a, expected, atol=1e-15, rtol=0)

    def test_gamma_drive_product_energy_and_perturbation(self):
        period = 2 * math.pi / self.base.omega
        for case in self.cases:
            with self.subTest(gamma=case["gamma"], option=case["local_rocker_option"]):
                self.assertAlmostEqual(case["gamma_times_T_drive"], case["gamma"] * period, places=14)
                e = case["energy_ledger"]
                self.assertGreater(e["passive_loss_to_bookkeeping_bath"], 0)
                self.assertGreaterEqual(e["total_positive_actuator_work"], 0)
                self.assertGreaterEqual(e["total_negative_actuator_work_magnitude"], 0)
                self.assertAlmostEqual(e["net_actuator_work"], e["total_positive_actuator_work"] - e["total_negative_actuator_work_magnitude"], places=13)
                self.assertLess(abs(e["delta_H"] - e["net_actuator_work"] + e["passive_loss_to_bookkeeping_bath"]), 1e-7)
                self.assertLess(e["max_abs_sampled_ledger_residual"], 1e-7)
                self.assertLess(e["perturbed_max_abs_sampled_ledger_residual"], 1e-7)
                self.assertAlmostEqual(e["perturbed_delta_H"] - e["perturbed_net_actuator_work"] + e["perturbed_passive_loss_to_bookkeeping_bath"], 0.0, delta=1e-7)
                self.assertLess(abs(e["sampled_total_minus_ode_integrated_passive_loss"]), 0.01)
                if case["gamma"] == 0:
                    self.assertEqual(e["orb_passive_loss_sampled_trapezoid"], 0)
                else:
                    self.assertGreater(e["orb_passive_loss_sampled_trapezoid"], 0)
                d = case["perturbation"]
                self.assertAlmostEqual(d["initial_physical_state_difference_norm"], math.hypot(.04, .04))
                self.assertEqual(set(d["difference_norm_at_cycle_ends"]), {"1"})
                self.assertGreaterEqual(d["maximum_sampled_difference_norm"], d["final_sampled_difference_norm"])
                r = case["final_cycle_sampled_response"]
                self.assertLessEqual(r["theta_min"], r["theta_max"])
                fractions = [r["fraction_below_or_equal_lower_threshold"], r["fraction_between_thresholds"], r["fraction_above_or_equal_upper_threshold"]]
                self.assertAlmostEqual(sum(fractions), 1.0)
                json.dumps(case, allow_nan=False)

    def test_tighter_solver_rerun_is_diagnostic_only(self):
        p = probe.parameters_for_gamma(0.16, self.base)
        original_run = model.run
        with patch.object(probe.model, "run", wraps=original_run) as spy:
            case = probe.evaluate_case("B", p)
        self.assertEqual(len(spy.call_args_list), 3)
        self.assertEqual(spy.call_args_list[2].args, ("B", p))
        self.assertEqual(spy.call_args_list[2].kwargs, probe.TIGHT_SOLVER)
        self.assertLess(case["numerical_comparison"]["max_abs_physical_state_difference"], 1e-6)
        self.assertLess(case["numerical_comparison"]["max_abs_gross_ledger_state_difference"], 1e-5)

    def test_invalid_law_is_rejected(self):
        with self.assertRaises(ValueError):
            probe.evaluate_case("archived A", self.base, convergence=False)


if __name__ == "__main__":
    unittest.main(verbosity=2)
