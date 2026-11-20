#!/usr/bin/env python3
"""Independent numerical consistency tests; no network or pytest required."""
import math
import unittest
from dataclasses import replace

import numpy as np
import model as m


class ToyModelTests(unittest.TestCase):
    def test_options_identical_below_transition_and_odd(self):
        for theta in (-.10, -.07, 0., .07, .10):
            self.assertAlmostEqual(m.rocker_force(theta, 'A'), m.rocker_force(theta, 'B'), places=14)
            self.assertAlmostEqual(m.rocker_energy(theta, 'A'), m.rocker_energy(theta, 'B'), places=14)
        for theta in (.14,.24,.5):
            self.assertAlmostEqual(m.rocker_force(theta, 'B'), -m.rocker_force(-theta, 'B'), places=14)
            self.assertAlmostEqual(m.rocker_energy(theta, 'B'), m.rocker_energy(-theta, 'B'), places=14)
            self.assertGreater(abs(m.rocker_force(theta,'B')), abs(m.rocker_force(theta,'A')))
        self.assertAlmostEqual(m.rocker_force(.4,'B')/m.rocker_force(.4,'A'),2.5,places=14)

    def test_restoring_potential_is_integral_of_force(self):
        h=1e-6
        for option in ('A','B'):
            for theta in (-.35,-.24,-.17,-.10,0,.10,.17,.24,.35):
                derivative=(m.rocker_energy(theta+h,option)-m.rocker_energy(theta-h,option))/(2*h)
                self.assertAlmostEqual(derivative,-m.rocker_force(theta,option),delta=2e-9)

    def test_transition_continuity_of_force_and_first_two_derivatives(self):
        # Evaluate force, first and second derivatives from either side. Nonlinear
        # smoothstep is C2, so limits exist even at +/- lo and +/- hi.
        h=1e-5
        for x in (-.24,-.10,.10,.24):
            f=lambda x:m.rocker_force(x,'B')
            dleft=(3*f(x)-4*f(x-h)+f(x-2*h))/(2*h)
            dright=(-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h)
            d2left=(2*f(x)-5*f(x-h)+4*f(x-2*h)-f(x-3*h))/h**2
            d2right=(2*f(x)-5*f(x+h)+4*f(x+2*h)-f(x+3*h))/h**2
            self.assertLess(abs(dleft-dright),1e-6)
            self.assertLess(abs(d2left-d2right),.002)

    def test_all_conservative_forces_are_energy_gradient(self):
        q=np.array([.18,.12,-.1,.07,.2,.1,.3,.02])
        h=1e-6; v=np.zeros(8)
        for option in ('A','B'):
            numerical=[]
            for i in range(8):
                d=np.zeros(8); d[i]=h
                numerical.append(-(m.energy(q+d,v,option)-m.energy(q-d,v,option))/(2*h))
            np.testing.assert_allclose(m.conservative_force(q,option),numerical,rtol=2e-9,atol=2e-10)

    def test_bilateral_reciprocity_and_no_one_way_forcing(self):
        q=np.array([.18,.12,-.1,.07,.2,.1,.3,.02]); h=1e-6
        for j in (2,3):
            d=np.zeros(8); d[j]=h
            orb_sensitivity=(m.conservative_force(q+d,'B')[1]-m.conservative_force(q-d,'B')[1])/(2*h)
            d=np.zeros(8); d[1]=h
            fulcrum_sensitivity=(m.conservative_force(q+d,'B')[j]-m.conservative_force(q-d,'B')[j])/(2*h)
            self.assertAlmostEqual(orb_sensitivity,fulcrum_sensitivity,delta=1e-10)
            self.assertGreater(orb_sensitivity,0.)

    def test_target_quadrature_and_smaller_lower_trim(self):
        target_phases=np.degrees(m.P.target_phase)
        np.testing.assert_allclose(target_phases,[0,90,180,270],atol=1e-12)
        self.assertAlmostEqual((target_phases[1]-target_phases[3])%360,180.)
        self.assertLess(m.P.target_amplitude[3], min(m.P.target_amplitude[:3]))
        period=2*math.pi/m.P.omega
        np.testing.assert_allclose(m.targets(.3),m.targets(.3+period),atol=1e-14)

    def test_instantaneous_energy_ledger_independent_gradient(self):
        p=m.P; z=m.initial_state(); z[:8]=[.18,.1,-.11,.14,.2,.3,.24,-.01]
        z[8:16]=[.08,-.11,.12,-.04,.14,-.07,.05,.02]
        step=1e-6
        for option in ('A','B'):
            f=m.rhs(.7,z,option,p)
            # Independent directional energy derivative: include q and v variation.
            dH=(m.energy(z[:8]+step*f[:8],z[8:16]+step*f[8:16],option)
                -m.energy(z[:8]-step*f[:8],z[8:16]-step*f[8:16],option))/(2*step)
            booked=sum(f[m.POS])-sum(f[m.NEG])-f[m.LOSS]
            self.assertAlmostEqual(dH,booked,delta=2e-10)
            self.assertGreaterEqual(f[m.LOSS],0.)

    def test_lossless_unforced_and_damped_unforced(self):
        p=replace(m.P,cycles=1,samples_per_cycle=120,damping=(0.,)*8)
        for option in ('A','B'):
            _,z=m.run(option,p,enabled=False,rtol=1e-10,atol=1e-12,max_step=.05)
            H,r=m.ledger(z,option,p)
            self.assertLess(np.max(abs(H-H[0])),1e-8)
            self.assertEqual(z[-1,m.LOSS],0.)
            self.assertEqual(np.sum(z[-1,m.POS])+np.sum(z[-1,m.NEG]),0.)
            pd=replace(p,damping=m.P.damping)
            _,zd=m.run(option,pd,enabled=False,rtol=1e-10,atol=1e-12,max_step=.05)
            Hd,rd=m.ledger(zd,option,pd)
            self.assertLess(np.max(np.diff(Hd)),1e-9)
            self.assertGreater(zd[-1,m.LOSS],0.)
            self.assertLess(np.max(abs(rd)),1e-8)

    def test_driven_ledger_and_tighter_solver_agree(self):
        p=replace(m.P,cycles=1,samples_per_cycle=120)
        for option in ('A','B'):
            t,z=m.run(option,p)
            _,tight=m.run(option,p,rtol=2e-11,atol=2e-13,max_step=.04)
            H,r=m.ledger(z,option,p)
            self.assertLess(np.max(abs(r)),3e-7)
            self.assertLess(np.max(abs(z[:,:16]-tight[:,:16])),5e-7)
            stages=m.stage_rows(t,z,option,p)
            self.assertEqual(len(stages),4)
            self.assertLess(max(abs(row['ledger_residual']) for row in stages),3e-7)
            self.assertAlmostEqual(sum(row['delta_H'] for row in stages),H[-1]-H[0],places=12)


if __name__ == '__main__':
    unittest.main(verbosity=2)
