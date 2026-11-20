#!/usr/bin/env python3
"""Mathematical, numerical, and on-disk reproducibility checks; run after experiment."""
import csv
import hashlib
import json
import math
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import hybrid_experiment as m
import reference_model as ref


class HybridMathematics(unittest.TestCase):
    def test_reference_source_hash(self):
        path=Path(__file__).resolve().parent/'reference_model.py'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
            'f067a35c52327bbfaab842b710dd2087377cf51b64a60aa8806518640494dc2f')

    def test_exact_endpoints_and_below_threshold_all_weights(self):
        for th in (-.4,-.24,-.17,-.1,-.03,0.,.03,.1,.17,.24,.4):
            for w,option in ((0.,'A'),(1.,'B')):
                self.assertEqual(m.rocker_force(th,w),ref.rocker_force(th,option))
                self.assertEqual(m.rocker_energy(th,w),ref.rocker_energy(th,option))
            for w in m.WEIGHTS:
                self.assertAlmostEqual(m.rocker_force(-th,w),-m.rocker_force(th,w),places=14)
                self.assertAlmostEqual(m.rocker_energy(-th,w),m.rocker_energy(th,w),places=14)
                if abs(th)<=.1:
                    self.assertAlmostEqual(m.rocker_force(th,w),ref.rocker_force(th,'A'),places=15)
                    self.assertAlmostEqual(m.rocker_energy(th,w),ref.rocker_energy(th,'A'),places=15)

    def test_explicit_convex_formula_and_high_angle_stiffness(self):
        for w in m.WEIGHTS:
            for th in (-.35,-.19,-.05,.05,.19,.35):
                a,b=ref.rocker_force(th,'A'),ref.rocker_force(th,'B')
                ua,ub=ref.rocker_energy(th,'A'),ref.rocker_energy(th,'B')
                self.assertAlmostEqual(m.rocker_force(th,w),(1-w)*a+w*b,places=15)
                self.assertAlmostEqual(m.rocker_energy(th,w),(1-w)*ua+w*ub,places=15)
            self.assertAlmostEqual(m.rocker_force(.3,w)/ref.rocker_force(.3,'A'),1+1.5*w,places=14)

    def test_rocker_gradient_across_regimes_and_boundaries(self):
        h=1e-6
        for w in m.WEIGHTS:
            for th in (-.35,-.24,-.17,-.10,-.01,0.,.10,.17,.24,.35):
                derivative=(m.rocker_energy(th+h,w)-m.rocker_energy(th-h,w))/(2*h)
                self.assertAlmostEqual(derivative,-m.rocker_force(th,w),delta=2e-9)

    def test_transition_has_matching_first_and_second_force_derivatives(self):
        h=1e-5
        for w in m.WEIGHTS:
            for x in (-.24,-.10,.10,.24):
                f=lambda y:m.rocker_force(y,w)
                dl=(3*f(x)-4*f(x-h)+f(x-2*h))/(2*h)
                dr=(-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h)
                ddl=(2*f(x)-5*f(x-h)+4*f(x-2*h)-f(x-3*h))/h**2
                ddr=(2*f(x)-5*f(x+h)+4*f(x+2*h)-f(x+3*h))/h**2
                self.assertLess(abs(dl-dr),1e-6)
                self.assertLess(abs(ddl-ddr),.002)

    def test_all_eight_force_components_equal_potential_gradient(self):
        h=1e-6
        q=np.array([0.,.12,-.1,.07,.2,.1,.3,.02]); v=np.zeros(8)
        for w in m.WEIGHTS:
            for th in (-.24,-.17,-.1,0.,.1,.17,.24):
                q[0]=th; numerical=[]
                for i in range(8):
                    d=np.zeros(8); d[i]=h
                    numerical.append(-(m.energy(q+d,v,w)-m.energy(q-d,v,w))/(2*h))
                np.testing.assert_allclose(m.conservative_force(q,w),numerical,rtol=3e-9,atol=3e-10)

    def test_endpoint_vector_fields_and_potential_exact(self):
        for th in (-.3,-.17,-.1,0.,.1,.17,.3):
            z=m.initial_state(); z[:8]=[th,.1,-.11,.14,.2,.3,.24,-.01]
            z[8:16]=[.08,-.11,.12,-.04,.14,-.07,.05,.02]
            for w,opt in ((0.,'A'),(1.,'B')):
                np.testing.assert_array_equal(m.conservative_force(z[:8],w),ref.conservative_force(z[:8],opt))
                self.assertEqual(m.energy(z[:8],z[8:16],w),ref.energy(z[:8],z[8:16],opt))
                np.testing.assert_array_equal(m.rhs(.7,z,w),ref.rhs(.7,z,opt))

    def test_instantaneous_energy_ledger_independent_directional_derivative(self):
        h=1e-6
        for w in m.WEIGHTS:
            for th in (-.3,-.17,-.1,0.,.1,.17,.3):
                z=m.initial_state(); z[:8]=[th,.1,-.11,.14,.2,.3,.24,-.01]
                z[8:16]=[.08,-.11,.12,-.04,.14,-.07,.05,.02]
                dz=m.rhs(.7,z,w)
                dh=(m.energy(z[:8]+h*dz[:8],z[8:16]+h*dz[8:16],w)-
                    m.energy(z[:8]-h*dz[:8],z[8:16]-h*dz[8:16],w))/(2*h)
                booked=sum(dz[m.POS])-sum(dz[m.NEG])-dz[m.LOSS]
                self.assertAlmostEqual(dh,booked,delta=3e-10)
                self.assertTrue(np.all(dz[16:]>=0.))

    def test_pairs_and_weights_are_controlled(self):
        base=m.initial_state(); plus=m.initial_state(1); minus=m.initial_state(-1)
        d=np.zeros(25); d[0]=d[1]=.04
        np.testing.assert_allclose(plus-base,d,atol=1e-17)
        np.testing.assert_allclose(minus-base,-d,atol=1e-17)
        np.testing.assert_array_equal(plus,ref.initial_state(True))
        for w in (-.1,1.1,float('nan'),float('inf')):
            with self.assertRaises(ValueError): m.rocker_force(.2,w)
        with self.assertRaises(ValueError): m.initial_state(2)
        period=2*math.pi/m.P.omega
        np.testing.assert_allclose(ref.targets(.3),ref.targets(.3+period),atol=1e-14)

    def test_short_unforced_lossless_and_passive_damping(self):
        p=replace(m.P,cycles=1,samples_per_cycle=120,damping=(0.,)*8)
        for w in m.WEIGHTS:
            _,z=m.run(w,p,enabled=False); H,R=m.ledger(z,w,p)
            self.assertLess(np.max(np.abs(H-H[0])),1e-8)
            self.assertEqual(float(np.sum(z[-1,16:])),0.)
            pd=replace(p,damping=m.P.damping)
            _,z=m.run(w,pd,enabled=False); H,R=m.ledger(z,w,pd)
            self.assertLess(np.max(np.diff(H)),1e-9)
            self.assertGreater(z[-1,m.LOSS],0.)
            self.assertLess(np.max(np.abs(R)),1e-8)

    def test_short_endpoint_integrations_recover_reference_exactly(self):
        p=replace(m.P,cycles=1,samples_per_cycle=120)
        for w,opt in ((0.,'A'),(1.,'B')):
            t,z=m.run(w,p); tr,zr=ref.run(opt,p)
            np.testing.assert_array_equal(t,tr); np.testing.assert_array_equal(z,zr)

    def test_short_midpoint_tight_solver_and_stage_energy(self):
        p=replace(m.P,cycles=1,samples_per_cycle=120)
        t,z=m.run(.5,p); _,zt=m.run(.5,p,tight=True)
        self.assertLess(np.max(np.abs(z[:,:16]-zt[:,:16])),5e-7)
        self.assertLess(np.max(np.abs(z[:,16:]-zt[:,16:])),1e-4)
        H,R=m.ledger(z,.5,p); stages=m.stage_rows(t,z,.5,p)
        self.assertEqual(len(stages),4)
        self.assertLess(np.max(np.abs(R)),1e-7)
        self.assertLess(max(abs(s['ledger_residual']) for s in stages),1e-7)
        self.assertAlmostEqual(sum(s['delta_H'] for s in stages),H[-1]-H[0],places=12)


class RecordedOutput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=Path(__file__).resolve().parent/'output'
        if not (cls.out/'summary.json').exists():
            raise RuntimeError('Run hybrid_experiment.py before testing the generated output')
        cls.summary=json.loads((cls.out/'summary.json').read_text())

    def test_sha256_inventory(self):
        lines=(self.out/'SHA256SUMS.txt').read_text().splitlines()
        recorded={}
        for line in lines:
            digest,name=line.split('  ',1); recorded[name]=digest
            self.assertEqual(hashlib.sha256((self.out/name).read_bytes()).hexdigest(),digest,name)
        actual={p.name for p in self.out.iterdir() if p.is_file()}-{'SHA256SUMS.txt'}
        self.assertEqual(set(recorded),actual)

    def test_all_complete_saved_runs_and_paired_results(self):
        s=self.summary
        self.assertEqual(s['endpoint_reference_trajectory_bitwise_equal'],{'A':True,'B':True})
        with np.load(self.out/'trajectories.npz',allow_pickle=False) as a:
            expected={'time'}
            for w in m.WEIGHTS:
                kk=m.key(w)
                expected.update(kk+'__'+v for v in m.VARIANTS)
                if w in (0.,1.): expected.add(kk+'__reference_endpoint')
                base=a[kk+'__baseline_standard']; self.assertEqual(base.shape,(4801,25))
                np.testing.assert_array_equal(base[0],m.initial_state())
                np.testing.assert_array_equal(a[kk+'__plus_standard'][0],m.initial_state(1))
                np.testing.assert_array_equal(a[kk+'__minus_standard'][0],m.initial_state(-1))
                for name in m.VARIANTS:
                    z=a[kk+'__'+name]
                    self.assertTrue(np.isfinite(z).all())
                    diag=s['integrity']['monotonicity_by_run'][kk+'__'+name]
                    work=np.diff(z[:,16:24],axis=0)
                    passive=np.diff(z[:,24],axis=0)
                    self.assertEqual(float(np.min(work)),diag['minimum_work_increment'])
                    self.assertEqual(int(np.sum(work < -1e-9)),diag['work_decreases_exceeding_1e_9'])
                    self.assertEqual(int(np.sum(passive < -1e-9)),diag['passive_decreases_exceeding_1e_9'])
                    _,res=m.ledger(z,w); self.assertLess(float(np.max(np.abs(res))),1e-7)
                for name in ('baseline','plus'):
                    normal=a[kk+'__'+name+'_standard']; tight=a[kk+'__'+name+'_tight']
                    c=s['comparison'][kk]['normal_vs_tight_'+name]
                    self.assertEqual(float(np.max(np.abs(normal[:,:16]-tight[:,:16]))),c['max_physical_state_component'])
                    self.assertEqual(float(np.max(np.abs(normal[:,16:]-tight[:,16:]))),c['max_ledger_state_component'])
                    self.assertLess(c['max_physical_state_component'],5e-7)
                    self.assertLess(c['max_ledger_state_component'],1e-4)
                for sign in ('plus','minus'):
                    d=np.linalg.norm(a[kk+'__'+sign+'_standard'][:,:16]-base[:,:16],axis=1)
                    self.assertAlmostEqual(float(d[-1]),s['comparison'][kk]['perturb_'+sign]['end_cycle_20'],places=16)
                if w in (0.,1.): np.testing.assert_array_equal(base,a[kk+'__reference_endpoint'])
            self.assertEqual(set(a.files),expected)
            self.assertEqual(a['time'].shape,(4801,))
            self.assertTrue(np.all(np.diff(a['time'])>0.))
        diags=s['integrity']['monotonicity_by_run'].values()
        self.assertEqual(s['integrity']['work_monotone_with_tolerance_1e_9'],
                         all(d['work_decreases_exceeding_1e_9']==0 for d in diags))
        diags=s['integrity']['monotonicity_by_run'].values()
        self.assertEqual(s['integrity']['passive_monotone_with_tolerance_1e_9'],
                         all(d['passive_decreases_exceeding_1e_9']==0 for d in diags))

    def test_csv_trajectories_stage_counts_and_ledgers(self):
        with np.load(self.out/'trajectories.npz',allow_pickle=False) as a:
            for w in m.WEIGHTS:
                kk=m.key(w); raw=np.genfromtxt(self.out/(kk+'_baseline.csv'),delimiter=',',skip_header=1)
                self.assertEqual(raw.shape,(4801,28))
                np.testing.assert_array_equal(raw[:,0],a['time'])
                np.testing.assert_array_equal(raw[:,1:26],a[kk+'__baseline_standard'])
                H,R=m.ledger(a[kk+'__baseline_standard'],w)
                np.testing.assert_array_equal(raw[:,26],H); np.testing.assert_array_equal(raw[:,27],R)
        with (self.out/'stages.csv').open(newline='') as f: rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),400)
        for w in m.WEIGHTS:
            selected=[r for r in rows if float(r['lambda_weight'])==w]
            self.assertEqual(len(selected),80)
            self.assertLess(max(abs(float(r['ledger_residual'])) for r in selected),1e-7)
        with (self.out/'metrics.csv').open(newline='') as f: self.assertEqual(len(list(csv.DictReader(f))),5)


if __name__=='__main__': unittest.main(verbosity=2)
