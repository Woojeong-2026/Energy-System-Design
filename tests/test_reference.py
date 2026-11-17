"""Run with python -m unittest discover -s tests -p 'test_*.py'."""
import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reference'))
from charging_reference import Thermo, solve, audit, case_list


class ChargingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.th=Thermo(); cls.base=solve(th=cls.th)

    def test_prescribed_cases(self):
        for label,c in case_list():
            with self.subTest(case=label):
                rows,_=audit(solve(c,self.th),self.th)
                self.assertTrue(all(r['Pass'] for r in rows))

    def test_hand_inventory_and_electron_count(self):
        p=self.base
        self.assertAlmostEqual(p['n'][3][1],20/7,places=10)
        self.assertAlmostEqual(p['n'][10][4],1)
        self.assertAlmostEqual(p['n'][12][2],1)
        self.assertEqual(p['metrics']['I'],385940000)
        self.assertEqual(p['metrics']['Q_link'],0)
        self.assertAlmostEqual(p['metrics']['P_in']/1e6,637.958716,delta=1e-5)

    def test_scale_invariance(self):
        for f in (.01,2,10):
            p=solve(dict(nCO2_feed=f),self.th)
            self.assertAlmostEqual(p['metrics']['eta_ch'],self.base['metrics']['eta_ch'],places=9)
            self.assertAlmostEqual(p['metrics']['P_in']/self.base['metrics']['P_in'],f,places=8)

    def test_thermal_link_boundary(self):
        a=solve(dict(TBch_C=869.999),self.th); b=solve(dict(TBch_C=870),self.th)
        self.assertEqual(a['metrics']['Q_link'],0)
        self.assertGreater(b['metrics']['Q_link'],0)
        self.assertTrue(all(r['Pass'] for r in audit(b,self.th)[0]))

    def test_zero_heat_recovery(self):
        p=solve(dict(eff_hex=0,f_link=0),self.th)
        self.assertEqual(p['metrics']['Q_hex'],0)
        self.assertAlmostEqual(p['T'][2],298.15,places=6)
        self.assertGreater(p['metrics']['P_in'],self.base['metrics']['P_in'])

    def test_invalid_inputs(self):
        for c in [dict(X_soec=0),dict(X_soec=1),dict(nCO2_feed=0),dict(P_bar=-1),
                  dict(eff_hex=1.1),dict(Tcell_C=float('nan')),dict(p0=101325)]:
            with self.subTest(config=c),self.assertRaises(ValueError): solve(c,self.th)

    def test_faults_are_detected(self):
        for fault in ('carbon','current','link','voltage','heat'):
            p=copy.deepcopy(self.base)
            if fault=='carbon': p['n'][10][4]=.9
            if fault=='current': p['metrics']['I']/=1000
            if fault=='link': p['metrics']['Q_link']=1e6
            if fault=='voltage': p['metrics']['V']-=.2
            if fault=='heat': p['metrics']['P_in']-=p['metrics']['Q_hex']
            with self.subTest(fault=fault): self.assertFalse(all(r['Pass'] for r in audit(p,self.th)[0]))

    def test_malformed_export_is_rejected(self):
        for fault in ('shape','nan','missing_config'):
            p=copy.deepcopy(self.base)
            if fault=='shape': p['n']=p['n'][:-1]
            if fault=='nan': p['T'][0]=float('nan')
            if fault=='missing_config': del p['cfg']['P_bar']
            with self.subTest(fault=fault),self.assertRaises(ValueError): audit(p,self.th)


if __name__=='__main__': unittest.main()
