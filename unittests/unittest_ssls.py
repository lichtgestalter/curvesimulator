from unittest import TestCase
from cs_physics import CurveSimPhysics


class TestSsls(TestCase):
    # Tests may not run in the order of their appearance!

    @classmethod
    def setUpClass(cls):
        # automatically executed before the first test in this class
        pass

    @classmethod
    def tearDownClass(cls):
        # automatically executed after the last test in this class
        pass

    def setUp(self):
        # automatically executed before every single test in this class
        self.ea, self.e, self.ma = 0.0, 0.0, 0.0

    def tearDown(self):
        # automatically executed after every single test in this class
        pass

    @staticmethod
    def does_not_do_anything():
        print("This will not be printed because the function does not start with \'test\'")

    def test_keplers_equation_0(self):
        result = CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)
        self.assertEqual(result, 0)

    def test_keplers_equation_rangecheck_ea(self):
        self.ea = -7.0
        with self.assertRaises(ValueError):
            CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)
        self.ea = 7.0
        with self.assertRaises(ValueError):
            CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)

    def test_keplers_equation_rangecheck_e(self):
        self.e = 1.0
        with self.assertRaises(ValueError):
            CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)
        self.e = -0.1
        with self.assertRaises(ValueError):
            CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)

    def test_keplers_equation_rangecheck_ma(self):
        self.ma = -7.0
        with self.assertRaises(ValueError):
            CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)
        self.ma = 7.0
        with self.assertRaises(ValueError):
            CurveSimPhysics.kepler_equation(self.ea, self.e, self.ma)
