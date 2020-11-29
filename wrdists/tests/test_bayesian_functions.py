"""Test cases for bayesian calculation class"""

import wrdists.bayesian_functions as bc
import unittest
import numpy as np
import matplotlib.pyplot as plt
import wrdists.collated_functions as cf

class DistribTestCase(unittest.TestCase):
    """
    Test all functions associated with the bayesian Distribution object
    in bayesian_cal.py
    """
    
    def create_dist(self, r=15000, n=15000, par=0.1, err=0.01):
        """
        Install distribution.
        """
        test_r = np.linspace(0, r, n)
        test_dist = bc.Distribution(par, err, test_r)
        # Create distribution object. 

        return test_dist

    def test_gauss_simple(self):
        """
        Test gauss_simple produces a gaussian as expected.
        """
        test_dist = self.create_dist()
        result1 = test_dist.gauss_simple(0.5, 1, 0)
        # Generate gaussian result for parameters given above.

        self.assertEqual(np.around(result1, decimals=10), 0.3520653268)
        # Test gaussian compared to manual calculation.
        result2 = test_dist.gauss_simple(-0.5, 1, 0)
        # Test negative distribution.
        self.assertEqual(result1, result2)

    def test_correct_gaia(self):
        """
        Test correction for mas, zero point and uncertainties is working.
        """
        test_dist = self.create_dist(par=0.1, err=0.01)
        test_dist.correct_gaia(10, uwu=False, zpt=-0.029)

        self.assertTrue(np.isclose(test_dist.err, 1.24465e-5, rtol=1e-7))
        self.assertTrue(np.isclose(test_dist.dpt, 1.29e-4))

    def test_radial_dist_height(self):
        """
        Test radial distance from galactic centre and height is correct.
        """
        test_dist = self.create_dist()
        b = np.degrees(-np.arctan(test_dist.sun_height/test_dist.gal_cent_dist))
        rg, z = test_dist.radial_dist_height(0, b, 8122)

        self.assertEqual(rg, 0.0)
        self.assertEqual(np.around(z, decimals=1), 0.0)

        # Test WR61, compared to Rosslowe and Crowther 2015:
        test_dist.sun_height = 20
        test_dist.gal_cent_dist = 8000
        rg, z = test_dist.radial_dist_height(311.28, -3.91, 6320)
        # Change input parameters to be the same as the paper.

        self.assertEqual(np.around(rg/1e3, decimals=2), 6.10)
        self.assertEqual(np.around(z, decimals=0), -411)

    def test_total_dust(self):
        """
        Test dust generation model
        """
        
        test_dist = self.create_dist()

        # Test horizontal dust vs manually calculated results:
        d, mol, atom, molh, atomh = test_dist.total_dust(2000, 0)
        self.assertEqual(np.around(mol, decimals=4), 0.9209)
        self.assertEqual(np.around(atom, decimals=4), 0.0800)

        # Test vertical dust vs manually calculated results:
        d, mol, atom, molh, atomh = test_dist.total_dust(15000, 500)
        self.assertEqual(np.around(atomh, decimals=1), 0.5)
        d, mol, atom, molh, atomh = test_dist.total_dust(15000, 100)
        self.assertEqual(np.around(atomh, decimals=4), 0.9696)
        # Atomic gas.
        d, mol, atom, molh, atomh = test_dist.total_dust(10000, 500)
        self.assertEqual(np.around(molh, decimals=6), 2.23e-4)
        d, mol, atom, molh, atomh = test_dist.total_dust(10000, 100)
        self.assertEqual(np.around(molh, decimals=4), 0.4333)
        # Molecular gas. 

    def test_errs(self):
        """
        Test process to generate 68% credible intervals
        """
        test_dist = self.create_dist(r=2)
        std_test, avg = 0.2, 1
        gauss = test_dist.gauss_simple(test_dist.r_range, std_test, avg)
        bign, bigr, n, interval = test_dist.errs(gauss, 0.68)

        stds = np.array([np.abs(interval[0]-avg), np.abs(interval[1]-avg)])
        self.assertEqual(np.around(stds[0], decimals=1), std_test)
        self.assertEqual(np.around(stds[1], decimals=1), std_test)
        # Check standard deviations are equal to std_test above.

    def test_dust_model(self):
        """
        Test that the extinctions produced along the line of sight are sensible.
        """
        l, b = 2.38, 1.41
        # Coordinates for WR102 (galactic centre direction)
        test_dist = self.create_dist(r=2639)
        ai, delta_p, dtot, dmol = test_dist.dust_model(l, b)
        self.assertEqual(np.around(ai[-1], decimals=1), 1.4)

        # E(B-V)~0.91 from http://argonaut.skymaps.info/, AI=1.35, based on 
        # interpolation from 12 to 12.5 to get best fit gradient 
        # (0.750-1.113)/0.5.

        # Check galactic centre:
        l = 0
        test_dist = self.create_dist(r=8122)
        b = np.degrees(-np.arctan(test_dist.sun_height/test_dist.gal_cent_dist))
        ai, delta_p, dtot, dmol = test_dist.dust_model(l, b)
        self.assertEqual(np.around(ai[-1], decimals=1), 
                         np.around(0.48*32, decimals=1))

        # Check very close to the sun:
        test_dist = self.create_dist(r=1)
        b = np.degrees(-np.arctan(test_dist.sun_height/test_dist.gal_cent_dist))
        ai, delta_p, dtot, dmol = test_dist.dust_model(l, b)
        self.assertEqual(np.around(ai[-1], decimals=2), 0)
        

    def test_known_value(self):
        """
        Test full calculation with a known result.
        """
        wr_dist = self.create_dist(par=0.5, err=0.02)
        # Distance = 2kpc.

        wr_dist.dpt, wr_dist.err = wr_dist.dpt*1e-3, wr_dist.err*1e-3
        # Make data mas scale...
    
        l, b = cf.conv_to_galactic(180, 0)

        n, nl, nb, sig, mu = wr_dist.hii_numbers(l, b)
        # Run depending on HII regions along line of sight. 

        av, delta_p, dtot, dmol = wr_dist.dust_model(l, b)
        # Dust model. 

        normedg, maximum_r = wr_dist.apply_wr_prior(mu=mu, sigma=sig, delta=delta_p)
        # Get unnormalised posterior from gaussian prior. 

        self.assertEqual(np.around(maximum_r/1e3, decimals=2), 2.00)


if __name__ == '__main__':
    unittest.main()