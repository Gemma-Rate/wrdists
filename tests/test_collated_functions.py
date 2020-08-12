"""Test vital wrapper_functions"""

import unittest
import numpy as np
import wrdists.collated_functions as cf

class WrapTestCase(unittest.TestCase):
    """
    Test key functions in wrapper_functions.py
    """

    def test_to_decimal(self):
        """
        Test to_decimal. 
        """
        ra, dec = '00 43 28.40', '+64 45 35.40' 
        # RA and DEC for WR1 in hms, dms. 
        ra_new, dec_new = cf.to_decimal(ra, dec)

        self.assertEqual(dec_new[0], '64.76')
        self.assertEqual(ra_new[0], '10.87')


    def test_create_cs(self):
        """
        Test create_cs, which generates coordinate list. 
        """
        ra, dec = '00 43 28.40', '+64 45 35.40' 
        # RA and DEC for WR1 in hms, dms. 
        ra_new, dec_new, c = cf.create_cs([ra], [dec])

        self.assertEqual(np.around(dec_new[0].value, decimals=4), 64.7598)
        self.assertEqual(np.around(ra_new[0].value, decimals=4), 10.8683)


    def test_conv_to_galactic(self):
        """
        Test conv_to_galactic coordintes conversion.
        """  
        ra, dec = 271.18193183732967, -21.158487056124294
        # WR106
        l, b = cf.conv_to_galactic(ra, dec)

        self.assertEqual(np.around(l, decimals=4), 8.8957)
        self.assertEqual(np.around(b, decimals=4), 0.1952)

    
    def test_standardize_coords(self):
        """
        Test RA and DEC conversion to hms and dms format. 
        """

        ra, dec = '00 43 28.39', '+64 45 35.4'
        # RA and DEC for WR1 in hms, dms. 
        ra_new, dec_new = cf.to_decimal(ra, dec, pres=7)
        ra_proc, dec_proc = cf.standardize_coords(ra_new, dec_new)

        self.assertEqual(ra, ra_proc[0])
        self.assertEqual(dec, dec_proc[0])


    def test_gaia_flag(self):
        """
        Test flagging of Gaia astrometry. 
        """

        # Negative parallax flag:
        par, parerr, ast = -0.671, 0.004, 0
        f = cf.gaia_flag(par, parerr, ast)
        self.assertEqual(f, ' n ')   

        # Large parallax error flag:
        par, parerr, ast = 0.671, 0.700, 0
        f = cf.gaia_flag(par, parerr, ast)
        self.assertEqual(f, ' e ')       

        # High astrometric noise:
        par, parerr, ast = 0.671, 0.040, 1.3
        f = cf.gaia_flag(par, parerr, ast)
        self.assertEqual(f, ' a ')     

        # Good value:
        par, parerr, ast = 0.671, 0.040, 0
        f = cf.gaia_flag(par, parerr, ast)
        self.assertEqual(f, ' g ')    


    def test_most_recent(self):
        """
        Test obtaining most recent spectral type.
        """
        
        wr_type = 'WC4; WC5; WC4pd+O9.6'
        ref = 'SS90; CDB98; WH96'
        # Information for WR19. 
        wrecent = cf.most_recent(wr_type, ref)

        self.assertEqual(wrecent, 'WC5')


if __name__ == '__main__':
    unittest.main()