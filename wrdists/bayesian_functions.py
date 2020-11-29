"""Class for obtaining the WR star distance distributions"""

import numpy as np

class Distribution(object):
    def __init__(self, data, error, r=None):
        self.dpt = data
        # Parallax data point, which is the average point of the distribution.
        self.err = error
        # Errors of parallax data points, which are the standard deviations of the
        # distribution.
        self.r_range = r
        # Range of possible values over which to plot output
        # distribution.

        self.gal_cent_dist = 8122
        # Galactic centre distance in pc.
        self.sun_height = 20.8
        # Height of sun above galactic plane in pc.
        # From Bennett and Bovy, (2019).

        self.dist = None
        # Distribution output probabilities .

        self.prior = None
        self.dust_component = None
        self.hii_regions = None
        # Prior components.

    @staticmethod
    def gauss_simple(x, sigma, avg):
        """
        Simple gaussian helper function.

        Parameters
        -------
        x : array
            Input array for the gaussian.
        sigma : float
            Standard deviation.
        avg : float
            Average of the gaussian.

        Returns
        -------
        gauss : array
            Gaussian probability distribution for the data.

        """
        gauss = float(1) / np.sqrt(float(2) * np.pi * sigma ** float(2)) * \
                np.exp(float(-1) / (float(2) * sigma ** float(2)) * (x - avg) \
                       ** float(2))

        return gauss

    def likelihood(self, avg_change=1, range_dist=4):

        """
        Plot gaussian likelihood of the data distribution.

        Parameters
        -------
        avg_change : float
            Alter distribution average position (to illustrate sampling
            from distribution).
        range_dist : float
            Change distribution +/- range.

        Returns
        -------
        gauss : array
            Probability distribution for the data, given the error and
            average.

        """

        navg = avg_change * self.dpt
        # Alter average position.

        x = np.linspace(navg - (range_dist * self.err),
                        navg + (range_dist * self.err))
        # Make x values for plot.

        gauss = self.gauss_simple(x, self.err, navg)
        gauss_norm = gauss / np.trapz(gauss, x)
        # Normalise gaussian.

        return x, gauss_norm

    def exp_decrease(self, L):
        """
        Assign exponentially decreasing prior to r range and parallax.
        Best for objects >2kpc away (Bailer-Jones (2015)).

        Parameters
        -------
        L : float
            Length scale of exponential decrease (pc).
        r : float
            Range of distances for the WR distribution.

        Returns
        -------
        modes : array
            Distance to WR in pc.

        """

        L = L * 1e3
        # L in pc.

        r_pos = self.r_range[self.r_range > 0]
        unnorm_pos = r_pos ** 2 * np.exp(-r_pos / L) / (self.err) * np.exp(
            -1 / (2 * self.err ** 2) * (self.dpt - 1 / r_pos) ** 2)  #
        # Positive values of r
        r_neg = self.r_range[self.r_range <= 0]
        unnorm_neg = 0 * r_neg
        # Negative values of r

        unnorm = np.concatenate((unnorm_neg, unnorm_pos))
        # Unnormalised posterior.

        Z = np.trapz(unnorm, self.r_range)
        # Integrated normalisation.

        normed = unnorm / Z
        # Normalised posterior. Area under curve sums to one.

        coeffs = [1 / L, -2, self.dpt / self.err ** 2, -1 / self.err ** 2]
        roots = np.roots(coeffs)
        # Calculate the modes by taking roots of differentiated
        # distribution, according to Bailer-Jones (2015).

        # Select modes according to real roots:
        real_roots = roots[np.isreal(roots)]

        try:
            roots[1]
            if self.dpt >= 0:
                mode = np.min(real_roots)
            elif self.dpt < 0:
                mode = real_roots[real_roots > 0]
                # Select the roots.
        except IndexError:
            mode = real_roots[0]
            # Select only real root.

        full_mode = np.real(mode)  # /1000
        # Selected max likelihood in pc.

        return full_mode, real_roots, normed, Z

    def errs(self, normed, cred_range):
        """
        Get r68 errors and scale to 1 sigma errors (Astraatmadja + Bailer-Jones (2016)), using credible intervals.

        Parameters
        -------
        normed : array
            Probabilities.
        cred_range : float
            Credible interval value to use.

        Returns
        -------
        greater_n : array
            Values of normed greater than the credible interval.
        greater_r : array
            Corresponding values of r to greater_n.
        n : array
            Value of normed at the credible interval.
        intervals : list (of tuples)
            Minimum, maximum distance for intercepts of each peak (around each mode). For a single peak distribution,
            (e.g the WR prior) this corresponds to the upper and lower limits of the distance uncertainties.

        """
        cred = []
        # Empty array to fill.

        for n in sorted(list(normed), reverse=True):
            # Loop through indices backwards.

            greater_n = normed[normed > n]
            # Values greater than n for each value of normed.
            greater_r = self.r_range[normed > n]
            # r values in this range.

            inter = np.trapz(greater_n, greater_r)
            # Integrate in range of normed data up to ith value.

            if inter > cred_range:
                # Stop when cred_range probability is reached.
                break

        # Get all intercepts:

        linear = np.around(self.r_range[1] - self.r_range[0])
        # Get linear difference between consecutive r.
        diffs = np.around(np.diff(greater_r))
        diffs_min = np.where(diffs > linear)

        try:
            diffs_max = diffs_min[0] + 1
            # Find where differences are particularly large.
            # e.g, difference of 99 between 55 and 56th elements. Need to
            # select both 55 and 56. diffs_min are first r around a curve
            # diffs_max are second encountered r.
            # Case where there are only two intercepts (i.e diffs_min is an
            # empty array).

            int_min, int_max = [greater_r[0]] + list(greater_r[diffs_min]), \
                               list(greater_r[diffs_max]) + [greater_r[-1]]

        except TypeError:
            # Case where there are only two intercepts (diffs_min is an
            # empty array)
            int_min, int_max = greater_r[0], greater_r[-1]

        intervals = [int_min[0], int_max[0]]

        return greater_n, greater_r, n, intervals
          

    def correct_gaia(self, gphot, uwu=False, zpt=-0.029):
        """
        Account for random errors and parallax error underestimates
        (compared to external catalogues).

        Parameters
        -------
        gphot : float
            Gaia G band required for choosing right correction to apply.
        uwu : boolean or float
            Either False for calculation of weighting to be applied using
            external catalogue reweighting (table 1 in Arenou et al. 2018)
            from the Gaia catalogue validation paper, or float for user
            supplied reweighting values.
        zpt : float
            Either Default is to apply global zero point of -0.029 from Gaia
            astrometric solution paper (Lindegren et al. 2018) or other float
            for user supplied zero point.
        """

        self.dpt = (self.dpt - zpt) * 1e-3
        # Correct parallax data point for global zero point error.

        if uwu:
            self.err = self.err * uwu * 1e-3
            # Apply user supplied reweighting.
        else:
            # Apply function for weighting error.
            uwu = -0.01319 * gphot + 1.37552 + \
                  self.gauss_simple(gphot, 1.35373, 14.58974) * 1.1
            self.err = self.err * uwu * 1e-3

    def hii_numbers(self, l, b):
        """
        Get number of Hii regions at galactic latitude/longitude, via fits to
        the graphs of synthetic catalogue in figure 6 of Paladini et al (2003).
        For normalisation, assume number at a given latitude and longitude are
        spread evenly throughout the latitude and longitude.

        Parameters
        -------
        l : float
            Galactic longitude coordinate.
        b : float
            Galactic latitude coordinate.

        Returns
        -------
        nl : float
            Estimated fraction of Hii regions at longitutde coordinate.
        nb : float
            Estimated fraction of Hii regions at latitude coordinate.
        n : float
            Estimated number of Hii regions at coordinate.

        """

        # Galactic latitude:
        gauss_normed = self.gauss_simple(b, 0.4, -0.25) * 293 * 0.8
        gauss_normed_2 = self.gauss_simple(b, 1.6, 0.15) * 293 * 0.75
        nb = (gauss_normed + gauss_normed_2) / 14
        # Normalised number at latitude.

        # Galactic longitude:
        if l > 180:
            l = l - 360
        # Convert to -180 to 180 from 0 to 360.

        gauss_normed = self.gauss_simple(l, 20, 30) * 76 * 20
        gauss_normed_2 = self.gauss_simple(l, 30, -20) * 76 * 20
        gauss_normed_3 = self.gauss_simple(l, 3, 79.5) * 76 * 7
        gauss_normed_4 = self.gauss_simple(l, 10, 110) * 76 * 3
        gauss_normed_5 = self.gauss_simple(l, 10, -165) * 76 * 4
        gauss_normed_6 = self.gauss_simple(l, 20, 155) * 76 * 2

        nl = (gauss_normed + gauss_normed_2 + gauss_normed_3 + gauss_normed_4 + gauss_normed_5 + gauss_normed_6) / 360
        # Normalised number at longitude.

        # Get total number:
        n = nl * nb
        # Get estimated number density on line of sight.

        # Rescale sigma:
        sigma = 3000 / (n + 1)

        # Relocate mu:
        mu = 3000

        if l < 86 and l > 73 and b > -3 and b < 4:
            mu = 1400
            sigma = 1400 / n + 1
        # Cygnus x region peak.

        return n, nl, nb, sigma, mu

    def radial_dist_height(self, l, b, r):
        """
        Get radial distance from galactic centre and height above plane.

        Parameters
        -------
        l : float
            Galactic longitude coordinate.
        b : float
            Galactic latitude coordinate.
        r : array or float
            Distance of object from sun.

        Returns
        -------
        re : float
            Radial distance from the sun.
        rg : float
            Distance from galactic centre.
        z : float
            Height above plane.
        """
        # For this calculation, assume we are in the galactic plane,
        # and then add 20.8pc height at:

        gal_cent_angle = np.degrees(-np.arctan(self.sun_height / self.gal_cent_dist))
        ro = self.gal_cent_dist * np.cos(np.radians(gal_cent_angle))
        # Radial distance of galactic centre from the sun in pc.

        re = r * np.cos(np.radians(b))
        # Radial distance from sun (projected).
        rg = np.sqrt(ro ** 2 + re ** 2 - 2 * re * ro * np.cos(np.radians(l)))
        # Radial distance from galactic centre.
        z = r * np.sin(np.radians(b)) + float(self.sun_height)

        return rg, z

    def total_dust(self, rg, z):
        """
        Calculate atomic and molecular dust in midplane and height.
        Then get total dust. From the model in Rosslowe and Crowther (2015).

        Parameters
        -------
        rg : array or float
            Radial distances from galactic centre (pc).
        z : array or float
            Height above plane (pc).

        Returns
        -------
        Dtot : array or float
            Total dust content at each grid point along the line of sight.
        Dmol : array or float
            Molecular gas in the plane at each distance point from galactic
            centre (rg).
        Datom : array or float
            Atomic gas in the plane at each distance point from galactic
            centre (rg).
        Dmolh : array or float
            Molecular gas at height z above the plane.
        Datomh : array or float
            Atomic gas at height z above the plane.
        """

        pmol, patom = np.poly1d(np.polyfit([0, 10000], [25, 90], 1)), \
                      np.poly1d(np.polyfit([0, 15000], [100, 500], 1))
        # Polynomial fit to molecular and atomic Z1/2 data.
        z12mol, z12atom = pmol(rg), patom(rg)
        # Get scale heights.

        nmol, natom = 10, 0.08  # cm^-3
        # Molecular number density.
        Dmol = nmol * 1 / (np.cosh(rg / 800)) ** 2 + np.exp(-(rg - 4300) ** 2 / (2 * 2500 ** 2))
        # Midplane density of molecular gas.
        Datom = natom * (1 + 1.3 / (1 + np.exp(-(rg - 6500) / 200)) - 1.3 / (1 + np.exp(-(rg - 13200) / 550)))
        # Midplane density of atomic gas.
        Dmolh = 1 / (np.cosh(np.log(1 + np.sqrt(2)) * z / z12mol)) ** 2
        Datomh = 1 / (np.cosh(np.log(1 + np.sqrt(2)) * z / z12atom)) ** 2
        # Vertical dependence of both molecular and atomic gas.
        Dtot = (Dmol * Dmolh) + (Datom * Datomh)
        # Dust function.

        return Dtot, Dmol, Datom, Dmolh, Datomh

    def gal_cent_dust(self):
        """
        Get total dust integrated to galactic centre (for normalisation).
        """
        b = np.degrees(-np.arctan(self.sun_height / self.gal_cent_dist))
        # Get direction of galactic centre.
        r_gal = np.linspace(0, self.gal_cent_dist, self.gal_cent_dist)
        # Get points along line of sight.
        rg, z = self.radial_dist_height(0, b, r_gal)
        Dtot, Dmol, Datom, Dmolh, Datomh = self.total_dust(rg, z)
        # Dust along line of sight.

        bh = []
        for i, ri in enumerate(r_gal):
            bh.append(np.trapz(Dtot[0:i], r_gal[0:i]))
            # Total dust (integrated along line of sight)

        centre_dust = bh[-1]
        # Dust at galactic centre.

        return centre_dust

    def dust_model(self, l, b):
        """
        Apply dust model from Rosslowe and Crowther (2015) and then calibrate
        to extinction in I band at galactic centre.

        These extinction magnitudes are then applied to the flux differences
        using 2.512**-AI (where flux is proportional to probability).
        This generates a list of corrections to apply to each r
        value in the prior.

        Parameters
        -------
        l : float
            Galactic longitude coordinate.
        b : float
            Galactic latitude coordinate.

        Returns
        -------
        nl : float
            Estimated fraction of Hii regions at longitutde coordinate.
        nb : float
            Estimated fraction of Hii regions at latitude coordinate.
        n : float
            Estimated number of Hii regions at coordinate.

        """
        # Convert distance from sun to galactic centre distance:

        rg, z = self.radial_dist_height(l, b, self.r_range)
        Dtot, Dmol, Datom, Dmolh, Datomh = self.total_dust(rg, z)

        centre_dust = self.gal_cent_dust()
        # Get dust at galactic centre (8.122kpc).

        cI, gcV = 0.48, 32
        # Conversion from Av to AI and V band exctinction at galactic centre.

        kI = gcV * cI / centre_dust
        # Normalise galactic dust at galactic centre, for I band extinction at
        # centre.

        AI = [kI * np.trapz(Dtot[0:i], self.r_range[0:i]) for i in \
              range(len(self.r_range))]
        # Integrate dust along line of sight and multiply to get
        # extinction at each point.
        AI = np.array(AI)

        delta_p = 2.512 ** (-AI)
        # Factor by which probablity is supressed.
        # Modelled on delta flux = 2.512**-m

        return AI, delta_p, Dtot, np.max(Dtot)

    def apply_wr_prior(self, mu=3000, sigma=1100, delta=np.ones((15000,)),
                       min_dist=300):
        """
        Applies WR prior from mu and sigma (outputted from the hii_numbers function), to likelihood parallax.
        This produces the posterior and determines the most likely distance for the star.

        Parameters
        -------
        mu : float
            Specify average of prior distribution.
        sigma : float
            Specify sigma of prior distribution.
        min_dist : float
            'Cutoff' distance below which we expect no WR.
        delta : array
            The 'delta' factor supressing the HII region likelihood, obtained from the dust distribution.

        Returns
        -------
        normed : array
            Normalised posterior distribution.
        maxpt : float
            Most likely distance.

        """
        r_pos = self.r_range[self.r_range > min_dist]
        # Get r greater than minimum distance only.
        mu_prior, sigma_prior = mu, sigma
        # unnorm_pos = 1/(2*np.pi*self.err*sigma_prior)*np.exp(-1/2*((self.dpt-1/r_pos)**2/self.err**2+
        # (r_pos-mu_prior)**2/sigma_prior**2))
        # Same as below, just in same form as paper...

        prior1 = self.gauss_simple(r_pos, sigma_prior, mu_prior)
        prior = prior1 * delta[self.r_range > min_dist]
        # Add in dust scaling.
        likelihood = self.gauss_simple(self.dpt, self.err, 1 / r_pos)

        unnorm_pos = prior * likelihood

        # Positive values of r
        r_neg = self.r_range[self.r_range <= min_dist]
        unnorm_neg = 0 * r_neg
        # Negligable values of r (below min dist).

        """Get the prior over the whole r range and normalise"""
        unnorm_prior = np.concatenate((unnorm_neg, prior))
        # Unnormalised prior (for plotting, analysis etc.)
        # Prior is 0 below min_dist.
        unnorm_dust = np.concatenate((unnorm_neg,
                                      delta[self.r_range > min_dist]))
        unnorm_hii = np.concatenate((unnorm_neg, prior1))
        Z_prior = np.trapz(unnorm_prior, self.r_range)
        # Integrated normalisation.
        normed_prior = unnorm_prior / Z_prior
        # Normalised posterior. Area under curve sums to ~one.

        normed_dust = unnorm_dust / np.trapz(unnorm_dust, self.r_range)
        normed_hii = unnorm_hii / np.trapz(unnorm_hii, self.r_range)

        """Get the posterior over the whole r range and normalise"""
        unnorm = np.concatenate((unnorm_neg, unnorm_pos))
        # Unnormalised posterior.
        Z = np.trapz(unnorm, self.r_range)
        # Integrated normalisation.
        normed = unnorm / Z
        # Normalised posterior. Area under curve sums to ~one.

        """Select maximum point"""

        maxpt = self.r_range[np.argmax(normed)]
        # Maximum point from numerical array.

        """Save prior"""

        self.dist = normed
        # Save the distribution in the object!
        self.prior = normed_prior
        self.dust_component = normed_dust
        self.hii_regions = normed_hii
        # Save the prior.

        return normed, maxpt



