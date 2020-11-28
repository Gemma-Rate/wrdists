"""Side functions that are not bound to classes.
   Can be used on any data."""

import astropy.coordinates as co
import numpy as np
import astropy.units as u
import pandas as pd
import wrdists.bayesian_functions as bc
import matplotlib.pyplot as plt

def from_decimal(ra, dec):
    """
    Wrapper function to convert ra and dec from decimal format to hms,
    dms.  

    Parameters
    ----------
    ra : float
        RA ICRS in decimal format. 
    dec : float
        DEC ICRS in decimal format. 
        
    Returns
    ----------
    ra: str
        RA ICRS in hms format.
    dec : str
        DEC ICRS in dms format. 
    """
    gs = co.SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
    # Sky coordinate system to use.
    ra = gs.ra.to_string(decimal=False, unit=u.hourangle, sep=' ',
                         precision=2)
    dec = gs.dec.to_string(decimal=False, unit=u.degree, sep=' ', 
                           precision=1)
    # Convert gaia coordinates back to ra and dec hms, dms format
    # from decimals. 

    if len(ra)<11:
        ra = '0'+ra
    if (dec[0].isnumeric() and len(dec)<10):
        dec = '0'+dec
    elif (not dec[0].isnumeric() and len(dec)<11):
        dec = dec[0]+'0'+dec[1:]
    # Make all coordinates same length. 

    if dec[0].isnumeric():
        dec = '+'+dec
    elif dec[0]=='-':
        dec = dec.replace('-', '$-$')
    # Add + to the front of declination coordinates (for catalogue format).

    return ra, dec


def to_decimal(ra, dec, pres=2):
    """
    Wrapper function to convert ra and dec from hms dms to decimal format.

    Parameters
    ----------
    ra: str or list
        RA ICRS in hms format.
    dec : str or list
        DEC ICRS in dms format. 
        
    Returns
    ----------
    ra : float
        RA ICRS in decimal format. 
    dec : float
        DEC ICRS in decimal format. 
    """
    ra_new = []
    dec_new = []
    cs_list = []

    if isinstance(ra, (str,)):
        ra = [ra]
    if isinstance(dec, (str,)):
        dec = [dec]
    # If string, convert to list. 
    
    for r, d in zip(ra, dec):
        try:
            r, d = str(r, 'utf-8'), str(d, 'utf-8')
            # Python 2 byes conversion to string.
        except:
            pass
        try:
            if unicode(d[0], 'utf-8').isnumeric():
            # Python2 version.
                d = '+'+d
        except NameError:
            # Python 3 version.
            if d[0].isnumeric():
                    d = '+'+d
                    
        cs = co.SkyCoord(r+' '+d, unit=(u.hourangle, u.deg))
        new_r = cs.ra.to_string(decimal=True, sep=' ', precision=pres)
        new_d = cs.dec.to_string(decimal=True, sep=' ', precision=pres)
        ra_new.append(new_r)
        dec_new.append(new_d)

    return ra_new, dec_new


def conv_to_galactic(ra, dec):
    """
    Convert RA and DEC to galactic coordinates. 

    Parameters
    ----------
    ra : float
        Gaia RA ICRS in decimal format. 
    dec : float
        Gaia DEC ICRS in decimal format. 
         
    Returns
    ----------
    l: float
        Galactic reference frame longitude.
    b : str
        Galactic reference frame latitude.
    """

    try:
        c_icrs = co.SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree,
                             frame='icrs')
    except ValueError: 
        if dec[0].isnumeric():
                dec = '+'+dec
        c_icrs = co.SkyCoord(ra+' '+dec, unit=(u.hourangle, u.deg))

    c_galactic = c_icrs.galactic
    l, b = c_galactic.l.value, c_galactic.b.value
    return l, b


def conv_from_galactic(l, b):
    """
    Convert l and b from galactic coordinates to RA and DEC. 

    Parameters
    ----------
    l: float
        Galactic reference frame longitude.
    b : str
        Galactic reference frame latitude.
         
    Returns
    ----------
    ra : float
        Gaia RA ICRS in decimal format. 
    dec : float
        Gaia DEC ICRS in decimal format. 
    """
    c_galactic = co.SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')

    c_icrs = c_galactic.icrs
    ra, dec = c_icrs.ra.value, c_icrs.dec.value
    return ra, dec


def standardize_coords(ra, dec):
    """
    Convert RA and DEC to hms, dms format. 

    Parameters
    ----------
    ra : list
        List of strings containing RA ICRS in decimal format. 
    dec : list
        List of strings containing DEC ICRS in decimal format.  
         
    Returns
    ----------
    ra_proc : list
        List of strings containing RA ICRS in hms format.
    dec_proc : list
        List of strings containing DEC ICRS in dms format. 
    """

    ra_proc, dec_proc = [], []

    for r, d in zip(ra, dec):

        try: 
            r, d = float(r), float(d)
            ra_hms, dec_dms = from_decimal(r, d)
            ra_proc.append(str(ra_hms)), dec_proc.append(str(dec_dms))
        except ValueError:
            ra_proc.append(r), dec_proc.append(d)

    return ra_proc, dec_proc


def create_cs(ra, dec):
    """
    Create skycoordinates from ra and dec.  

    According to the Gaia website 
        The reference epoch for Gaia DR2 is J2015.5 (compared to the
        J2015.0 epoch for Gaia DR1). Positions and proper motions are 
        referred to the ICRS, to which the optical reference frame 
        defined by Gaia DR2 is aligned. The time coordinate for Gaia 
        DR2 results is the barycentric coordinate time (TCB).

    According to astropy:
        The J2000 equatorial is 'within 10s mas' of the ICRS frame 
        used in the SkyCoord definition here. 

    Parameters
    ----------
    ra : list
        List of strings containing RA ICRS in hms or decimal format. 
    dec : list
        List of strings containing DEC ICRS in dms or decimal format.  

    Returns
    -------
    ra : list
        List of skycoord ra. 
    dec : list
        List of skycoord dec.
    cs_list : list
        List of full skycoords. 

    """ 

    ra_new, dec_new, cs_list = [], [], []
    # Make an empty list to fill with ra, dec, icrs coordinates. 

    try: 
        ra.append('b')
        ra.pop()
    except TypeError:
        print('Database not loaded...Run class.load_data and continue')
    except AttributeError:
        print('Ensure input is a list of string coordinates, in the ra, dec J2000 ICRS format')
        raise

    for r, d in zip(ra, dec):

        try:
            if d[0].isnumeric():
                d = '+'+d
            # Ensure conformity to astropy by checking input for numeric
            # first value and adding + if required. 
            cs = co.SkyCoord(r+' '+d, unit=(u.hourangle, u.deg))
            cs_list.append(cs)
            # Get sky coordinates from database. Convert galactic coords
            # to icrs for gaia search. 
            ra_new.append(cs.ra)
            dec_new.append(cs.dec)
        except (ValueError, IndexError, TypeError, AttributeError):
            try:
            # Try with decimal values:
                cs = co.SkyCoord(r, d, unit=(u.deg, u.deg))
                cs_list.append(cs), ra_new.append(cs.ra)
                dec_new.append(cs.dec)
            except (ValueError, IndexError):
                # Except non string results:
                ra_new.append(' ')
                dec_new.append(' ')
                cs_list.append(' ')
                # Ignore results that are not strings. 

    return ra_new, dec_new, cs_list

def run_dist(pars, parserr, phots, ra, dec, ast, name, wdust=True, werr=True, 
             md=300, zpt_data=-0.029, err_sig=0.68, plot_image=False, save_distributions=False):
    """
    Get all distances for a list of objects using bayesian method.

    Parameters
    ----------
    pars : numpy array
        Array of parallaxes from Gaia. 
    parserr : numpy array
        Array of parallax errors from Gaia. 
    phots : numpy array
        G band photometry from Gaia
    ra : list of str or float
        Right asenscion in either hms or decimal format.
    dec : list of str or float
        Declination in either dms or decimal format. 
    ast : numpy array
        Astrometric excess noise.    
    md : int
        Minimum distance for the prior in pc. Default=300pc.
    zpt_data : float or list
        Apply zero point correction to parallaxes. Default=-0.029mas.
    err_sig : float (0-1)
        Specify the region of credible interval. Default=0.68, 1 sigma (e.g alternative 0.95=2 sigma).
    plot_image : bool
        Plot the output image distribution of the prior and posterior.
    wdust : bool
        Include the dust extinction in the prior.
    werr : bool
        Include an increase in the parallax error (used for DR2 data only).
        
    Returns
    ----------
    np.array(dists) : numpy array
        Distances at maximum probability.
    np.array(upper) : numpy array
        Upper 68% interval. 
    np.array(lower) : numpy array
        Lower 68% interval. 
    flags : list
        Flags applied to data. 
    """
       
    max_dist, upper, lower, heights, heights_upper, heights_lower, flags, omega, omega_err = [], [], [], [], [], [], \
                                                                                             [], [], []

    for i in range(len(pars)):
    
        try: 
        # Applying individual different zero points to each WR star.
            len(zpt_data)
            zpt = zpt_data[i]
        except TypeError:
        # Single zero point to apply to all data.
            zpt = zpt_data

        maximum_r, interval, height, height_interval, flagstr, fail, distribution = run_dist_single(pars[i], \
                                                     parserr[i], phots[i],\
                                                     ra[i], dec[i], ast[i], 
                                                     name[i], wdust=wdust, zpt=zpt,
                                                     werr=werr, md=md,
                                                     plot_image=plot_image)
        # Calculate the distance for each individual star in the list.
        max_dist.append(maximum_r)
        upper.append(interval[1]), lower.append(interval[0])
        flags.append(flagstr)

        heights.append(height)
        heights_upper.append(height_interval[1])
        heights_lower.append(height_interval[0])

        if save_distributions:
            fname = save_distributions + '\\' + name[i] + '_posterior.csv'
            save_dist_dict = {'Distance (pc)': distribution.r_range, 'Probability': distribution.dist}
            # Obtain distribution distance range and normalized posterior.
            df = pd.DataFrame(data=save_dist_dict)
            df.to_csv(fname, index=False)

        if werr:
            omega_err.append(distribution.err*1e3)
            # Save updated parallax error if it is calculated.

        omega.append(distribution.dpt*1e3)
        # Always save zero point corrected parallax.

    return np.array(max_dist), np.array(upper), np.array(lower), np.array(heights), np.array(heights_upper), \
           np.array(heights_lower), omega, omega_err, flags


def run_dist_single(pars, parserr, phots, ra, dec, ast, name, r_num=15000, 
                    wdust=True, werr=True, md=300, zpt=-0.029,
                    err_sig=0.68, plot_image=False):
    """
    Get distance for a single object using bayesian method.

    Parameters
    ----------
    pars : float
        Parallax from Gaia (mas).
    parserr : float
        Parallax error from Gaia (mas).
    phots : float
        G band photometry from Gaia (mag).
    ra : str or float
        Right asenscion in either hms or decimal format.
    dec : str or float
        Declination in either dms or decimal format. 
    ast : float
        Astrometric excess noise.
    name : str
        Name of object with distance calculation.      
    r_num : float
        Maximum range of data.
    wdust : boolean
        Include dust in distance estimate.  
    werr : boolean 
        Include inflated Gaia errors. 
    md : int 
        Minimum distance for the prior in pc. Default=300pc.
    zpt : float
        Zero point of the parallax (mas).
    err_sig : float
        Percentage coverage of credible intervals (e.g default 0.68 is one sigma, 0.95 is two sigma etc).
    plot_image : boolean
        Plot distribution of likelihood, prior and posterior, together with most likely distance and credible intervals.
    wdust : bool
        Include the dust extinction in the prior.
    werr : bool
        Include an increase in the parallax error (used for DR2 data only).
        
    Returns
    ----------
    maximum_r : float
        Distance at maximum probability.
    upper : float
        Upper 68% interval. 
    lower : float
        Lower 68% interval. 
    flag : str
        Flags applied to data. 
    dist : distribution object
        The full distribution object. 
    """

    r = np.linspace(0, 15000, r_num)
    # Set up distribution of distances r in pc. 
    
    fail = ' '
    # Empty fail to start (this will be added to at the end).

    try:

        """Create distribution object: """

        dist = bc.Distribution(pars, parserr, r)
        # Create distribution object. 
        if werr == True:
            dist.correct_gaia(phots, uwu=False, zpt=zpt)   
            # Correct WR distance, with inflated Gaia errors.
        elif werr == False:
            dist.dpt = dist.dpt - zpt
            # Apply the zero point.
            dist.dpt, dist.err = dist.dpt*1e-3, \
                                         dist.err*1e-3
            # Apply the zero point and adjust the mas to arcsec, without inflated errors 
            # (which should be used for DR2 only). 

        """Calculate posterior"""

        l, b = conv_to_galactic(ra, dec)
        n, nl, nb, sig, mu = dist.hii_numbers(l, b)
        # Run depending on HII regions along line of sight. 
        av, delta_p, dtot, dmol = dist.dust_model(l, b)
        # Dust model. 

        if wdust:
            normedg, maximum_r = dist.apply_wr_prior(mu=mu, sigma=sig, delta=delta_p,
                                                                min_dist=md)
            # Get unnormalised posterior from gaussian prior. 
        else:
            normedg, maximum_r = dist.apply_wr_prior(mu=mu, sigma=sig, min_dist=md)

        dist.maxr = maximum_r
        # Add max r to the class attributes.

        bign, bigr, n, interval = dist.errs(normedg, err_sig)

        height = maximum_r * np.sin(np.radians(b)) + dist.sun_height
        # Height at max likelihood (including approx. earth height above galactic plane).

        height_68_lower, height_68_upper = interval[0] * np.sin(np.radians(b)) + dist.sun_height, \
                                           interval[1] * np.sin(np.radians(b)) + dist.sun_height
        height_interval = [height_68_lower, height_68_upper]

        # Flags are for results with large parallax error, astrometric 
        # excess noise and negative parallaxes:
        flagstr = gaia_flag(dist.dpt, dist.err, ast)

        if plot_image:
            # Plot prior, likelihood and posterior distributions.

            fig = plt.figure(figsize=(18, 0.6 * 18))

            ax = fig.add_subplot(111)
            # Add subplot to figure.

            ax.plot(r, normedg, '-b', label='Posterior', linewidth=3.0)

            ax.plot(maximum_r, np.max(normedg), color='b', linestyle=' ',
                    marker='*', markersize=25, markeredgewidth=0,
                    label='Most likely distance')

            ax.plot(r, dist.prior, 'k-', label='Total prior', linewidth=3.0)
            ax.plot(r, dist.dust_component, 'r-', label='Dust component', linewidth=3.0)
            ax.plot(r, dist.hii_regions, 'r:', label='H$_\mathrm{{II}}$ region component', linewidth=4.0)

            ax.fill_between(bigr, 0, y2=bign, facecolor='blue', alpha=0.5)
            # Plot on area of 68% errors.

            fsize = 20
            # Set font size.

            ax.set_ylabel('P(r|$\Psi$,$\sigma_\Psi$)', fontsize=fsize)
            ax.set_xlabel('d (pc)', fontsize=fsize)
            ax.tick_params(labelsize=fsize)

            figname = plot_image +'\\'+ name + '_distance_distribution.pdf'

            fig.savefig(figname)

    except: # np.linalg.LinAlgError:
            fail = 'Invalid distance for {}'.format(name)
            # Return failure. 
            maximum_r = np.nan
            interval = [np.nan, np.nan]
            flagstr = ' '
            height = np.nan
            height_interval = [np.nan, np.nan]

    return maximum_r, interval, height, height_interval, flagstr, fail, dist

def process_load(data):
    """
    Load in data from file in numeric format  .  
    
    Parameters
    ----------
    data : pandas series
        Pandas series from dataframe.       
        
    Returns
    ----------
    data_num : numpy array
        Data in numpy float 64 array. 
    """

    data_vals_nans = data.replace(to_replace=' ', value=np.nan).values
    data_num = pd.to_numeric(data_vals_nans)

    return data_num

def err_matrix(upper, lower, dist):
    """
    Get upper and lower errors suitable for plotting with matplotlib errorbar.     
    
    Parameters
    ----------
    upper : numpy array
        Upper error bars. 
    lower : numpy array
        Lower error bars      
    dist : numpy array
        Measured value (e.g distance or height)
        
    Returns
    ----------
    err_array : 2D numpy array
        2D data array with errors. 
    """

    max_68 = upper-dist
    min_68 = dist-lower
    err_array = np.stack((min_68, max_68))

    return err_array


def gaia_flag(par, parerr, ast):
    """
    Flag data according to Gaia parameters.
    
    Parameters
    -------
    par : float
        Gaia modified parallax. 
    parerr : float
        Gaia modified uncertainty.
    ast : float
        Astrometric excess noise. 

    Returns
    -------
    flagstr : string
        Flags
    
    """  
    flagstr = ' '
    if par<0:
        flagstr+='n '
        # Parallax is negative.
    if np.abs(parerr/par)>1:
        flagstr+= 'e '
        # Relative parallax error is large (above 1).
    if ast>=1:
        flagstr+='a '
        # Astrometric excess noise is large.
    if all([ast<1, np.abs(parerr/par)<1, par>0]):
        flagstr+='g '
        # All good parameters.  

    return flagstr
