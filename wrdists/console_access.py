"""Argument parser for obtaining WR distances"""

import argparse
import wrdists.collated_functions as cf
import pandas as pd

def main():
    """Incorporate into function to be run using the command line when installed."""

    parser = argparse.ArgumentParser()

    """Set Pandas dataframes to show all elements"""

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)

    """Necessary arguments to calculate distances"""

    parser.add_argument('-p', help='Gaia parallax (mas) when in single mode (float). Column number containing parallax \
                                    list when in file mode (int).', action='store', dest='par', type=float, required=True)
    parser.add_argument('-pe', help='Gaia parallax error (mas) when in single mode (float). Column number containing \
                                     parallax error list when in file mode (int).', action='store', dest='parerr',
                               type=float, required=True)

    parser.add_argument('-g', help='Gaia G band magnitude when in single mode (float). Column number containing G mag list \
                                    when in file mode (int).', action='store', dest='g', type=float, required=True)
    parser.add_argument('-ra', help='Gaia right ascension (RA) when in single mode (float). Column number containing RA \
                                     list when in file mode (int).', action='store', dest='ra', type=float, required=True)
    parser.add_argument('-dec', help='Gaia declination (DEC) when in single mode (float). Column number containing DEC \
                                      list when in file mode (int).', action='store', dest='dec', type=float, required=True)

    parser.add_argument('-ast', help='Gaia astrometric excess noise when in single mode (float). Column number containing \
                                    excess noise list when in file mode (int).', action='store', dest='ast', type=float,
                              required=True)

    parser.add_argument('-n', help='Star name or identifier (e.g WR1) when in single mode (str). Column number containing \
                                    list of names when in file mode (int).', action='store', dest='name', type=str,
                              required=True)

    """Optional arguments"""

    # Load in a list of results:
    parser.add_argument('-fin', help='File path from which to load csv containing parameters when executing for \
                                      lists of stars', action='store', dest='list_mode_fin', type=str, default=False)
    parser.add_argument('-fout', help='File path to store csv ouput to when executing with a file input',
                        action='store', dest='list_mode_fout', type=str, default=False)
    parser.add_argument('-ph', help='Preserve the header if the file input contains one',
                        action='store_true', dest='header', default=False)
    parser.add_argument('-dmt', help='Specify a delimiter for the csv', action='store', dest='delimit', type=str,
                        default=',')

    # Other options:
    parser.add_argument('-zpt', help='Set the zero point of the parallaxes (mas) to an alternative value \
                                     (default = -0.029 mas).',  action='store',  default=-0.029, type=float, dest='zpt')
    parser.add_argument('-md','--minimum_dist', help='Set the minimum distance of the prior (pc), which is useful for \
                                                      constraining the prior.',  action='store', default=300, type=float,
                                                      dest='md')
    parser.add_argument('-es','--error_sigma', help='Set the credible interval coverage range.',  action='store',
                                               default=0.68, type=float, dest='esig')
    # Save plots and/or posterior distribution:
    parser.add_argument('-pt', '--plot', help='Plot the output distributions of the prior, likelihood and posterior, \
                                               along with the credible intervals (uncertainty bounds) and most likely \
                                               distance (default = False). The input string should be the path to save \
                                               the plotted image(s).', action='store', default=False, type=str,
                                               dest='plot_data')
    parser.add_argument('-dist', '--distribution', help='Saves the posterior distance distribution as a numpy array for '
                                                        'which can be loaded and used in another python program. The '
                                                        'input string should be the path to save the distribution data.',
                                                        action='store', default=False, type=str, dest='save_distribution')
    # Exclude dust distribution or parallax resizing:
    parser.add_argument('-ed','--exclude_dust', help='Exclude dust from the prior (use HII regions only), which may be \
                                                      useful to compare the effects of different priors (default = False).',
                                                      action='store_true', default=False, dest='dust_exc')
    parser.add_argument('-ee','--exclude_err', help='Exclude resizing of parallax errors (compared to external catalogues, \
                                                     Arenou et al. 2018) and zero point correction. May be useful for data \
                                                     comparison or application to non Gaia parallaxes (e.g Hipparcos) \
                                                     (default = False)', action='store_true', default=False,
                                                     dest='err_exc')


    """Run the code to get distances"""

    args = parser.parse_args()

    if args.list_mode_fin:

        args.par = int(args.par)
        args.parerr = int(args.parerr)
        args.g = int(args.g)
        args.ra = int(args.ra)
        args.dec = int(args.dec)
        args.ast = int(args.ast)
        args.name = int(args.name)
        # Convert params to integer types.

        if args.header:
            data = pd.read_csv(args.list_mode_fin, delimiter=args.delimit)
        else:
            data = pd.read_csv(args.list_mode_fin, header=None, delimiter=args.delimit)
        # Load in the file.

        pars = data.iloc[:, args.par].values
        parserrs = data.iloc[:, args.parerr].values
        phots = data.iloc[:, args.g].values
        ras = data.iloc[:, args.ra].values
        decs = data.iloc[:, args.dec].values
        asts = data.iloc[:, args.ast].values
        names = data.iloc[:, args.name].values
        # Slice out columns with parameters.

        maxr, upper, lower, heights, heights_upper, heights_lower, omega, omega_err, flags = cf.run_dist(pars, parserrs,
                                                                                              phots, ras, decs, asts, names,
                                                                                              wdust= not args.dust_exc,
                                                                                              werr= not args.err_exc,
                                                                                              md=args.md, zpt=args.zpt,
                                                                                              err_sig=args.esig,
                                                                                              plot_image=args.plot_data,
                                                                                              save_distributions=args.save_distribution)
        # Calculate the distances for all data points in the list.

        data_dict = {'Distance (pc)':maxr, 'Upper distance (pc)':upper,
                     'Lower distance (pc)':lower,
                     'Flags for distance':flags,
                     'Distance from plane (|z|) (pc)':heights,
                     '|z| upper bound (pc)':heights_upper,
                     '|z| lower bound':heights_lower}

        if not args.err_exc:
            data_dict.update({'Omega (zero point corrected parallax) (mas)':omega,
                              'Sigma omega (increased error) (mas)':omega_err})

        df = pd.DataFrame(data=data_dict, index=names)
        # Turn the results into a dataframe.

        if args.list_mode_fout:

            df.to_csv(args.list_mode_fout)

            print('                              ')
            print('##############################')
            print('Distributions for file ' + args.list_mode_fin)
            print('Output saved to path ' + args.list_mode_fout)
            print('##############################')
            print('                              ')

        else:
            print('                              ')
            print('##############################')
            print('Distributions for file ' + args.list_mode_fin)
            print('##############################')
            print(df)
            print('                              ')

    else:
        maximum_r, interval, height, height_interval, flagstr, fail, dist = cf.run_dist_single(args.par, args.parerr,
                                                                                              args.g, args.ra, args.dec,
                                                                                              args.ast, args.name,
                                                                                              r_num=15000,
                                                                                              wdust= not args.dust_exc,
                                                                                              werr= not args.err_exc,
                                                                                              md=args.md, zpt=args.zpt,
                                                                                              err_sig=args.esig,
                                                                                              plot_image=args.plot_data)
        if args.save_distribution:
            fname = args.save_distribution + '\\' + args.name + '_posterior.csv'
            save_dist_dict = {'Distance (pc)': dist.r_range, 'Probability': dist.dist}
            # Obtain distribution distance range and normalized posterior.
            df = pd.DataFrame(data=save_dist_dict)
            df.to_csv(fname, index=False)

        print('                              ')
        print('##############################')
        print('Distribution for '+str(args.name))
        print('##############################')
        print('Distance: {:.0f} pc'.format(maximum_r))
        print('Upper interval: +{:.0f} pc ({:.0f} pc)'.format(interval[1]-maximum_r, interval[1]))
        print('Lower interval: -{:.0f} pc ({:.0f} pc)'.format(maximum_r-interval[0], interval[0]))
        print('Flags for distance: '+flagstr)
        # Warnings for negative parallaxes, large errors and high astrometric excess noise.
        if 'a' in flagstr:
            print('-----Warning: High (>1mas) astrometric excess noise! Input parallax and output distance may be unreliable-----')
        elif 'n' in flagstr:
            print('-----Warning: Negative parallax. Distance dominated by prior-----')
        elif 'e' in flagstr:
            print('-----Warning: Large uncertainty on input parallax. Increased prior influence on distance-----')
        print('                              ')
        print('Distance from plane (|z|): {:.2f} pc'.format(height))
        print('|z| upper bound: {:.2f} pc'.format(height_interval[1]))
        print('|z| lower bound: {:.2f} pc'.format(height_interval[0]))
        print('                              ')
        if not args.err_exc:
            # Print out the new parallax (with zero point added) and error.
            print('Omega (zero point corrected parallax): {:.4f} mas'.format(dist.dpt*1e3))
            print('Sigma omega (increased error): {:.4f} mas'.format(dist.err*1e3))
            print('                              ')
