# wrdists
Python program to calculate the distances to Galactic Wolf-Rayet stars, using *Gaia* DR2 parallaxes and a Bayesian method. The calculation methods (and much of the code) are the same as those used in the publication [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract). 

Our Bayesian prior is based on a combination of HII regions (to simulate massive star locations within the Galaxy) and dust extinction models (to simulate the distribution that can be observed by Gaia). The code can optionally account for underestimations in the DR2 parallax uncertainty and apply a zero point correction. Full posterior distributions can also optionally be viewed (as matplotlib images) or saved (as arrays, for example to use in other plots).

Further methodology information is available our in publication outlined above.

The code can be accessed either through the command line, or imported into your own python program.

## Contents
- [Installation](#installation)
- [Useage](#useage)
  - [Command line](#command-line)
    - [Required Arguments](#required-arguments)
    - [Optional Arguments](#optional-arguments)
  - [Use for DR3 and later releases](#use-for-dr3-and-later-releases)
  - [Imported into another program](#imported-into-another-program)
- [Dependencies](#dependencies)
- [Attribution and acknowledgements](#attribution-and-acknowledgements)
- [Update history](#update-history)

## Installation

To install, clone or download the zip file containing the code into your python or virtual environment Lib folder. In the wrdists directory that is then created (which may be wrdists-master) and contains setup.py, run
```pip install .```
at the command line. This will install the module and any dependencies that may be missing. 

## Useage

### Command line

After installation, wrdists can be run directly from the command line. For example, for WR1 run:  
```wrdists -p 0.285 -pe 0.032 -g 9.79 -ra 10.87 -dec 64.76 -ast 0 -n WR1```

This will produce the following output:
```
############################## 
Distribution for WR1 
############################## 
Distance: 3152 pc 
Upper interval: +463 pc (3615 pc) 
Lower interval: -360 pc (2792 pc) 
Flags for distance:  g 

Distance from plane (|z|): 125.39 pc 
|z| upper bound: 140.75 pc 
|z| lower bound: 113.44 pc 

Omega (zero point corrected parallax): 0.3140 mas 
Sigma omega (increased error): 0.0399 mas 
```

It's also possible to load in the data from a file. Here, instead of entering the values directly, enter the numbers of the columns you are loading from, as well as the path to load from/save to. For example, from this file (test.csv):
| WR number | Gaia parallax (mas) | Parallax err (mas) | G (mag) | Gaia RA (deg) | Gaia DEC (deg) | Astrometric excess noise (mas) |
|:---------:|:-------------------:|:------------------:|:-------:|:-------------:|:--------------:|:------------------------------:|
| WR1       | 0.285               | 0.032              | 9.79    | 10.868        |  64.760        |            0                   |
| WR2       |                     |                    | 11.00   | 16.346        |  60.422        |           0.27                 |
| WR3       | 0.313               | 0.041              | 10.58   | 24.732        |  58.156        |           0.10                 |
| WR4       | 0.229               | 0.041              | 9.68    | 40.300        |  56.730        |           0.06                 |

This code will load data from test.csv and save the distance data into the file test_2.csv at the path specified.  
```wrdists -p 1 -pe 2 -g 3 -ra 4 -dec 5 -ast 6 -n 0 -fin \directorypath\test.csv -fout \directorypath\test2.csv -ph``` 

The star names are in the first column and so -n is set to zero (python indexing). The parallax is the second column, and so -p has value 1 and so on. 
The -ph command is included to avoid  loading in the header.

#### Required arguments
```
-p        Parallax (mas) or number of the column containing parallaxes.
-pe       Uncertainty of parallax (mas) or number of the column containing parallax uncertainties.
-g        Gaia G band magnitude (mag) or number of the containing G magnitudes.
-ra       Right Ascension (deg) or number of the column containing RA.
-dec      Declination (deg) or number of the column containing DEC.
-ast      Gaia Astrometric excess noise (mas) or number of the column containing astrometric excess noises.
-n        Star name or number of the column containing names of the stars.
```
Also required when loading data from a file:
```
-fin      File path to load data.
-fout     File path to save new data with distances. 
```

#### Optional arguments
```
-zpt      Set the zero point of the parallaxes (mas) (default = -0.029 mas).
-md       Set the minimum distance of the prior (pc), below which the probability is zero (default = 300pc). 
-es       Set the credible interval coverage range (default = 0.68). 
-pt       Plot the output distributions of the prior, likelihood and posterior, along with the credible intervals 
          (uncertainty bounds) and most likely distance. The input string should be the path to save the plotted 
          image(s) (default = False). 
-dist     Saves the posterior distance distribution as a csv which can be loaded and used in another python 
          program. The input string should be the path to save the distribution data. (default = False). 
-ed       Include to exclude dust from the prior (use HII regions only), which may be useful to compare the 
          effects of different priors. 
-ee       Include to disregard resizing of parallax errors (compared to external catalogues, Arenou et al. 2018) . 
          Required for Gaia DR3 and later data releases. May also be useful for data comparison and application to non Gaia parallaxes 
          (e.g Hipparcos) and later Gaia data releases after DR2 (default = False).
```
For data loading: 
```
-ph       Include if the file input contains a header, to avoid issues with loading data. 
-dmt      Specify a delimiter for the input file (default = ',').
-zpt_list  Include if loading in a column of zero point data from a file.
```

### Use for DR3 and later releases

The default inputs to this code (e.g the zero point and the application of inflated uncertainties) are valid for DR2 data. To use this code for DR3 and later releases, we can disregard the adjustment of underestimated uncertainties by including -ee  in the command line argument. We should also update the zero point, -zpt.
Changing both these options will not modify the prior or method used to calculate the distances, as these two options only alter the input parallax data and uncertainty. 

The global zero point in Gaia DR3 is (at the time of the [EDR3 release](https://www.cosmos.esa.int/web/gaia/earlydr3)) -0.017 mas. However, for objects with very small parallaxes, a change in the zero point can have a large effect on the resulting distance. As DR3 results will come with a [correction formula](https://www.cosmos.esa.int/web/gaia/earlydr3) for the zero point, it may be worth calculating individual zero points for each star, rather than using the global value. 
 
For example, if the DR3 zero point for WR1 is found to be 0.015 mas, a DR3 calculation for WR1 would look like:   
```wrdists -p 0.285 -pe 0.032 -g 9.79 -ra 10.87 -dec 64.76 -ast 0 -n WR1 -ee -zpt 0.015```   
(though replacing the DR2  parallax, position and other parameters with DR3 data).  The zero point and expanded uncertainties can also be changed if 
loading data from a file. Using the same test.csv file above (and again considering DR3 data), we can either apply the same zero point (e.g -0.017) to all the data:  
```wrdists -p 1 -pe 2 -g 3 -ra 4 -dec 5 -ast 6 -n 0 -fin \directorypath\test.csv -fout \directorypath\test2.csv -ph -zpt -0.017 -ee```

Alternatively, if we wish to apply a different zero point to each star,  we can add a column to the file, like below:
| WR number | Gaia parallax (mas) | Parallax err (mas) | G (mag) | Gaia RA (deg) | Gaia DEC (deg) | Astrometric excess noise (mas) | Zero Point (mas)|
|:---------:|:-------------------:|:------------------:|:-------:|:-------------:|:--------------:|:------------------------------:|:-------------------:|
| WR1       | 0.285               | 0.032              | 9.79    | 10.868        |  64.760        |            0                   | -0.032             |
| WR2       |                     |                    | 11.00   | 16.346        |  60.422        |           0.27                 | -0.029             |
| WR3       | 0.313               | 0.041              | 10.58   | 24.732        |  58.156        |           0.10                 | -0.025             |
| WR4       | 0.229               | 0.041              | 9.68    | 40.300        |  56.730        |           0.06                 | -0.017             |

This data can be loaded and run with similar syntax to above:   
```wrdists -p 1 -pe 2 -g 3 -ra 4 -dec 5 -ast 6 -n 0 -fin \directorypath\test.csv -fout \directorypath\test2.csv -ph -ee -zpt_list 7```   
In this instance, the zero points listed in the final (seventh) column are applied to the corresponding WR stars. If -zpt and -zpt_list are accidentally used together, then the code will use the values in the column specified from -zpt_list.

### Imported into another program

Import as ```wrdists.module```. E.g, to access the bayesian distribution class, use ```import wrdists.bayesian_functions as bc```.

## Dependencies

Python 3.0 or newer. Matplotlib, Pandas, Astropy and Numpy modules.

## Attribution and acknowledgements

If you have found this code useful, then please cite it as [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract).

This work has made use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.

## Update history

28/11/2020: Update to version 1.1 with the following major changes:  

- The run_dist_single function in collated_functions.py has been updated to apply the zero point correction independent to the inflated DR2 parallaxes.
  Additionally, a bug preventing the zero point from being properly incorporated into the calculation has been corrected.  
- The run_dist function in collated_functions.py has been updated to allow for a list of zero points to be used. Additionally, a bug preventing the zero point 
  from being updated has been corrected.  
- The console_access function has been updated to include the -zpt_list argument (allowing a user to apply a list of zero points from a file) and propagate
  this to later functions. The zero point corrected parallax and list of zero points (when used) are now also saved in the output file.  
- This Readme has been updated to include sections on using the code with DR3 (and later release) data and to list the update history.   
- Added further updates to the test modules. These are unittests of the distance calculation code, rather than console_access.py tests and can be run by using 
```python -m unittest``` in the wrdists main directory.
- Removed extraneous functions from collated_functions.py.  
                        
