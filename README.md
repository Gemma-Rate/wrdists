# wrdists
Python program to calculate the distances to Galactic Wolf-Rayet stars, using *Gaia* DR2 parallaxes and a Bayesian method. The calculation methods (and much of the code) are the same as those used in the publication [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract). 

Our Bayesian prior is based on a combination of HII regions (to simulate massive star locations within the Galaxy) and dust extinction models (to simulate the distribution that can be observed by Gaia). The code can optionally account for underestimations in the DR2 parallax uncertainty and apply a zero point correction. Full posterior distributions can also optionally be viewed (as matplotlib images) or saved (as arrays which may be used, for example, as priors for DR3 parallaxes). 

Further methodology information is available our publication outlined above.

The code can be accessed either through the command line (to calculate distances without modifying the internal processes), or imported into your own python program.

## Installation

Installation uses pip direct from this github page. 
installing into a python virtual environment

## Useage

### Command line

After installation, wrdists can be run directly from the command line. For example, for WR1 (where > is the command line cursor):  
```> wrdists -p 0.285 -pe 0.032 -g 9.79 -ra 10.87 -dec 64.76 -ast 0 -n WR1```

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
```> wrdists -p 1 -pe 2 -g 3 -ra 4 -dec 5 -ast 6 -n 0 -fin \directorypath\test.csv -fout \directorypath\test2.csv``` 

The star names are in the first column and so -n is set to zero (python indexing). The parallax is the second column, and so -p has value 1 and so on. 

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
-ee       Include to disregard resizing of parallax errors (compared to external catalogues, Arenou et al. 2018) 
          and zero point correction. May be useful for data comparison or application to non Gaia parallaxes 
          (e.g Hipparcos) (default = False).
```
For data loading: 
```
-ph       Include if the file input contains a header, to avoid issues with loading data. \
-dmt      Specify a delimiter for the input file (default = ',').
```

### Imported into another program

To access the bayesian distribution class, import as:
```import wrdists.bayesian_functions as bc```

## Dependencies

Requires: 
Python 3.0 or better. 

## Attribution and acknowledgements
If you have found this code useful, then please cite it as [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract).

This work has made use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.
