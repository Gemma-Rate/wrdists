# wrdists
Python program to calculate the distances to Galactic Wolf-Rayet stars, using Gaia DR2 parallaxes and a Bayesian method. The calculation methods (and much of the code) are the same as those used in the publication [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract). 

Our Bayesian prior is based on a combination of HII regions (to simulate massive star locations within the Galaxy) and dust extinction models (to simulate the distribution that can be observed by Gaia). The code can optionally account for underestimations in the DR2 parallax uncertainty and apply a zero point correction. Full posterior distributions can also optionally be viewed (as matplotlib images) or saved (as arrays which may be used, for example, as priors for DR3 parallaxes). 

Further methodology information is available our publication outlined above.

The code can be accessed either through the command line (to calculate distances without modifying the internal processes), or imported into your own python program.

---

## Installation

Installation requires python and uses pip direct from this github page. 
installing into a python virtual environment

## Useage

### Command line

After installation, wrdists can be run directly from the command line. For example, for WR1 (where > is the command line cursor):
```> wrdists -p 0.285 -pe 0.032 -g 9.79 -ra 10.87 -dec 64.76 -ast 0 -n WR1```

This will produce the following output:




It's also possible to load in the data from a file. Here, instead of entering the values directly, enter the numbers of the columns you are loading from, as well as the path to load from/save to. For example, from this file (test.csv):
| WR number | Gaia parallax (mas) | Parallax err (mas) | G (mag) | Gaia RA (deg) | Gaia DEC (deg) | Astrometric excess noise (mas) |
|:---------:|:-------------------:|:------------------:|:-------:|:-------------:|:--------------:|:------------------------------:|
| WR1       | 0.285               | 0.032              | 9.79    | 10.868        |                |            0                   |
| WR2       |                     |                    | 11.00   | 16.346        |                |           0.27                 |
| WR3       | 0.313               | 0.041              | 10.58   |               |                |                                |
| WR4       | 0.229               | 0.041              | 9.68    |               |                |           0.06                 |

This code will load 
and save the distance data into the file test_2.csv at the path specified.

```> wrdists -p 1 -pe 2 -g 3 -ra 4 -dec 5 -ast 6 -n 0 -fin \directorypath\test.csv -fout directorypath\test2.csv``` 



#### Required arguments

-p Gaia DR2 parallax (mas) or column number.\
-pe Uncertainty of Gaia DR2 parallax (mas) or column number.\
-g Gaia G band magnitude (mag) or column number.\
-ra Gaia Right Ascension (deg) or column number.\
-dec Gaia Declination (dec) or column number.\
-ast Gaia Astrometric excess noise (mas) or column number.\
-n Star name or column number.

Also required when loading data from a file:

-fin File path to load csv containing data.
-fout File path to save new csv with distances. 


#### Optional arguments

-zpt Set the zero point of the parallaxes (mas) (default = -0.029 mas).


### Imported into another program

It is possible to use the model as part of another code. For instance, to access the bayesian distribution class, import as:
```import wrdists.bayesian_functions as bc```

## Dependencies

## Attribution and acknowledgements
If you have found this code useful, then please cite it as [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract).

This work has made use of data from the European Space Agency (ESA) mission Gaia ([https://www.cosmos.esa.int/gaia] (https://www.cosmos.esa.int/gaia)), processed by the Gaia Data Processing and Analysis Consortium (DPAC, [https://www.cosmos.esa.int/web/gaia/dpac/consortium](https://www.cosmos.esa.int/web/gaia/dpac/consortium)). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.
