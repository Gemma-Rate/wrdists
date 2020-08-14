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

Aft
can be run directly from the command line. For example, for WR1:
```> wrdists -p 0.285 -pe 0.032 -g 9.79 -ra 10.87 -dec 64.76 -ast 0 -n WR1```

It's also possible to load in the data from a file. Here, instead of entering the values directly, enter the numbers of the columns you are loading from. For example, from this file:
| WR number | Gaia parallax (mas) | Parallax err (mas) | G (mag) | Gaia RA (deg) | Gaia DEC (deg) | Astrometric excess noise (mas) |
|:---------:|:-------------------:|:------------------:|:-------:|:-------------:|:--------------:|:------------------------------:|
| WR1       | 0.285               | 0.032              | 9.79    |               |                |                                |
| WR2       |                     |                    | 11.00   |               |                |                                |
| WR3       | 0.313               | 0.041              | 10.58   |               |                |                                |
| WR4       | 0.229               | 0.041              | 9.68    |               |                |                                |

#### Required arguments





### Imported into another program

For instance, to access the bayesian distribution class, import as:
```import wrdists.bayesian_functions as bc```

## Example outputs

## Dependencies

## Attribution and acknowledgements
If you have found this code useful, then please cite it as [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract).

This 
uses data from the Gaia DR2 

