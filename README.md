# wrdists
Python program to calculate the distances to Galactic Wolf-Rayet stars, using Gaia DR2 parallaxes and a Bayesian method. The calculation methods (and much of the code) are the same as those used in the publication [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract). 

Our Bayesian prior is based on a combination of HII regions (to simulate massive star locations within the Galaxy) and dust extinction models (to simulate the distribution that can be observed by Gaia). The code can optionally account for underestimations in the DR2 parallax uncertainty and apply a zero point correction. Full posterior distributions can also optionally be viewed (as matplotlib images) or saved (as arrays which may be used, for example, as priors for DR3 parallaxes). 

Further methodology information is available our publication outlined above.

The code can be accessed either through the command line (to calculate distances without modifying the internal processes), or imported into your own python program.

## Installation

Installation uses pip direct from this github page:

## Useage

### Command line

Aft
can be run directly from the command line. For example, for WR1:
`> wrdists -p 0.285 -pe 0.032 -g 9.79 -ra 10.87 -dec 64.76 -ast 0 -n WR1`



## Example outputs

## Dependencies

## Attribution 
If you have found this code useful, then please cite it as [Rate & Crowther, 2020 (MNRAS, 493, 1512)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1512R/abstract).

