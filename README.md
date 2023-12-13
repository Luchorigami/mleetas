![event](logo_mleetas.png)
# Package MleETAS 

A python package to fit and simulate the Epidemic Type Aftershock model in seismology (Ogata, 1988).

The parameter estimation rely on the maximization of a likelyhood function trough scipy L-BFGS-B optimization.

## Overview
The mleetas package contain 3 modules:

1. **mleetas.etas**:
    * A module to fit the classic ETAS model. 5 parameters: (A,c,p,al,mu)
3. **mleetas.etasi**:
    * A module to fit the ETASI model. A modified version of the classic ETAS to take into account a rate dependent incompletness effect (Hainzl; 2021). 7 parameters: (A,c,p,al,mu,b,Tb)
5. **mleetas.simulation**:
    * A module to generate synthetics ETAS and ETASI catalogs, stationary background catalogs, Gutenberg-Richter magnitude distribution and more.

for more detail, refer to the code documentation inside functions and the 2 example python file in the dir. "example/"

## Package installation
Download the repo wherever you like (Under the directory name mleetas/)

    git clone https://github.com/Luchorigami/mleetas.git

Install the package with pip

    pip install mleetas/

## Example
You can find examples to simulate and fit etas with the MleETAS package in the example dir.

## References
- Hainzl, S., 2021. ETAS-Approach Accounting for Short-Term Incompleteness of Earthquake Catalogs. Bulletin of the Seismological Society of America. https://doi.org/10.1785/0120210146

- Zhuang, J., Harte, D., Werner, M.J., Hainzl, S., Zhou, S., 2012. Basic models of seismicity: Temporal models. Community Online Resource for Statistical Seismicity Analysis Theme V.

Luc Moutote

