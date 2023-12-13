# Package mleetas 

A python package to fit and simulate the Epidemic Type Aftershock model in seismology (Ogata, 1988).

## Overview
The mleetas package contain 3 modules:
    - etas       : A module to fit the classic ETAS model. (5 parameters A,c,p,al,mu)
    - etasi      : A module to fit the ETASI model. A modified version of the classic ETAS to take into account a rate dependent incompletness effect (Hainzl; 2021)
    - simulation : A module to generate synthetics ETAS and ETASI model, background rate, Gutenberg-Richter magnitude distribution and more.

For now refer to the inner function documentations and the two example file in the dir. "example/"

## Package installation
Download the repo

    git clone https://github.com/Luchorigami/mleetas.git

Install the package with pip

    pip install mleetas/

## Example
You can find examples to simulate and fit etas with this package in the example dir.

## References
- Hainzl, S., 2021. ETAS-Approach Accounting for Short-Term Incompleteness of Earthquake Catalogs. Bulletin of the Seismological Society of America. https://doi.org/10.1785/0120210146

- Zhuang, J., Harte, D., Werner, M.J., Hainzl, S., Zhou, S., 2012. Basic models of seismicity: Temporal models. Community Online Resource for Statistical Seismicity Analysis Theme V.

Luc Moutote

