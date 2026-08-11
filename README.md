# CurveSimulator

[![GitHub Wiki](https://img.shields.io/badge/docs-Wiki-red)](https://github.com/lichtgestalter/curvesimulator/wiki)
[![PyPI version](https://badge.fury.io/py/curvesimulator.svg)](https://badge.fury.io/py/curvesimulator)
[![Python Versions](https://img.shields.io/pypi/pyversions/curvesimulator.svg)](https://pypi.org/project/curvesimulator/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[//]: # (![PyPI Downloads]&#40;https://static.pepy.tech/badge/curvesimulator&#41;)


_CurveSimulator_ is a reliable n-Body Python package for fast and 
  flexible Orbital Parameter Determination (MCMC).   

It also generates highly configurable videos for teaching, public outreach, 
and scientific insight.   

If _CurveSimulator_ could be useful for your next research project feel 
free to reach out to us. 


### Fast and easy to use

* Just run five lines of Python code that never change.
* Everything else is controlled through one user-friendly configuration file.
* The software is continuously developed and improved.


### Robust and Fast Planetary System Parameter Estimation

_CurveSimulator_ determines orbital parameter using n-Body simulations.

Markov Chain Monte Carlo (MCMC) sampler (emcee package) build the parameter 
posteriors and obtain parameter uncertainties.

_CurveSimulator_ can process different types of observations: Flux, 
transit times and radial velocity.

For fast preliminary results, run least-squares fits (LMfit package) on 
observed transit times.

#### Output
* JSON files contain best-fit parameters, per-transit parameters, fit 
  quality, fit diagnostics and much more.
* API: Interact directly with Python objects which contain all relevant 
  input and output.   
* Videos (see below)
* Plots, e.g. parameter histograms, observed vs. computed data, sampling 
  diagnostics 

#### Performance Example
* Fitting 50 parameters using 20000 flux data points and integrating 
  body movements over 7 years
* 5000 iterations of 100 MCMC Walkers on an average Windows PC under 3 hours
* That's 45 simulation runs per second


### Videos of planetary systems

_CurveSimulator_ generates videos of the movements and eclipses of 
celestial bodies. The video can simultaneously display: 
* a view of the star system from above 
* a view from Earth (edge view) 
* an animated plot of the system's total luminosity (lightcurve) 
* an animated plot of the star's radial velocity over time     

Video generation is very fast: One minute playing time takes only seconds to 
generate.   
The videos are highly customizable and use very little disk space - only 0.5 
MB per minute of playing time.   
