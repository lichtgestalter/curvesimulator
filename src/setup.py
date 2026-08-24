import os
from setuptools import setup, find_packages

# Read project info from ../README.md and put it on pypi.org
current_directory = os.path.abspath(os.path.dirname(__file__))  # get the directory containing this file
current_directory = os.path.dirname(current_directory)          # remove the src directory from the path
readme_path = os.path.join(current_directory, 'README.md')      # construct the path to the README file

with open(readme_path, "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="curvesimulator",
    version="0.6.3.2",
    packages=find_packages(),
    install_requires=[
        "colorama",             # prints colored text in console
        "configparser",         # reads config file
        "corner",               # generates corner plot for mcmc results
        "emcee",                # Markov Chain Monte Carlo
        "lmfit",                # local minimization fits
        "matplotlib",           # plots
        "numpy",                # numerical (vector) operations
        "pandas",               # tabular flux / TT / RV data
        "rebound>=5.1.1",       # N-body integration
        "scipy",                # statistics and optimization helpers
        "tqdm",                 # emcee progression bar
        # "matplotlib==3.10.0",   # plots
        # "numpy>=1.5,<2.3",      # example showing how to constrain versions
        # "lightkurve",           # TESS light-curve download/processing
    ],
    author="Uli Scheuss",
    description="CurveSimulator is a n-body library for orbital parameter determination and visualization of exoplanet systems.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lichtgestalter/curvesimulator",
    classifiers=[
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
