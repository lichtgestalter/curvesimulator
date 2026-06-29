import os
from setuptools import setup, find_packages

# Get the directory containing this file
current_directory = os.path.abspath(os.path.dirname(__file__))
# remove the src directory from the path
current_directory = os.path.dirname(current_directory)
# Construct the path to the README file
readme_path = os.path.join(current_directory, 'README.md')

with open(readme_path, "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="curvesimulator",
    version="0.6",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
        "colorama",
        "configparser",
        "corner",
        "emcee",
        # "json", commented out because it is a standard lib
        "matplotlib",
        # "multiprocessing", commented out because it is a standard lib
        "numpy",
        # "numpy>=1.5,<2.3",  example showing how to constrain versions
        "rebound"
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
