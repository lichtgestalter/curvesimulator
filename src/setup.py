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
    version="0.2.14.1",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
        "numpy",
        "matplotlib",
        "configparser",
    ],
    author="Uli Scheuss",
    description="CurveSimulator generates a video of the movements and eclipses of celestial bodies and the resulting lightcurve.",
    long_description=long_description,
    long_description_conte8t_type="text/markdown",
    url="https://github.com/lichtgestalter/curvesimulator",
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
