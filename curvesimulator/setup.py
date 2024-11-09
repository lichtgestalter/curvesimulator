from setuptools import setup, find_packages

setup(
    name="curvesimulator",
    version="0.1.5",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
        "numpy",
        "matplotlib",
        "configparser",
        # "math",
        # "time",
        # "sys",
    ],
    author="Uli Scheuss",
    description="Curvesimulator calculates the movements and eclipses of celestial bodies and produces a video of this and of the resulting lightcurve.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/lichtgestalter/curvesimulator",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
