# curvesimulator
[![GitHub Wiki](https://img.shields.io/badge/docs-Wiki-red)](https://github.com/lichtgestalter/curvesimulator/wiki)
[![PyPI version](https://badge.fury.io/py/curvesimulator.svg)](https://badge.fury.io/py/curvesimulator)
[![Python Versions](https://img.shields.io/pypi/pyversions/curvesimulator.svg)](https://pypi.org/project/curvesimulator/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Star System and Lightcurve Simulator

Curvesimulator produces a video of the movements and eclipses of celestial bodies and of the 
resulting lightcurve.
The video shows simultanously a view of the star system from the top and from the side and
the lightcurve of the system's total luminosity over time.

Check the documentation on the [GitHub Wiki](https://github.com/lichtgestalter/curvesimulator/wiki).

Curvesimulator is fast and the videos use little disk space. A video takes about the same time 
to produce as its playing time and uses less than 0.5 MB disc space per minute.

Run just 2 lines of python code to produce the video.
Specify mass, radius, orbital elements and other properties of some stars and planets in a 
configuration file. The meaning of all configuration parameters 
is documented in the wiki and in an example config file.

Curvesim uses ffmpeg to convert the data into a video. 
Download ffmpeg from https://www.ffmpeg.org/download.html.
Extract the downloaded zip file and (on Windows) add "yourdriveandpath\FFmpeg\bin" to the 
environment variable PATH.

For questions and comments just open an issue on https://github.com/lichtgestalter/curvesim/issues to 
get my attention :)

