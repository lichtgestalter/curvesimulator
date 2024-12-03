# CurveSimulator
[![GitHub Wiki](https://img.shields.io/badge/docs-Wiki-red)](https://github.com/lichtgestalter/curvesimulator/wiki)
[![PyPI version](https://badge.fury.io/py/curvesimulator.svg)](https://badge.fury.io/py/curvesimulator)
[![Python Versions](https://img.shields.io/pypi/pyversions/curvesimulator.svg)](https://pypi.org/project/curvesimulator/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![PyPI Downloads](https://static.pepy.tech/badge/curvesimulator)

CurveSimulator: A Star System and Lightcurve Simulator

CurveSimulator generates a video of the movements and eclipses of celestial bodies and the 
resulting lightcurve.
The video simultanously displays a view of the star system from the top and from the side alongside
the lightcurve of the system's total luminosity over time.

It takes just 2 lines of python code to produce the video.

In a configuration file, specify the physical properties of the stars and planets in your system. 
Also, provide some parameters of the video you want to make.
The meaning of all configuration parameters is documented in CurveSimulator's 
[wiki](https://github.com/lichtgestalter/curvesimulator/wiki) and in an example config file.

CurveSimulator is fast and the videos use little disk space. A video takes about the same time 
to produce as its playing time and uses less than 0.5 MB disc space per minute.

CurveSimulator uses ffmpeg to convert the data into a video. 
Download an executable version of ffmpeg from [ffmpeg.org](https://www.ffmpeg.org/download.html).
Extract the downloaded zip file and (on Windows) add "yourdriveandpath\FFmpeg\bin" to the 
environment variable PATH.

For questions and comments open an issue on [GitHub](https://github.com/lichtgestalter/curvesimulator/issues).

