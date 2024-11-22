# curvesimulator
[![GitHub Wiki](https://img.shields.io/badge/docs-Wiki-blue)](https://github.com/lichtgestalter/curvesimulator/wiki)
[![PyPI version](https://badge.fury.io/py/curvesimulator.svg)](https://badge.fury.io/py/curvesimulator)
[![Python Versions](https://img.shields.io/pypi/pyversions/curvesimulator.svg)](https://pypi.org/project/curvesimulator/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<br>
A Star System and Lightcurve Simulator<br>
<br>
Curvesimulator produces a video of the movements and eclipses of celestial bodies and of the resulting lightcurve.<br>

Check the documentation on the [GitHub Wiki](https://github.com/lichtgestalter/curvesimulator/wiki).
<br>

Curvesimulator is fast and the videos use little disk space. A video takes about the same time to produce as its playing time and uses less than 0.5 MB disc space per minute.<br>
<br>
Specify mass, radius, orbital elements and other properties of some stars and planets in a configuration file. Then run 2 lines of python code to produce the video.

The video shows simultanously a view of the star system from the top and from the side and
the lightcurve of the system's total luminosity over time.<br>
<br>
Control the program's outcome with the config file. The meaning of all program parameters is documented in the wiki and in an example config file.<br>
<br>
Curvesim uses ffmpeg to convert the data into a video. <br> 
Download ffmpeg from https://www.ffmpeg.org/download.html. <br>
Extract the downloaded zip file and (on Windows) add "yourdriveandpath\FFmpeg\bin" to the environment variable PATH.<br>
<br>
For questions and comments just open an issue on https://github.com/lichtgestalter/curvesim/issues to get my attention :)<br>
