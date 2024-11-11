# curvesimulator
[![PyPI version](https://badge.fury.io/py/curvesimulator.svg)](https://badge.fury.io/py/curvesimulator)
[![Python Versions](https://img.shields.io/pypi/pyversions/curvesimulator.svg)](https://pypi.org/project/curvesimulator/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<br>
A Star System and Lightcurve Simulator<br>
<br>
Curvesimulator produces a video of the movements and eclipses of celestial bodies and of the resulting lightcurve.<br>
<br>
Curvesimulator is fast and the videos use little disk space. A video takes about the same time to produce as its playing time and uses less than 0.5 MB disc space per minute.<br>
<br>
Specify mass, radius, orbital elements and other properties of some stars and planets in a configuration file.<br>
Then run this code to produce the video:
```python
from curvesimulator import curvesim
parameters, bodies, lightcurve = curvesim("MyConfigFileName.ini")
```
The video shows simultanously a view of the star system from the top and from the side and
the lightcurve of the system's total luminosity over time.<br>
<br>
Usually you do not need to look at or even modify the python code. Instead control the program's
outcome with the config file. The meaning of all program parameters is documented in the example config file (see folder 'configurations').<br>
<br>
Curvesim uses ffmpeg to convert the data into a video. <br> 
Download ffmpeg from https://www.ffmpeg.org/download.html. <br>
Extract the zip file and (on Windows) add "<yourdriveandpath>\FFmpeg\bin" to the environment variable PATH.<br>
<br>
For questions and comments just open an issue on https://github.com/lichtgestalter/curvesim/issues to get my attention :)<br>
