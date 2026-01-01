# CHANGELOG

### 0.2.15  Text "lichtgestalter/CurveSimulator" added to video.
- Debugging-Code (especially for state vector calculation) moved to backup 
files in folder "debug".
### 0.2.16  Lightcurve in video now contains all data.
### 0.2.17  Calculates for each body how many video frames are needed to complete one orbit.
### 0.2.18  Using Kepler's Third Law, calculates Period P from semi-major axis a and vice versa.
### 0.2.19  Calculates the initial velocity of the primary body based on the principles of
- conservation of momentum and the center of mass motion.
### 0.2.20  Correct handling of secondary stars
### 0.2.21  Configfile errorproofing
### 0.2.22  Example videos
- TIC470710327 (with hacked velocity correction)
- TOI-4504 / TIC349972412 (1 star, 2 transiting planets, 1 non-transiting planet)
### 0.3.0   Limb darkening alternative models/parameters
### 0.3.1   Check if ffmpeg is available
### 0.3.2   Parameter output to JSON-File
- b (Transit Impact parameter): minimum projected separation at the time of transit, in stellar radii
    Bug Gefahr: der relative Radius bezieht sich auf den verdeckten Koerper? Das waere bei sekundaeren Eclipses wohl unerwuenscht
- Transit times:
    - T1, T2, TT, T3, T4, T14
    - TT (Time of Transit): The time of minimum projected separation between the star and planet, as seen by the observer
  - T12 (Ingress) is the period when an exoplanet begins to cross in front of its host star. It consists of two parts:
      - T1 Exterior ingress (first contact): The moment when the planet's disk first touches the edge of the star's disk.
      - T2 Interior ingress (second contact): The point at which the entire planet has moved onto the stellar disk.
  - T23 is called full transit
  - T34 (Egress) is the period when an exoplanet exits its transit across the star. It also has two parts:
      - T3 Interior egress (third contact): The moment when the planet's disk touches the opposite edge of the star's disk as it begins to exit.
      - T4 Exterior egress (fourth contact): The point at which the planet has completely moved off the star's disk, concluding the transit.
  - T14 is called total transit
  - The accuracy of all transit times depends on the parameter dt, the real time difference between simulation iterations.
  - By contrast the timestamps and luminosity values of the lightcurve minima have been slightly adjusted using quadratic interpolation.
### 0.3.3   Improved orbit simulation using Verlet integration. Orbits are now numerically stable.
### 0.3.4   Result files will no longer be overwritten. Instead, a new file with a running number will be created.
### 0.3.5   Result file usage examples (result_evaluation.py)
### 0.3.6   Postponed! state_vector_to_keplerian_elements
### 0.3.7   Result output enhancements
- Angles are now stored in the result file not only in radians but also in degrees.
- Comment-Field in Config-File -> Result file, plots
- Depth of each transit added to result file
- TOI-4504 Research
    - Save plots
    - Sensitivity analysis:  What happens wenn parameters are changed + or - 3 standard deviations?
        - Checked for mass, e, i, P of planets c and d.
        - Checked for Omega and omega of planet c.
        - Result:
            - For some changes of Omega (e.g. +3 or +6 standard deviations) the decline
            - of planet c's transit durations is significantly delayed.
            - In each scenario, planets c and d have alternately longer phases with transits.
     - Manually extracted transit parameters (T1, T2, TT, T3, T4, depth) with high accuracy from TESS data.
### 0.3.8 Artificial Star System
- Under which circumstances are planets alternately transiting?
    - Larger and closer second planets with different inclinations change the inclination of the first planet faster.
    - The effect changes somewhat smoothly with changing with increasing mass and decreasing distance of the second planet.
    - This supports the thesis that the alternately transiting behaviour is real and not just an artefact.
### 0.3.9 Result output enhancements
- More flexible result evaluation
### 0.4.0 Constant energy check
- Added function energy() to check if the total energy in the star system stays constant as it should.
- Results:
    - Precision: Error remains below ±2×10⁻⁷ relative to the mean (1 part in 5 million).
    - Periodicity of error: Repeats about every 5,500 iterations (Δt=600 s per iteration) in TOI-4504
    - This indicates a highly stable system with minimal energy drift over time.
### 0.4.1 Fit simulation to d-Transits
- Planet d's mean anomaly adjusted in the config file until d's transit at 2460736.63 happens at the right time.
- But then the transit before was a whole day off (2460696.46 iso 2460695.53).
- - -> postponed
### 0.4.2 Trying out Rebound and other software
- N-Body programs:
    - Exostriker: does not run correctly on Windows (had to rename python.exe, still had problems)
    - Rebound: Excellent package. I consider using Rebound internally in Curvesim, replacing some of my own functions.
- Programs for estimating planet system parameters:
    - Juliet: Project inactive?
    - exoplanet: Project inactive?
    - Exofast? Seems professionell. Doubtful if it works on Windows .
### 0.4.3 Fitting TOI-4504
### 0.5 Partial Rebound integration
- CurveSimulator now uses Rebound for integration and for calculating results (T1, T2, TT, T3, T4, depth, impact).
### 0.5.2
- Searched with BoxLeastSquared after transits. Confirmed TOI4504b- and c-transits. No evidence for d-transits in sectors 1 - 69 :(
### 0.5.3 MCMC
- Function corresponding_flux() converts TESS-data into data corresponding to a simulated lightcurve.
- Added cs_flux_data.py to CurveSimulator
- First rudimentary integration of emcee (MCMC library) with just one Parameter (Planet-b Period).
- Introduced silent mode (so Curvesimulator does not write into the console during MCMC)
### 0.5.4 Flexible dt and intervals
- Now intervals can be defined, inside of which body positions and flux is recorded.
- Outside of these intervals the integration is calculated with a (potentially very large) default dt.
- This speeds the integration up significantly.
- dt can be set to different values in each time interval.
### 0.5.5   MCMC integration, part 2
- MCMC can now change every parameter of every body (if it shall be fitted)
- CurveSimulator can now save a simulated lightcurve plus white noise as a csv-file.
- This file can be used as "measured" data to check if mcmc can generate meaningful results in this configuration.
- MCMC results and plots are now saved.
### 0.5.6 Fitting TOI-4504
- MCMC uses multiprocessing now.
- Input and output of rebound (body orbit parameters) are now Jacobian iso 
  Cartesian!
### 0.5.7 Automated Flux Data Download
- action get_flux downloads observed flux according to tt_file