# Example 1: Fitting Transit Times with LMFit

## What it does
Running `run_example.py` reads `config.ini`, which instructs the program to:
- read `flux.csv`, `rv.csv`, and `tt.csv`
- run the integration once
- generate `video.mp4`, showing the movements of the bodies, the light curve 
  and the radial velocity curve
- add Gaussian white noise to the simulated flux and save it as `sim_flux.csv`
- calculate all transit events and save the parameters of each individual 
  transit, along with extensive additional information, in a JSON file
- generate numerous plots comparing the observed and computed flux, radial 
  velocity, and transit times
- save the results of each run in a separate numbered subdirectory

The individual components can be disabled in the configuration file by 
commenting out the corresponding file names.   
For example, you can disable the generation of the video or the 
simulated flux file and, instead of using the three data sources, only read 
in the flux.

## Files in this example

| File                           | Description                                                                                           |
|--------------------------------|-------------------------------------------------------------------------------------------------------|
| run_example.py                 | Execute this Python script to run this example yourself. It uses a Windows-safe main entrypoint.      |
| **Input:**                     |                                                                                                       |
| config.ini                     | This configuration file controls what CurveSimulator does.                                            |
| flux.csv, rv.csv, tt.csv       | Observations                                                                                          |
| **Output:**                    | Created by CurveSimulator in this example.                                                            |
| sim_flux.csv                   | Simulated flux including white noise.                                                                 |
| video.mp4                      | Video of the simulation.                                                                              |
| videoHQ.mp4                    | Backup video in high quality. Will not be overwritten when this example runs.                         |
| Subdirectories 0000, 0001, ... | Contain result plots and a JSON file with transit parameters, body parameters, fit quality, and more. |