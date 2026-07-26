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

| File                           | Description                                                                                                  |
|--------------------------------|--------------------------------------------------------------------------------------------------------------|
| run_example.py                 | Execute this Python script to run this example yourself. It uses a Windows-safe main entrypoint.             |
| **Input:**                     |                                                                                                              |
| config.ini                     | This configuration file controls what CurveSimulator does. Click the link to learn about the file structure. |
| flux.csv                       | Flux observations. Click the link to learn about the file structure.                                         |                                                                                     |
| rv.csv                         | Radial velocity observations. Click the link to learn about the file structure.                              |
| tt.csv                         | Tansit times observations. Click the link to learn about the file structure.                                 |                                                                             |
| **Output:**                    | Created by CurveSimulator in this example.                                                                   |
| sim_flux.csv                   | Simulated flux including white noise.                                                                        |
| video.mp4                      | Video of the simulation.                                                                                     |
| videoHQ.mp4                    | Backup video in high quality. Will not be overwritten when this example runs.                                |
| Subdirectories 0000, 0001, ... | Contain result plots and a JSON file with transit parameters, body parameters, fit quality, and more.        |


## Flux Files

Set the parameter `flux_file` in the config file to the path and name of the 
flux file.

### Format
Comma-separated values (CSV)  
Decimal separator is "."

### Columns

| Name     | Unit   | Description                                                                                           |
|----------|--------|-------------------------------------------------------------------------------------------------------|
| time     | [days] | Time of observation, typically expressed as a Barycentric Julian Date (TDB); chronologically ordered. |
| flux     | [1]    | Observed relative flux. The flux median outside of transits should be 1.                              |
| flux_err | [1]    | Standard deviation of this particular observation.                                                    |

### Example
```csv
time,flux,flux_err
2460001.899,1.000222,4.06e-04
2460001.905,1.000458,4.06e-04
2460001.911,0.998324,4.06e-04
```


## Radial Velocity Files

Set the parameter `rv_file` in the config file to the path and name of the 
radial velocity file.

### Format
Comma-separated values (CSV)  
Decimal separator is "."

### Columns

| Name   | Unit   | Description                                                                                           |
|--------|--------|-------------------------------------------------------------------------------------------------------|
| time   | [days] | Time of observation, typically expressed as a Barycentric Julian Date (TDB); chronologically ordered. |
| rv_obs | [m/s]  | Observed radial velocity.                                                                             |
| rv_err | [m/s]  | Standard deviation of this particular observation.                                                    |

### Example
```csv
time,rv_obs,rv_err
2460001.926,100.4,0.5
2460002.547,100.2,0.3
2460003.561,102.8,0.9
```

### Offset and Jitter
Specify the parameters `rv_offset` and `rv_jitter` for the star in the config file.


## Transit Time Files

Set the parameter `tt_file` in the config file to the path and name of the 
transit time file.

### Format
Comma-separated values (CSV)  
Decimal separator is "."

### Columns

| Name     | Unit   | Description                                                                                           |
|----------|--------|-------------------------------------------------------------------------------------------------------|
| eclipser | -      | xxx                                                                                                   |
| tt       | [days] | Time of observation, typically expressed as a Barycentric Julian Date (TDB); chronologically ordered. |
| tt_err   | [days] | Standard deviation of this particular observation.                                                    |
| nr       | [1]    | Index of the primary transits of this planet; starting from 0; chronologically ordered.               |

The transit times file can have many columns more, e.g. T1 and T4 which can 
be used to automatically extract transits from flux data. But this is still 
experimental.

### Example
```csv
eclipser,tt,tt_err,nr
InnerPlanet,2460002.699,0.00458,0
InnerPlanet,2460013.746,0.00481,1
OuterPlanet,2460002.885,0.00320,0
```
