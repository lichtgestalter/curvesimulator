This file is in Markdown format. Make sure to view it with an appropriate 
program so it appears properly formatted. 

- Body parameters in the config file may be specified either as a single 
number/expression or as a tuple of four numbers.  
- For single runs, only the single value or the first value of the tuple is 
relevant.  
- When fitting, parameters specified as single numbers are treated as fixed 
values, while parameters specified as tuples of four numbers (separated by 
commas) are treated as fitting parameters.  
  - The first item of the tuple is the initial value.  
  - The second and third items are the lower and upper bounds.  
  - The fourth item is the standard deviation used to generate random, normally 
  distributed starting values around the initial value.  


| File                           | Description                                                                                           |
|--------------------------------|-------------------------------------------------------------------------------------------------------|
| run_example.py                 | Execute this Python script to run this example yourself. It uses a Windows-safe main entrypoint.     |
| **Input:**                     |                                                                                                       |
| config.ini                     | This configuration file controls what CurveSimulator does.                                            |
| flux.csv, rv.csv, tt.csv       | Observations                                                                                          |
| **Output:**                    | Created by CurveSimulator in this example.                                                            |
| sim_flux.csv                   | Simulated flux including white noise.                                                                 |
| video.mp4                      | Video of the simulation.                                                                              |
| Subdirectories 0000, 0001, ... | Contain result plots and a JSON file with transit parameters, body parameters, fit quality, and more. |