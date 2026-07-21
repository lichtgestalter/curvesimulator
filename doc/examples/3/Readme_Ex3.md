# Example 2: Fitting Transit Times with LMFit

## What it does
Running `run_example.py` reads `config.ini`, which instructs the program to:
- read `tt.csv`
- perform a least-squares fit with a single fitting parameter: the mean anomaly of the outer planet
- save the results of each run in a separate numbered subdirectory in 
  `lmfit_best_fits.txt`

## Files in this example
| File                           | Description                                                                                           |
|--------------------------------|-------------------------------------------------------------------------------------------------------|
| run_example.py                 | Execute this Python script to run this example yourself.                                              |
| **Input:**                     |                                                                                                       |
| config.ini                     | This configuration file controls what CurveSimulator does.                                            |
| tt.csv                         | Transit Time Observations                                                                             |
| **Output:**                    | Created by CurveSimulator in this example.                                                            |
| Subdirectories 0000, 0001, ... | Contain result plots and a JSON file with transit parameters, body parameters, fit quality, and more. |


## How to specify fixed values, initial values and bounds for fitting
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
  - If you want (almost) uniformly distributed starting values, just specify 
    a standard deviation several times larger than the interval breadth 
    [lower bound; upper bound]. 
