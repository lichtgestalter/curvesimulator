# Example 2: Just generating a video

## What it does
Running `run_example.py` reads `config.ini`, which instructs the program to:
- run the integration once
- generate `video.mp4`, showing the movements of the bodies, the light curve 
  and the radial velocity curve

The four individual plots of the video can be enabled and disabled in the 
configuration file by setting the parameters `show_left_plot`, 
`show_right_plot`, `show_lc_plot` and `show_rv_plot` to `True` or `False`. 

## Files in this example

| File                           | Description                                                                                           |
|--------------------------------|-------------------------------------------------------------------------------------------------------|
| run_example.py                 | Execute this Python script to run this example yourself. It uses a Windows-safe main entrypoint.     |
| **Input:**                     |                                                                                                       |
| config.ini                     | This configuration file controls what CurveSimulator does.                                            |
| **Output:**                    | Created by CurveSimulator in this example.                                                            |
| video.mp4                      | Video of the simulation.                                                                              |