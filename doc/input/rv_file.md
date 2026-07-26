# Radial Velocity File

Set the parameter `rv_file` in the config file to the path and name of the 
radial velocity file.

### Format
Comma-separated values (CSV)  
Decimal separator is "."

### Columns

| Name   | Unit   | Description                                                                                           |
|--------|--------|-------------------------------------------------------------------------------------------------------|
| time   | [days] | Time of observation, typically expressed as a Barycentric Julian Date (TDB); chronologically ordered. |
| rv | [m/s]  | Observed radial velocity.                                                                             |
| rv_err | [m/s]  | Standard deviation of this particular observation.                                                    |

### Example
```csv
time,rv,rv_err
2460001.926,100.4,0.5
2460002.547,100.2,0.3
2460003.561,102.8,0.9
```

### Offset and Jitter
Specify the parameters `rv_offset` and `rv_jitter` for the star in the config file.

