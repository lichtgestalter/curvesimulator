# Flux File

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

