# Transit Time File

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
