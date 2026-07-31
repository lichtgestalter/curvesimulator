# Sector Parameters File

Set the parameter `sector_params_file` in the config file to the path and 
name of the sector parameters file.

### Format
Comma-separated values (CSV)  
Decimal separator is "."

### Columns

| Name          | Unit | Description                                                                              |
|---------------|------|------------------------------------------------------------------------------------------|
| sector        | [1]  | An id that is shared by all observations with the same offset and jitter.                |
| offset        | [1]  | Initial value for fitting. flux_corr = flux + offset.                                    |
| offset_low    | [1]  | Lower bound for offset while fitting.                                                    |
| offset_up     | [1]  | Upper bound for offset while fitting.                                                    |
| offset_spread | [1]  | Standard deviation of the normally distributed starting values around the initial value. |
| jitter        | [1]  | Initial value for fitting. flux_total_err = sqrt(flux_err^2 + jitter^2).                 |
| jitter_low    | [1]  | Lower bound for jitter while fitting.                                                    |
| jitter_up     | [1]  | Upper bound for jitter while fitting.                                                    |
| jitter_spread | [1]  | Standard deviation of the normally distributed starting values around the initial value. |

### Example
```csv
sector,offset,offset_low,offset_up,offset_spread,jitter,jitter_low,jitter_high,jitter_spread
0,-0.000088,-0.01,0.01,0.00001,0.001664,0,0.003,0.0001
1,-0.000017,-0.01,0.01,0.00001,0.001215,0,0.003,0.0001
2,-0.000019,-0.01,0.01,0.00001,0.001117,0,0.003,0.0001
```
